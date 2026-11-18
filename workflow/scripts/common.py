"""Shared helpers for the AGMP GWAS Catalog pipeline."""
from __future__ import annotations

import email.utils
import io
import logging
import sys
import time
from pathlib import Path

import requests
import yaml
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

USER_AGENT = "agmp-gwas-pipeline/1.0 (+https://agmp.h3abionet.org)"


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )
    return logging.getLogger(name)


def load_yaml(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def http_session(total_retries: int = 6, pool: int = 8) -> requests.Session:
    """Session with exponential back-off on rate limits and server errors."""
    s = requests.Session()
    retry = Retry(
        total=total_retries,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST", "HEAD"}),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=pool, pool_maxsize=pool)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers["User-Agent"] = USER_AGENT
    return s


def release_info(url: str, session: requests.Session | None = None) -> dict:
    """Release date (YYYY-MM-DD) and ETag of a file on the GWAS Catalog FTP."""
    s = session or http_session()
    r = s.head(url, timeout=60, allow_redirects=True)
    r.raise_for_status()
    lm = r.headers.get("Last-Modified")
    if not lm:
        raise RuntimeError(f"No Last-Modified header for {url}")
    dt = email.utils.parsedate_to_datetime(lm)
    return {
        "release": dt.strftime("%Y-%m-%d"),
        "last_modified": lm,
        "etag": r.headers.get("ETag", ""),
        "bytes": int(r.headers.get("Content-Length", 0)),
        "url": url,
    }


def fetch_bytes(url: str, log: logging.Logger, session: requests.Session | None = None) -> bytes:
    """Download a file straight into memory (no file written to disk)."""
    s = session or http_session()
    t0 = time.time()
    with s.get(url, stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        buf = io.BytesIO()
        for chunk in r.iter_content(chunk_size=1 << 20):
            buf.write(chunk)
    data = buf.getvalue()
    log.info("downloaded %s (%.1f MB in %.0fs)", url.rsplit("/", 1)[-1],
             len(data) / 1e6, time.time() - t0)
    return data


def read_previous_ids(paths: list[str], base: Path, log: logging.Logger) -> set[str]:
    """Variant IDs already imported into AGMP (column `id`)."""
    import polars as pl

    ids: set[str] = set()
    for p in paths or []:
        path = Path(p) if Path(p).is_absolute() else base / p
        if not path.exists():
            raise FileNotFoundError(
                f"previous_imports file not found: {path}\n"
                "Check the path given to --previous.")
        if path.suffix == ".parquet":
            col = pl.read_parquet(path, columns=["id"])["id"]
        else:
            col = pl.read_csv(path, columns=["id"], infer_schema=False)["id"]
        before = len(ids)
        ids.update(x for x in col.drop_nulls().to_list() if x)
        log.info("previous import %s: %d ids (+%d)", path.name, len(col), len(ids) - before)
    return ids
