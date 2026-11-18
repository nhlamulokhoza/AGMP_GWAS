// AGMP_GWAS documentation: navigation, scroll progress, section rail,
// reveal on scroll, side contents, theme and copy buttons.
(function () {
  "use strict";
  const root = document.documentElement;
  root.classList.add("js");
  try { const t = localStorage.getItem("agmp-gwas-theme"); if (t) root.dataset.theme = t; } catch (e) {}
  const calm = matchMedia("(prefers-reduced-motion: reduce)").matches;
  const page = location.pathname.split("/").pop() || "index.html";

  // ------------------------------------------------------------ top nav
  const items = [
    ["index.html#about", "About"], ["index.html#how", "How it works"], ["index.html#sources", "Sources"],
    ["index.html#columns", "Columns"], ["index.html#outputs", "Outputs"],
    ["guide.html", "Guide"], ["builder.html", "Command builder"], ["technical.html", "Technical"],
  ];
  const nav = document.createElement("header");
  nav.className = "topnav";
  nav.innerHTML = '<div class="in"><a class="brand" href="index.html"><i aria-hidden="true"></i>AGMP_GWAS</a>' +
    '<nav class="links" aria-label="Main">' +
    items.map(([h, t], i) => `<a href="${h}" data-href="${h}"><b>${String(i + 1).padStart(2, "0")}</b>${t}</a>`).join("") +
    '<a href="https://github.com/nhlamulokhoza/AGMP_GWAS" rel="noopener">GitHub</a>' +
    '<button type="button" class="theme" aria-label="Switch light or dark theme">Light / dark</button></nav></div>';
  document.body.prepend(nav);
  nav.querySelector(".theme").addEventListener("click", () => {
    const cur = root.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    root.dataset.theme = cur === "dark" ? "light" : "dark";
    try { localStorage.setItem("agmp-gwas-theme", root.dataset.theme); } catch (e) {}
  });
  const navLinks = [...nav.querySelectorAll("a[data-href]")];
  function navOn(href) {
    navLinks.forEach(a => a.classList.toggle("on", a.dataset.href === href));
    const hit = navLinks.find(a => a.dataset.href === href), box = nav.querySelector(".links");
    if (hit && (hit.offsetLeft < box.scrollLeft || hit.offsetLeft + hit.offsetWidth > box.scrollLeft + box.clientWidth))
      box.scrollTo({left: hit.offsetLeft - 40, behavior: calm ? "auto" : "smooth"});
  }

  // ------------------------------------------- scroll progress and rail
  const bar = document.createElement("div");
  bar.className = "progress"; bar.setAttribute("aria-hidden", "true");
  document.body.append(bar);
  const panels = [...document.querySelectorAll("section.panel[id]")];
  let rail = null;
  if (panels.length > 2) {
    rail = document.createElement("nav");
    rail.className = "rail"; rail.setAttribute("aria-label", "Sections on this page");
    rail.innerHTML = panels.map(p => `<a href="#${p.id}" aria-label="${p.dataset.label || p.id}"><span>${p.dataset.label || p.id}</span></a>`).join("");
    document.body.append(rail);
  }
  // Side contents on long pages (guide, technical)
  const toc = [...document.querySelectorAll(".toc a[href^='#']")];
  const tocTargets = toc.map(a => document.getElementById(a.getAttribute("href").slice(1))).filter(Boolean);

  let ticking = false;
  function onScroll() {
    ticking = false;
    const max = root.scrollHeight - innerHeight;
    bar.style.setProperty("--p", max > 0 ? Math.min(1, scrollY / max) : 0);
    const line = innerHeight * 0.4;
    if (panels.length) {
      let cur = panels[0];
      for (const p of panels) if (p.getBoundingClientRect().top <= line) cur = p;
      if (rail) rail.querySelectorAll("a").forEach(a => a.classList.toggle("on", a.getAttribute("href") === "#" + cur.id));
      navOn(navLinks.some(a => a.dataset.href === page + "#" + cur.id) ? page + "#" + cur.id : "");
    }
    if (tocTargets.length) {
      let cur = tocTargets[0];
      for (const t of tocTargets) if (t.getBoundingClientRect().top <= line) cur = t;
      toc.forEach(a => a.classList.toggle("on", a.getAttribute("href") === "#" + cur.id));
    }
  }
  addEventListener("scroll", () => { if (!ticking) { ticking = true; requestAnimationFrame(onScroll); } }, {passive: true});
  addEventListener("resize", onScroll);
  if (!panels.length) navOn(page);
  onScroll();

  // ---------------------------------------------------- reveal on scroll
  // Each panel's heading and its content items slide in, one after another.
  const groups = [".intro", ".stats > *", ".track > *", ".features > *", ".flow", ".flowcard",
                  ".filterbar", ".collist", ".filepills", ".filehint", ".doc > h2", ".builder fieldset"];
  document.querySelectorAll("section.panel, .doc, .builder").forEach(sec => {
    let i = 0;
    sec.querySelectorAll(groups.join(",")).forEach(el => {
      if (el.classList.contains("reveal")) return;
      el.classList.add("reveal");
      el.style.setProperty("--d", (Math.min(i++, 8) * 0.07) + "s");
    });
  });
  const reveals = [...document.querySelectorAll(".reveal")];
  if (calm || !("IntersectionObserver" in window)) reveals.forEach(el => el.classList.add("in"));
  else {
    const io = new IntersectionObserver(es => es.forEach(e => {
      if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); }
    }), {rootMargin: "0px 0px -8% 0px", threshold: 0.08});
    reveals.forEach(el => io.observe(el));
  }

  // ---------------------------------------------------------- footer, copy
  const foot = document.createElement("footer");
  foot.className = "site";
  foot.textContent = "AGMP_GWAS" + (window.AGMP && window.AGMP.version ? " " + window.AGMP.version : "") +
    " · data from the NHGRI-EBI GWAS Catalog, Ensembl, UniProt, PubMed and Europe PMC";
  document.body.append(foot);
  for (const pre of document.querySelectorAll("pre:not(.nocopy)")) {
    const b = document.createElement("button");
    b.className = "copy"; b.type = "button"; b.textContent = "Copy";
    b.addEventListener("click", async () => {
      const code = pre.querySelector("code") || pre;
      const text = [...code.childNodes].filter(n => !(n.classList && n.classList.contains("c")))
        .map(n => n.textContent).join("").replace(/[ \t]+$/gm, "").trim();
      try { await navigator.clipboard.writeText(text); b.textContent = "Copied"; }
      catch (e) { b.textContent = "Select and copy"; }
      setTimeout(() => (b.textContent = "Copy"), 1500);
    });
    pre.prepend(b);
  }
})();
