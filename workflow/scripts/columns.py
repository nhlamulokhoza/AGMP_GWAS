"""The 30 standard AGMP columns: order, meaning, the resources each one comes
from, and every transformation applied on the way.

This file is the single description of the output. The export (column
choice), the launcher (`agmp_gwas columns`), the run report and the
documentation all read it, so they always agree with each other and with the
code. If you change how a column is made, change its entry here too.
"""
from __future__ import annotations

# Resources the pipeline reads. "feeds" is shown in the documentation.
RESOURCES = [
    {"id": "gwas_assoc", "name": "GWAS Catalog associations", "short": "Associations",
     "what": "gwas-catalog-associations_ontology-annotated-full.zip: one row per reported "
             "variant-trait association.",
     "url": "https://www.ebi.ac.uk/gwas/docs/file-downloads"},
    {"id": "gwas_anc", "name": "GWAS Catalog ancestry", "short": "Ancestry",
     "what": "gwas-catalog-ancestry.tsv: one row per study population (each stage of each "
             "study), with ancestry, countries and sample descriptions.",
     "url": "https://www.ebi.ac.uk/gwas/docs/file-downloads"},
    {"id": "ensembl", "name": "Ensembl", "short": "Ensembl",
     "what": "REST API, GRCh38 and GRCh37: variants and genes at a position, gene names and "
             "descriptions. Its variants come from dbSNP.",
     "url": "https://rest.ensembl.org"},
    {"id": "vep", "name": "Ensembl VEP", "short": "VEP",
     "what": "Variant Effect Predictor: the consequence of a variant on each transcript.",
     "url": "https://rest.ensembl.org/documentation/info/vep_id_post"},
    {"id": "uniprot", "name": "UniProt", "short": "UniProt",
     "what": "Reviewed (Swiss-Prot) protein entries: accession and function, linked to "
             "Ensembl gene IDs.",
     "url": "https://rest.uniprot.org"},
    {"id": "regions", "name": "AGMP region table", "short": "Regions",
     "what": "workflow/config/regions.yaml: countries to AGMP regions, country-name aliases, "
             "and an ancestry-category fallback.",
     "url": ""},
    {"id": "pipeline", "name": "AGMP_GWAS", "short": "Pipeline",
     "what": "Fixed values required by the AGMP schema, and the run date.",
     "url": ""},
    {"id": "literature", "name": "PubMed and Europe PMC", "short": "Literature",
     "what": "Abstracts and open-access methods sections. Used only for the curation "
             "queue, never written into the AGMP columns.",
     "url": "https://europepmc.org"},
]

# name, group, meaning, sources {resource: field}, transformations (in order)
COLUMNS = [
    {"name": "id", "group": "Variant",
     "desc": "Variant ID (dbSNP rsID).",
     "sources": {"gwas_assoc": "SNPS, SNP_ID_CURRENT, MERGED, STRONGEST SNP-RISK ALLELE",
                 "ensembl": "/overlap/region (variation)"},
     "steps": ["SNPS is split on ';' and ',' into one record per variant; rsIDs are lower-cased.",
               "A merged rsID (MERGED = 1) is replaced by the current one, rs + SNP_ID_CURRENT.",
               "A synonym pair such as rs139025433/rs58533488 becomes the current rsID.",
               "A positional ID (chr7:107194965) becomes an rsID only if, on the paper's genome "
               "build, exactly one rsID starts at that position and carries the reported allele.",
               "Interaction terms (rs1 x rs2), HLA alleles and other non-variant IDs are removed "
               "and listed in unmapped_variants.tsv."]},
    {"name": "PUBMEDID", "group": "Publication", "desc": "PubMed ID of the publication.",
     "sources": {"gwas_assoc": "PUBMEDID"}, "steps": ["Unchanged."]},
    {"name": "phenotype", "group": "Association",
     "desc": "Trait or disease, as reported by the study.",
     "sources": {"gwas_assoc": "DISEASE/TRAIT"}, "steps": ["Unchanged."]},
    {"name": "origin_of_participants", "group": "Population",
     "desc": "Country or countries where participants were recruited.",
     "sources": {"gwas_anc": "COUNTRY OF RECRUITMENT"},
     "steps": ["Unchanged (a comma-separated list, as in the Catalog).",
               "Also decides inclusion: recruitment in an African country brings a study in."]},
    {"name": "geographical_region", "group": "Population",
     "desc": "AGMP region(s) of the study population, e.g. West Africa or "
             "African-American / Afro-Caribbean.",
     "sources": {"gwas_anc": "COUNTRY OF RECRUITMENT, BROAD ANCESTRAL CATEGORY",
                 "regions": "regions, aliases, ancestry_fallback"},
     "steps": ["Each recruitment country is matched to its AGMP region (case-insensitive, after "
               "aliases such as 'United Republic of Tanzania' to 'Tanzania').",
               "All regions found are joined in alphabetical order, e.g. 'East Africa, West Africa'.",
               "If no country maps, the broad ancestral category is used instead, e.g. "
               "'African American or Afro-Caribbean' gives 'African-American / Afro-Caribbean'."]},
    {"name": "mixed_population", "group": "Population",
     "desc": "True when participants come from more than one country.",
     "sources": {"gwas_anc": "COUNTRY OF RECRUITMENT"},
     "steps": ["'True' if the recruitment countries contain a comma, otherwise 'False'."]},
    {"name": "p-value", "group": "Association", "desc": "Association p-value.",
     "sources": {"gwas_assoc": "P-VALUE"},
     "steps": ["Converted to a number.",
               "When records are identical except for the p-value, the smallest is kept."]},
    {"name": "curated_gene_symbol", "group": "Gene",
     "desc": "Gene(s) the variant is mapped to: genes it falls in ('A, B'), or the nearest "
             "gene on each side ('A - B').",
     "sources": {"gwas_assoc": "MAPPED_GENE", "ensembl": "/overlap/region (gene)"},
     "steps": ["Taken from MAPPED_GENE.",
               "Converted positional variants have no Catalog mapping, so it is built the same "
               "way: the named genes overlapping the position ('A, B'), otherwise the nearest "
               "named gene on each side within 500 kb ('UP - DOWN')."]},
    {"name": "study_type", "group": "Fixed", "desc": "Always 'GWAS'.",
     "sources": {"pipeline": "fixed value"}, "steps": ["Set to 'GWAS'."]},
    {"name": "gene_symbol_mined", "group": "Gene", "desc": "Gene(s) named by the authors.",
     "sources": {"gwas_assoc": "REPORTED GENE(S)"}, "steps": ["Unchanged."]},
    {"name": "gene_name", "group": "Gene",
     "desc": "Ensembl gene ID (ENSG...) linked to the variant; one record per gene.",
     "sources": {"gwas_assoc": "SNP_GENE_IDS", "ensembl": "/overlap/region (gene)"},
     "steps": ["Taken from SNP_GENE_IDS; ';' separators become ','.",
               "Converted positional variants get the genes overlapping the position, named "
               "genes in preference to unnamed 'novel transcript' genes.",
               "Split into one record per gene.",
               "Placeholder values (nan, null, none) are cleared."]},
    {"name": "chromosome", "group": "Variant", "desc": "Chromosome.",
     "sources": {"gwas_assoc": "CHR_ID"},
     "steps": ["Taken from CHR_ID.",
               "For converted positional variants without one, read from the position."]},
    {"name": "uniprot", "group": "Gene", "desc": "UniProt accession of the gene's protein.",
     "sources": {"uniprot": "Swiss-Prot Ensembl cross-references, ID mapping, search"},
     "steps": ["Looked up for gene_name: Swiss-Prot entries that cross-reference the Ensembl "
               "gene, then UniProt ID mapping, then UniProt search.",
               "'unresolved' for intergenic records, and when a function was found but no "
               "accession."]},
    {"name": "function", "group": "Gene",
     "desc": "What the gene's protein does. Describes the gene, not the variant.",
     "sources": {"uniprot": "FUNCTION comment", "ensembl": "gene description"},
     "steps": ["UniProt's FUNCTION text, cut to its first 150 characters.",
               "If UniProt has none, the Ensembl gene description.",
               "'unresolved' for intergenic records, and when an accession was found but no "
               "function."]},
    {"name": "data_ac", "group": "Fixed", "desc": "Date this table was made (DDMMYYYY).",
     "sources": {"pipeline": "run date"}, "steps": ["The date of the run, as DDMMYYYY."]},
    {"name": "publication_type", "group": "Fixed", "desc": "Always 'PubMed'.",
     "sources": {"pipeline": "fixed value"}, "steps": ["Set to 'PubMed'."]},
    {"name": "title", "group": "Publication", "desc": "Publication title.",
     "sources": {"gwas_assoc": "STUDY"}, "steps": ["Unchanged."]},
    {"name": "publication", "group": "Publication",
     "desc": "PubMed ID, repeated for the AGMP schema.",
     "sources": {"gwas_assoc": "PUBMEDID"}, "steps": ["A copy of PUBMEDID."]},
    {"name": "publication_year", "group": "Publication", "desc": "Year of publication.",
     "sources": {"gwas_assoc": "DATE"},
     "steps": ["The year part of DATE (the text before the first '-')."]},
    {"name": "source", "group": "Fixed", "desc": "Always 'GWAS Catalog'.",
     "sources": {"pipeline": "fixed value"}, "steps": ["Set to 'GWAS Catalog'."]},
    {"name": "id_in_source", "group": "Publication",
     "desc": "GWAS Catalog study accession (GCST...).",
     "sources": {"gwas_assoc": "STUDY ACCESSION", "gwas_anc": "STUDY ACCESSION"},
     "steps": ["Unchanged. Associations and ancestry are joined on it."]},
    {"name": "variant_type", "group": "Fixed", "desc": "Always 'SNV'.",
     "sources": {"pipeline": "fixed value"}, "steps": ["Set to 'SNV'."]},
    {"name": "intergenic", "group": "Variant",
     "desc": "'intergenic' (outside any gene) or 'non-intergenic'.",
     "sources": {"gwas_assoc": "INTERGENIC", "ensembl": "/overlap/region (gene)"},
     "steps": ["INTERGENIC 1 becomes 'intergenic'; 0 or empty becomes 'non-intergenic'.",
               "Converted positional variants are intergenic when Ensembl finds no gene at the "
               "position."]},
    {"name": "gene_symbol_external", "group": "Gene",
     "desc": "Gene symbol for gene_name.",
     "sources": {"ensembl": "/lookup/id display_name", "gwas_assoc": "MAPPED_GENE, REPORTED GENE(S)"},
     "steps": ["Ensembl's name for the gene in gene_name (GRCh38, or GRCh37 for retired IDs).",
               "If there is none, curated_gene_symbol; failing that, gene_symbol_mined."]},
    {"name": "ethnicity", "group": "Population",
     "desc": "Additional ancestry or ethnicity description, when the study gives one.",
     "sources": {"gwas_anc": "ADDITIONAL ANCESTRY DESCRIPTION"},
     "steps": ["Unchanged. Never inferred: empty when the study gives none."]},
    {"name": "replication_sample_size", "group": "Population",
     "desc": "Replication sample size.",
     "sources": {"gwas_assoc": "REPLICATION SAMPLE SIZE"},
     "steps": ["Thousands separators removed, then the first number kept: "
               "'1,234 African American cases, 2,000 controls' gives 1234."]},
    {"name": "initial_sample_size", "group": "Population", "desc": "Discovery sample size.",
     "sources": {"gwas_assoc": "INITIAL SAMPLE SIZE"},
     "steps": ["Thousands separators removed, then the first number kept, as for "
               "replication_sample_size."]},
    {"name": "broad_ancestral_category", "group": "Population",
     "desc": "GWAS Catalog broad ancestral category, e.g. 'African American or "
             "Afro-Caribbean' or 'NR'.",
     "sources": {"gwas_anc": "BROAD ANCESTRAL CATEGORY"},
     "steps": ["Unchanged.",
               "Also decides inclusion: African or Greater Middle Eastern categories are "
               "explicit; NR and Other are kept as unresolved, for curation."]},
    {"name": "country_of_origin", "group": "Population",
     "desc": "Country of origin of the participants' ancestors.",
     "sources": {"gwas_anc": "COUNTRY OF ORIGIN"}, "steps": ["Unchanged."]},
    {"name": "context", "group": "Variant",
     "desc": "Variant consequence, e.g. intron_variant or missense_variant.",
     "sources": {"gwas_assoc": "CONTEXT", "vep": "most_severe_consequence"},
     "steps": ["Taken from CONTEXT.",
               "Where that is empty (all converted positional variants), Ensembl VEP's most "
               "severe consequence over all transcripts on the primary GRCh38 assembly.",
               "As in the Catalog, upstream_gene_variant and downstream_gene_variant become "
               "intergenic_variant."]},
]

# Applied to whole records rather than one column
RECORD_STEPS = [
    "Only study populations relevant to Africa and the diaspora are kept (four inclusion "
    "routes: explicit ancestry, NR/Other label, African recruitment country, African sample "
    "description).",
    "Associations are joined to those populations by study accession.",
    "One record per variant, per gene, per study population.",
    "Identical records are reduced to one, keeping the smallest p-value.",
]

NAMES = [c["name"] for c in COLUMNS]
