# Omics metadata + data harmonization: a reproducible worked example

Three single-cell PBMC "studies" that disagree on column names, vocabularies,
age formats, chemistries, cell-type granularity and gene identifiers are turned
into one validated, ontology-mapped, integrated AnnData object — with every
decision logged and every unresolvable item sent to a curator queue instead of
being guessed.

**The data are simulated on purpose.** Each problem is planted, so you can check
that the pipeline catches it. The ontologies are real (pinned CL, UBERON, MONDO,
EFO releases), so mapping, ambiguity and branch checks behave exactly as they
would on public data. Section 6 shows how to swap in real datasets.

---

## 1. Run it

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh          # ~1 min after the one-time ~150 MB ontology download
```

Two clean runs produce byte-identical TSV/JSON outputs and identical
counts and Harmony embeddings (seeded RNG, pinned ontologies with SHA-256
lockfile, pinned package versions).

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `01_fetch_ontologies.py` | Downloads pinned releases, checks `data-version`, writes `config/ontologies.lock.json` with checksums; refuses to run if a cached file changed |
| 0 | `00_simulate_studies.py` | Writes three studies to `data/raw/` with planted problems |
| 2a | `02_harmonize_metadata.py --initial` | Metadata pass **without** curator decisions → `curator_queue_initial.tsv` (what a curator is sent) |
| 2b | `02_harmonize_metadata.py` | Same pass **with** `config/curator_decisions.yaml` applied → validated `samples.tsv` |
| 3 | `03_harmonize_data.py` | Gene-ID harmonization, data-level QC, batch-aware HVGs, PCA, Harmony, metrics, plot, `harmonized.h5ad` |

## 2. Layout

```
config/
  ontologies.yaml         pinned ontology URLs + expected versions
  ontologies.lock.json    generated: versions + SHA-256
  target_schema.json      JSON Schema for one harmonized sample record (v1.0.0)
  crosswalks.yaml         per-study source column -> target field
  mapping_memory.yaml     reusable free-text -> ontology decisions + cell-type roll-up level
  curator_decisions.yaml  resolutions for queued items, each with a reason
  genes.tsv               gene table: Ensembl ID, symbol, aliases, previous symbols
scripts/
  obo.py                  ~100-line OBO reader (terms, is_a closure, synonyms, obsolescence)
  00_ … 03_ …             pipeline steps above
data/harmonized/          all outputs (see section 4)
```

## 3. Planted problem → where it is caught

| Planted problem | Caught by | Outcome in this run |
| --- | --- | --- |
| Different column names per study | `crosswalks.yaml`; script stops if a mapped column is missing | All 3 studies mapped to one schema |
| Sex as `M`, `male`, `Female` | `mapping_memory.yaml` → PATO; schema enum | 12/12 mapped |
| Age as `45y`, `52`, `40-50`, `adult` | `parse_age()` → min/max years, raw kept | `adult` kept as text, ages null, queued as COARSE — never turned into a number |
| Tissue `PBMC` / `peripheral blood` | mapping memory → UBERON:0000178 *blood* (PBMC is a suspension *from* blood) | mapped |
| Tissue `blood` | automatic exact-label match within UBERON anatomical entity branch | `auto_exact` |
| Tissue `PBMC, frozen` | no memory, no exact match | **UNMAPPED → queue**; resolved in `curator_decisions.yaml` with a note that "frozen" belongs in a separate preservation field |
| Disease `SLE` | exact synonym of MONDO:0007915, unique in branch | `auto_exact` |
| Disease `Lupus` | exact synonym of **two** MONDO terms (MONDO:0004670 lupus erythematosus, MONDO:0007915 SLE) | **AMBIGUOUS → queue**; this was not planted on purpose — it is how the real MONDO release behaves, and the reason auto-matching needs a curator in the loop |
| Three 10x chemistries in free text | mapping memory → EFO, validated under EFO:0010183 *single cell library construction* | 3 distinct EFO terms |
| Author labels at different granularity (`CD4 T` vs `T cells`) | mapped to CL at the granularity stated, then rolled up via the real `is_a` hierarchy to {T cell, B cell, monocyte, NK cell} | `CD4 T` → CL:0000624 → harmonized T cell; `T cells` stays CL:0000084 (not upgraded) |
| Study B: versioned Ensembl IDs | version stripped, matched on stable ID | 302/302 |
| Study C: `CD16`, `CD20` aliases, `C10orf54` previous symbol | lookup order: Ensembl → symbol → previous symbol → alias; keys claimed by two genes are dropped as ambiguous | 2 alias + 1 previous-symbol rescues |
| Study C: unmappable feature | unmapped, reported, removed | 1 feature dropped |
| Study C: GAPDH not measured | genes intersected; union-with-zeros is **not** used (absent ≠ zero) | GAPDH dropped and listed in provenance |
| Donor B_D02 declared female, genetically male | per-donor XIST vs RPS4Y1 | **SEX MISMATCH flagged; not auto-corrected** |
| Study C contains only SLE donors | study × disease cross-tab, Cramér's V | V = 0.50; study C flagged as a single-disease batch |
| Study-specific gene bias + capture efficiency | Harmony on `study_id`, evaluated before/after | see section 5 |

## 4. Outputs (`data/harmonized/`)

| File | Use |
| --- | --- |
| `samples.tsv` | One row per sample in the target schema, plus labels and a `valid` flag |
| `cell_labels.tsv` | Author label → CL term → harmonized term, with evidence and cell counts |
| `curation_log.tsv` | Every value decision: source column, raw value, normalized value, term, status (`curated_memory`, `auto_exact`, `curated_override`, `parsed_*`) — the provenance trail |
| `curator_queue_initial.tsv` / `curator_queue.tsv` | What needed a human before / after decisions (one coarse age remains) |
| `validation_report.json` | Per-sample schema errors, status counts, ontology versions |
| `gene_mapping.tsv` | Every source feature → stable ID + method |
| `data_qc.json` | Sex check, confounding, integration metrics, interferon score by study × disease |
| `harmonized.h5ad` | `layers['counts']` = primary raw counts; `X` = log-normalized; `obsm` = `X_pca`, `X_harmony`, both UMAPs; `uns['provenance']` = schema + ontology versions, config hashes, package versions, excluded samples, dropped genes |
| `umap_before_after.png` | Before/after figure |

## 5. Results from this run, and how to read them

| Metric (higher is better) | Before (PCA) | After (Harmony) |
| --- | --- | --- |
| Batch ASW per cell type (1 = studies mixed) | 0.927 | 0.986 |
| kNN batch entropy (1 = studies mixed) | 0.622 | 0.955 |
| Cell-type ASW (1 = types separated) | 0.547 | 0.554 |
| kNN cell-type purity | 0.832 | 0.858 |

Batch mixing improved sharply while cell-type separation held. Both halves must
be reported: a method that mixes batches by blurring cell types would also
raise the first two numbers.

What these numbers do **not** tell you:

- **Confounding.** Study C has no healthy donors. Correcting on `study_id`
  treats everything that distinguishes C as batch, so any C-specific disease
  signal is at risk of being removed from the embedding. Harmony changes only
  the embedding, not `X`, so the interferon score in `data_qc.json` is still
  higher in SLE within studies A and B. For disease comparisons, use raw counts
  with study (or donor) as a covariate, aggregated to pseudobulk per donor, and
  rely on studies A and B for the contrast; study C cannot separate disease
  from batch on its own.
- **Simulation is easy mode.** Batch effects here are a simple gene-wise shift.
  Real chemistry differences (3' vs 5', cells vs nuclei) also change which
  transcripts are captured and are not fully correctable.
- **Metric choice.** These are simplified versions of scIB metrics. For real
  work use `scib-metrics` (iLISI, cLISI, kBET, graph connectivity, isolated
  label scores) and compare more than one method.

## 6. Swapping in real data

1. Pull public datasets, e.g. from CELLxGENE Census
   (`cellxgene_census.get_anndata(...)`) for already-standardized data, or
   GEO/ArrayExpress for raw submissions. The raw GEO route is the realistic
   curation exercise: sample sheets in GEO are where messy metadata lives.
2. Write one block per study in `crosswalks.yaml`.
3. Run `02_harmonize_metadata.py --initial`, read `curator_queue_initial.tsv`,
   and add decisions with reasons to `curator_decisions.yaml`.
4. Replace `config/genes.tsv` with a full table built from a pinned GENCODE GTF
   (gene_id, gene_name) plus HGNC aliases and previous symbols. The table here
   covers 22 marker genes only; the 280 `SIMG…` genes are simulated.
5. If studies can be re-quantified from FASTQ with one pipeline and reference
   (e.g. nf-core/scrnaseq), do that first. Uniform reprocessing removes more
   variation than any batch-correction method.

## 7. What a production version would add

- `preservation_method`, `development_stage_ontology_term_id` (HsapDv) and
  `self_reported_ethnicity_ontology_term_id` (HANCESTRO) fields; the queued
  `adult` age would map to a HsapDv stage.
- Schema defined once in **LinkML**, generating the JSON Schema, Python classes
  and docs.
- **OAK (oaklib)** or **pronto** instead of `obo.py`, for full OBO/OWL support
  and relations beyond `is_a`; **text2term** or embedding matching for fuzzy
  suggestions with curator confirmation.
- A scheduled job that re-validates all records when an ontology release
  changes and lists newly obsolete terms with `replaced_by` candidates.
- Tests: one per planted problem above, asserting the expected status.

## 8. Platforms you will hear about in interviews

| Platform / tool | What it is for | Practical notes |
| --- | --- | --- |
| CELLxGENE schema + `cellxgene-schema` validator | Standard for single-cell AnnData submissions | Required `obs` ontology columns; `cellxgene-schema validate file.h5ad` checks a file. Census provides harmonized data via API (TileDB-SOMA) |
| HCA metadata standard + ingest | Atlas-scale, modular JSON Schemas | Spreadsheet templates generated from schemas; strong provenance model |
| LinkML | Modeling language | One YAML model → JSON Schema, SQL, Python, docs |
| OLS / Zooma / OAK / text2term | Ontology lookup and mapping | OLS for search and API, Zooma for curated text-to-term suggestions, OAK for programmatic access |
| scanpy + harmonypy / scvi-tools / scib-metrics | Integration and evaluation (Python) | scVI/scANVI for larger, more complex batches; scib-metrics for benchmarking |
| Seurat + harmony / sva (R) | Same in R | `RunHarmony()` on a Seurat object; `ComBat_seq()` for bulk RNA-seq counts |
| Nextflow nf-core (scrnaseq, rnaseq) | Uniform reprocessing | Containerized, versioned pipelines; record pipeline release and parameters |
| ISA tools, GA4GH (Phenopackets, DUO, DRS) | Multi-omics description, clinical data and access | Relevant when the repository holds controlled human data |
