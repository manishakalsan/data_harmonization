# Omics metadata + data harmonization: a reproducible worked example



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

## 3. Outputs (`data/harmonized/`)

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

