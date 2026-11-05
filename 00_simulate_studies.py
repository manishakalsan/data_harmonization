"""Simulate three PBMC scRNA-seq 'studies' that disagree in all the usual ways.

Why simulated: the point is to have KNOWN problems planted, so you can check
the pipeline catches each one. Swap in real public data later (see README).

Planted problems
  metadata  - different column names per study (crosswalk needed)
            - sex as M/F, male/female, Female/Male
            - age as "45y", 52, "40-50", "adult"
            - tissue as PBMC / peripheral blood / blood / "PBMC, frozen" (unmapped)
            - disease as healthy/SLE, control/Lupus, SLE only (study C is all SLE)
            - three different 10x chemistries
            - one donor in study B declared female but genetically male
  cell labels - different author vocabularies and granularity ("T cells" vs "CD4 T")
  genes     - A: symbols, B: versioned Ensembl IDs, C: symbols + aliases
              (CD16, CD20), a previous symbol (C10orf54), a gene missing (GAPDH),
              and one unmappable feature
  counts    - study-specific gene-wise bias + capture efficiency (batch effect)
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 20261104
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

CELL_TYPES = ["CD4 T", "CD8 T", "B", "CD14 Mono", "NK"]
PROPS = [0.35, 0.15, 0.15, 0.25, 0.10]
MARKERS = {
    "CD4 T": ["CD3E", "CD3D", "CD4", "IL7R"],
    "CD8 T": ["CD3E", "CD3D", "CD8A"],
    "B": ["MS4A1", "CD79A"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "CD4"],
    "NK": ["NKG7", "GNLY", "FCGR3A"],
}
ISG = ["ISG15", "IFI27", "MX1"]

# study -> (capture efficiency, chemistry string, donors)
# donor tuple: (true_sex, declared_sex, age_raw, tissue_raw, disease_raw)
STUDIES = {
    "A": dict(eff=0.6, assay="10x v2 3'", donors=[
        ("M", "M", "45y", "PBMC", "healthy"),
        ("F", "F", "38y", "PBMC", "healthy"),
        ("F", "F", "51y", "PBMC", "SLE"),
        ("M", "M", "29y", "PBMC", "SLE")]),
    "B": dict(eff=1.0, assay="Chromium Single Cell 3' v3", donors=[
        ("female", "female", 52, "peripheral blood", "control"),
        ("male", "female", 47, "peripheral blood", "control"),   # planted sex error
        ("female", "female", 35, "peripheral blood", "Lupus"),
        ("male", "male", 60, "peripheral blood", "Lupus")]),
    "C": dict(eff=0.8, assay="10x 5' v1", donors=[
        ("Female", "Female", "40-50", "blood", "SLE"),
        ("Female", "Female", "adult", "blood", "SLE"),
        ("Male", "Male", "61", "PBMC, frozen", "SLE"),          # unmapped tissue
        ("Male", "Male", "33", "blood", "SLE")]),
}
AUTHOR_LABELS = {
    "A": {"CD4 T": "CD4 T", "CD8 T": "CD8 T", "B": "B", "CD14 Mono": "Mono", "NK": "NK"},
    "B": {"CD4 T": "T cells", "CD8 T": "T cells", "B": "B cells",
          "CD14 Mono": "CD14+ Monocytes", "NK": "NK cells"},
    "C": {"CD4 T": "Tcell_c1", "CD8 T": "Tcell_c1", "B": "Bcell", "CD14 Mono": "Monocyte", "NK": "NKcell"},
}
COLS = {  # deliberately inconsistent column names
    "A": dict(sample="SampleID", donor="Donor", sex="Sex", age="Age", tissue="Tissue",
              disease="Condition", assay="Chemistry", cell_label="celltype"),
    "B": dict(sample="sample_name", donor="individual", sex="gender", age="age_years",
              tissue="source", disease="disease_state", assay="protocol", cell_label="cell_type_original"),
    "C": dict(sample="id", donor="patient", sex="sex", age="age", tissue="tissue",
              disease="diagnosis", assay="kit", cell_label="annotation"),
}
N_BACKGROUND = 280
BATCH_SD = 0.5        # sd of study-specific gene-wise log fold change
MARKER_LFC = 3.0      # marker gene log fold change in its cell type
CELLS_PER_SAMPLE = 250


def main() -> None:
    rng = np.random.default_rng(SEED)
    RAW.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(ROOT / "config" / "genes.tsv", sep="\t", dtype=str).fillna("")
    sim = pd.DataFrame({
        "ensembl_id": [f"SIMG{i:06d}" for i in range(1, N_BACKGROUND + 1)],
        "symbol": [f"SIMSYM{i}" for i in range(1, N_BACKGROUND + 1)],
        "aliases": "", "previous_symbols": "", "note": "simulated background gene",
    })
    sim.to_csv(RAW / "simulated_genes.tsv", sep="\t", index=False)
    allg = pd.concat([genes, sim], ignore_index=True)
    sym = allg["symbol"].tolist()
    idx = {s: i for i, s in enumerate(sym)}
    G = len(sym)

    base = rng.normal(-1.0, 1.0, G)                      # log mean expression
    base[idx["ACTB"]] = base[idx["GAPDH"]] = 2.5
    base[idx["PTPRC"]] = 1.5

    truth_rows = []
    for study, spec in STUDIES.items():
        batch_lfc = rng.normal(0, BATCH_SD, G)                # study-specific gene bias
        blocks, obs_rows, sample_rows = [], [], []
        for d, (true_sex, decl_sex, age, tissue, disease) in enumerate(spec["donors"], start=1):
            sid, did = f"{study}_S{d:02d}", f"{study}_D{d:02d}"
            is_female = true_sex.lower().startswith("f")
            is_sle = disease.lower() in {"sle", "lupus"}
            ct = rng.choice(CELL_TYPES, size=CELLS_PER_SAMPLE, p=PROPS)
            mu = np.tile(base + batch_lfc, (CELLS_PER_SAMPLE, 1))
            for c, t in enumerate(CELL_TYPES):
                rows = ct == t
                for g in MARKERS[t]:
                    mu[rows, idx[g]] += MARKER_LFC
            mu[:, idx["XIST"]] += 3.0 if is_female else -4.0
            mu[:, idx["RPS4Y1"]] += -4.0 if is_female else 3.0
            if is_sle:
                for g in ISG:
                    mu[:, idx[g]] += 1.5
            lib = rng.lognormal(0, 0.3, (CELLS_PER_SAMPLE, 1)) * spec["eff"]
            mean = np.exp(mu) * lib
            # negative binomial via gamma-Poisson, dispersion 0.3
            counts = rng.poisson(rng.gamma(1 / 0.3, mean * 0.3))
            blocks.append(counts)
            for k in range(CELLS_PER_SAMPLE):
                obs_rows.append({"cell_id": f"{sid}_c{k:04d}", COLS[study]["sample"]: sid,
                                 COLS[study]["cell_label"]: AUTHOR_LABELS[study][ct[k]]})
                truth_rows.append({"cell_id": f"{sid}_c{k:04d}", "true_cell_type": ct[k]})
            sample_rows.append({COLS[study]["sample"]: sid, COLS[study]["donor"]: did,
                                COLS[study]["sex"]: decl_sex, COLS[study]["age"]: age,
                                COLS[study]["tissue"]: tissue, COLS[study]["disease"]: disease,
                                COLS[study]["assay"]: spec["assay"]})
            truth_rows.append({"cell_id": f"__donor__{did}", "true_cell_type": "",
                               "true_sex": "female" if is_female else "male"})

        X = np.vstack(blocks).astype(np.int32)
        var_names = list(sym)
        keep = np.ones(G, bool)
        if study == "B":                                 # versioned Ensembl IDs
            var_names = [f"{e}.{rng.integers(1, 20)}" for e in allg["ensembl_id"]]
        if study == "C":                                 # aliases, previous symbol, gene missing
            ren = {"FCGR3A": "CD16", "MS4A1": "CD20", "VSIR": "C10orf54"}
            var_names = [ren.get(s, s) for s in sym]
            keep[idx["GAPDH"]] = False
        X, var_names = X[:, keep], [v for v, k in zip(var_names, keep) if k]
        if study == "C":                                 # one feature nobody can map
            X = np.hstack([X, rng.poisson(1, (X.shape[0], 1)).astype(np.int32)])
            var_names.append("UNMAPPED_LOC1")

        obs = pd.DataFrame(obs_rows).set_index("cell_id")
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=var_names))
        out = RAW / f"study_{study}"
        out.mkdir(exist_ok=True)
        adata.write_h5ad(out / "counts.h5ad")
        pd.DataFrame(sample_rows).to_csv(out / "samples.csv", index=False)
        print(f"[sim] study {study}: {adata.n_obs} cells x {adata.n_vars} features")

    pd.DataFrame(truth_rows).to_csv(RAW / "_simulation_truth.csv", index=False)
    print("[sim] ground truth (simulation only, never used by the pipeline): data/raw/_simulation_truth.csv")


if __name__ == "__main__":
    main()
