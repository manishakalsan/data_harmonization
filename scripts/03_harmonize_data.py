"""Data harmonization: gene IDs -> combine -> data-level QC -> integrate -> evaluate.

Inputs : data/harmonized/samples.tsv, cell_labels.tsv (from step 02), raw counts
Outputs: data/harmonized/
           harmonized.h5ad     raw counts in layers['counts'] (primary),
                               log-normalized X, X_pca, X_harmony (derived),
                               provenance in .uns
           gene_mapping.tsv    every source feature -> stable ID, with method
           data_qc.json        sex check, confounding, integration metrics
           umap_before_after.png
"""
import hashlib
import json
import re
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import harmonypy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_samples
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
RAW, OUT, CFG = ROOT / "data" / "raw", ROOT / "data" / "harmonized", ROOT / "config"
SEED = 20261104
N_HVG, N_PCS = 60, 20
ENSEMBL_VERSIONED = re.compile(r"^((?:ENSG|SIMG)\d+)\.\d+$")


# --------------------------------------------------------------------------- genes
def build_gene_index(genes: pd.DataFrame):
    """key -> stable ID, per method. A key claimed by two IDs is ambiguous and never used."""
    index = {m: {} for m in ["ensembl", "symbol", "alias", "previous_symbol"]}
    ambiguous = set()
    for _, g in genes.iterrows():
        entries = [("ensembl", g.ensembl_id), ("symbol", g.symbol)]
        entries += [("alias", a) for a in str(g.aliases).split("|") if a]
        entries += [("previous_symbol", p) for p in str(g.previous_symbols).split("|") if p]
        for method, key in entries:
            k = key.upper()
            if k in index[method] and index[method][k] != g.ensembl_id:
                ambiguous.add((method, k))
            index[method][k] = g.ensembl_id
    for method, k in ambiguous:
        index[method].pop(k, None)
    return index, ambiguous


def map_feature(name: str, index) -> tuple[str | None, str]:
    m = ENSEMBL_VERSIONED.match(name)
    key = (m.group(1) if m else name).upper()
    for method in ["ensembl", "symbol", "previous_symbol", "alias"]:   # most to least trustworthy
        if key in index[method]:
            return index[method][key], method + ("_version_stripped" if m and method == "ensembl" else "")
    return None, "unmapped"


# --------------------------------------------------------------------------- metrics
def integration_metrics(emb, batch, label, k=30):
    """Two questions, both must be answered: did batches mix, and did biology survive?"""
    # batch ASW per cell type (scIB-style): 1 - |silhouette of batch labels|, higher = better mixed
    basw = []
    for ct in np.unique(label):
        sel = label == ct
        if len(np.unique(batch[sel])) > 1 and sel.sum() > 10:
            basw.append(np.mean(1 - np.abs(silhouette_samples(emb[sel], batch[sel]))))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    nbr = nn.kneighbors(emb, return_distance=False)[:, 1:]
    nb = len(np.unique(batch))
    ent = []
    for i in range(len(batch)):
        _, c = np.unique(batch[nbr[i]], return_counts=True)
        p = c / c.sum()
        ent.append(-(p * np.log(p)).sum() / np.log(nb))
    purity = (label[nbr] == label[:, None]).mean()
    return {
        "batch_asw (1=mixed)": round(float(np.mean(basw)), 3),
        "knn_batch_entropy (1=mixed)": round(float(np.mean(ent)), 3),
        "celltype_asw (1=separated)": round(float((np.mean(silhouette_samples(emb, label)) + 1) / 2), 3),
        "knn_celltype_purity (1=pure)": round(float(purity), 3),
    }


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


# --------------------------------------------------------------------------- main
def main() -> None:
    samples = pd.read_csv(OUT / "samples.tsv", sep="\t")
    labels = pd.read_csv(OUT / "cell_labels.tsv", sep="\t")
    valid = samples[samples.valid]
    rejected = samples.loc[~samples.valid, "sample_id"].tolist()
    print(f"[data] using {len(valid)} valid samples; excluded until curated: {rejected or 'none'}")

    genes = pd.concat([pd.read_csv(CFG / "genes.tsv", sep="\t", dtype=str),
                       pd.read_csv(RAW / "simulated_genes.tsv", sep="\t", dtype=str)]).fillna("")
    index, ambiguous = build_gene_index(genes)
    symbol_of = dict(zip(genes.ensembl_id, genes.symbol))

    import yaml
    xw = yaml.safe_load((CFG / "crosswalks.yaml").read_text())["studies"]

    parts, gmap_rows = [], []
    for study, spec in xw.items():
        a = ad.read_h5ad(RAW / f"study_{study}" / "counts.h5ad")
        sc_, lc = spec["cells"]["sample_col"], spec["cells"]["label_col"]
        a = a[a.obs[sc_].isin(valid.sample_id)].copy()

        # map features -> stable IDs
        mapped = [map_feature(v, index) for v in a.var_names]
        for v, (gid, method) in zip(a.var_names, mapped):
            gmap_rows.append(dict(study=study, source_feature=v, ensembl_id=gid, method=method))
        ids = pd.Series([g for g, _ in mapped], index=a.var_names)
        a = a[:, ids.notna().values].copy()
        ids = ids.dropna()
        # two source features -> one ID: sum counts (and report), never keep both
        dup = ids[ids.duplicated(keep=False)]
        if len(dup):
            print(f"[genes] study {study}: {dup.nunique()} ID(s) hit by >1 feature, counts summed")
        X = pd.DataFrame(a.X.toarray() if hasattr(a.X, "toarray") else a.X,
                         index=a.obs_names, columns=ids.values).T.groupby(level=0).sum().T
        b = ad.AnnData(X=X.values.astype(np.float32), obs=pd.DataFrame(index=X.index),
                       var=pd.DataFrame(index=X.columns))
        b.obs["sample_id"] = a.obs[sc_].values
        lab = labels[labels.study == study].set_index("author_label")
        b.obs["author_cell_type"] = a.obs[lc].astype(str).values
        for col in ["cl_id", "cl_label", "harmonized_id", "harmonized_label"]:
            b.obs[col] = b.obs["author_cell_type"].map(lab[col]).values
        parts.append(b)

    gmap = pd.DataFrame(gmap_rows)
    gmap.to_csv(OUT / "gene_mapping.tsv", sep="\t", index=False)
    print("[genes] mapping methods:\n" + gmap.groupby(["study", "method"]).size().to_string())

    shared = sorted(set.intersection(*(set(p.var_names) for p in parts)))
    dropped = sorted(set.union(*(set(p.var_names) for p in parts)) - set(shared))
    print(f"[genes] {len(shared)} genes shared by all studies; dropped (not measured everywhere): "
          f"{[symbol_of.get(g, g) for g in dropped]}")
    adata = ad.concat([p[:, shared] for p in parts], join="inner", merge="same")
    adata.var["symbol"] = [symbol_of[g] for g in adata.var_names]
    meta_cols = ["sample_id", "study_id", "donor_id", "sex_ontology_term_id", "sex_label",
                 "tissue_ontology_term_id", "disease_ontology_term_id", "disease_label",
                 "assay_ontology_term_id", "assay_label", "age_min_years", "age_max_years", "age_raw"]
    adata.obs = adata.obs.merge(valid[meta_cols], on="sample_id", how="left").set_index(adata.obs_names)
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str).astype("category")
    adata.layers["counts"] = adata.X.copy()

    # ----------------------------------------------------------- normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ----------------------------------------------------------- data-level QC
    qc = {}
    sym2id = {s: g for g, s in zip(adata.var_names, adata.var.symbol)}
    xist = np.asarray(adata[:, sym2id["XIST"]].X).ravel()
    ry = np.asarray(adata[:, sym2id["RPS4Y1"]].X).ravel()
    donor_df = pd.DataFrame({"donor": adata.obs.donor_id.values, "xist": xist, "rps4y1": ry,
                             "declared": adata.obs.sex_label.values}).groupby("donor", observed=True).agg(
        xist=("xist", "mean"), rps4y1=("rps4y1", "mean"), declared=("declared", "first"))
    donor_df["inferred"] = np.where(donor_df.xist > donor_df.rps4y1, "female", "male")
    donor_df["mismatch"] = donor_df.inferred != donor_df.declared
    qc["sex_check"] = donor_df.round(2).reset_index().to_dict(orient="records")
    for d, r in donor_df[donor_df.mismatch].iterrows():
        print(f"[qc] SEX MISMATCH donor {d}: declared {r.declared}, expression says {r.inferred} "
              f"(XIST {r.xist:.2f}, RPS4Y1 {r.rps4y1:.2f}) -> query submitter; NOT auto-corrected")

    ct = pd.crosstab(valid.study_id, valid.disease_label)
    chi2 = chi2_contingency(ct)[0]
    cramers_v = float(np.sqrt(chi2 / (ct.values.sum() * (min(ct.shape) - 1))))
    single = [s for s in ct.index if (ct.loc[s] > 0).sum() == 1]
    qc["confounding"] = {"samples_study_x_disease": ct.to_dict(), "cramers_v": round(cramers_v, 3),
                         "studies_with_one_disease_level": single}
    print(f"[qc] study x disease (samples):\n{ct.to_string()}")
    print(f"[qc] Cramér's V = {cramers_v:.2f}; studies with a single disease level: {single}. "
          "Integration on 'study' can absorb disease signal for these.")

    # ----------------------------------------------------------- integration
    # HVGs chosen per batch and combined, so one study's bias doesn't pick the genes
    adata.layers["lognorm"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key="study_id", flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, mask_var="highly_variable", random_state=SEED)
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["study_id"],
                               random_state=SEED, verbose=False)
    Z = np.asarray(ho.Z_corr)
    adata.obsm["X_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
    adata.X = adata.layers.pop("lognorm")             # X stays log-normalized, not scaled

    batch = adata.obs.study_id.to_numpy(dtype=str)
    lab = adata.obs.harmonized_label.to_numpy(dtype=str)
    qc["integration"] = {"before (PCA)": integration_metrics(adata.obsm["X_pca"], batch, lab),
                         "after (Harmony)": integration_metrics(adata.obsm["X_harmony"], batch, lab)}
    print("[int] metrics:\n" + pd.DataFrame(qc["integration"]).to_string())

    # disease signal: does the SLE interferon difference survive integration?
    isg = [sym2id[g] for g in ["ISG15", "IFI27", "MX1"]]
    score = np.asarray(adata[:, isg].X).mean(1).ravel()
    adata.obs["isg_score"] = score
    qc["isg_score_by_study_disease"] = (pd.DataFrame({"study": batch, "disease": adata.obs.disease_label.to_numpy(dtype=str),
                                                      "isg": score.astype(float)}).groupby(["study", "disease"], observed=True)
                                        .isg.mean().round(3).reset_index().to_dict(orient="records"))

    # ----------------------------------------------------------- UMAPs
    emb = {}
    for rep in ["X_pca", "X_harmony"]:
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=SEED)
        sc.tl.umap(adata, random_state=SEED)
        emb[rep] = adata.obsm["X_umap"].copy()
    adata.obsm["X_umap_pca"], adata.obsm["X_umap_harmony"] = emb["X_pca"], emb["X_harmony"]
    del adata.obsm["X_umap"]
    plot(adata, emb, qc, rejected)

    # ----------------------------------------------------------- provenance + write
    lock = json.loads((CFG / "ontologies.lock.json").read_text())
    adata.uns["provenance"] = {
        "schema_version": "1.0.0",
        "ontologies": {k: v["data_version"] for k, v in lock.items()},
        "config_sha256": {p.name: sha(p) for p in sorted(CFG.glob("*")) if p.is_file()},
        "packages": {p: version(p) for p in ["scanpy", "anndata", "harmonypy", "numpy", "scikit-learn"]},
        "integration": {"method": "harmonypy", "batch_key": "study_id", "n_hvg": N_HVG, "n_pcs": N_PCS, "seed": SEED},
        "excluded_samples": rejected,
        "genes_dropped_not_shared": dropped,
        "note": "layers['counts'] is primary. X_harmony is a derived analysis product, not repository truth.",
    }
    adata.write_h5ad(OUT / "harmonized.h5ad")
    (OUT / "data_qc.json").write_text(json.dumps(qc, indent=2, default=str) + "\n")
    print(f"[out] {OUT.relative_to(ROOT)}/harmonized.h5ad  ({adata.n_obs} cells x {adata.n_vars} genes)")


def plot(adata, emb, qc, rejected):
    study_colors = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a"}
    ct_colors = {"T cell": "#2a78d6", "B cell": "#eb6834", "monocyte": "#1baf7a", "natural killer cell": "#eda100"}
    ink, quiet, grid = "#0b0b0b", "#52514e", "#e7e6e2"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(2, 2, figsize=(11, 10.5), facecolor="#fcfcfb")
    titles = {"X_pca": "Before (PCA only)", "X_harmony": "After Harmony"}
    rng = np.random.default_rng(0)
    order = rng.permutation(adata.n_obs)                  # random draw order so no batch hides another
    for r, rep in enumerate(["X_pca", "X_harmony"]):
        xy = emb[rep][order]
        m = qc["integration"]["before (PCA)" if rep == "X_pca" else "after (Harmony)"]
        for c, (key, colors, name) in enumerate([("study_id", study_colors, "study"),
                                                 ("harmonized_label", ct_colors, "harmonized cell type")]):
            ax = axes[r, c]
            ax.set_facecolor("#fcfcfb")
            vals = adata.obs[key].to_numpy(dtype=str)[order]
            ax.scatter(xy[:, 0], xy[:, 1], s=4, c=[colors[v] for v in vals], linewidths=0, alpha=0.8,
                       rasterized=True)
            if key == "harmonized_label":                 # direct labels: identity is never colour-only
                for v in colors:
                    cx, cy = np.median(xy[vals == v], axis=0)
                    ax.text(cx, cy, v, ha="center", va="center", fontsize=9, color=ink,
                            bbox=dict(boxstyle="round,pad=0.25", fc="#fcfcfb", ec=grid, alpha=0.9))
            metric = (f"batch mixing (kNN entropy) {m['knn_batch_entropy (1=mixed)']:.2f}" if c == 0
                      else f"cell-type kNN purity {m['knn_celltype_purity (1=pure)']:.2f}")
            ax.set_title(f"{titles[rep]} · coloured by {name}\n{metric}", loc="left", fontsize=10.5, color=ink)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color(grid)
            handles = [plt.Line2D([], [], marker="o", ls="", color=col, markersize=7, label=lab)
                       for lab, col in colors.items()]
            ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper left", labelcolor=quiet,
                      bbox_to_anchor=(0, 0), ncol=len(handles), handletextpad=0.2, columnspacing=1.2)
    fig.suptitle("Harmony mixes the three studies while cell types stay separate",
                 x=0.02, ha="left", fontsize=14, fontweight="bold", color=ink)
    excl = f"; excluded pending curation: {', '.join(rejected)}" if rejected else ""
    fig.text(0.02, 0.944, f"Simulated PBMC data: {adata.obs.sample_id.nunique()} samples from 3 studies{excl}. "
             f"UMAP on {N_PCS} PCs from {N_HVG} batch-aware HVGs.", fontsize=9.5, color=quiet)
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=3)
    fig.savefig(OUT / "umap_before_after.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
