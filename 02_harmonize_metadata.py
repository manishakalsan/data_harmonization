"""Metadata harmonization: crosswalk -> normalize -> map to ontology -> validate.

Inputs : data/raw/study_*/samples.csv, counts.h5ad (obs only), config/*
Outputs: data/harmonized/
           samples.tsv           one row per sample, target schema + labels + valid flag
           cell_labels.tsv       author label -> CL term -> harmonized (rolled-up) term
           curation_log.tsv      every value decision with its status (provenance)
           curator_queue.tsv     what a human must resolve; never auto-guessed
           validation_report.json
"""
import json
import re
from pathlib import Path

import anndata as ad
import jsonschema
import pandas as pd
import yaml

from obo import Ontology

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config"
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "harmonized"

# field -> (which loaded ontology file to search, branch root the term must sit under)
FIELD_RULES = {
    "sex":       ("EFO", None),               # PATO terms ship inside EFO; schema enum restricts further
    "tissue":    ("UBERON", "UBERON:0001062"),  # anatomical entity
    "disease":   ("MONDO", "MONDO:0000001"),    # disease (PATO:0000461 'normal' handled below)
    "assay":     ("EFO", "EFO:0010183"),        # single cell library construction
    "cell_type": ("CL", "CL:0000000"),          # cell
}
NORMAL = "PATO:0000461"


def norm(text) -> str:
    return re.sub(r"\s+", " ", str(text).strip().lower())


def load_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text()) or {}


def parse_age(raw):
    """Return (min_years, max_years, status). Never invent a number."""
    s = norm(raw)
    if m := re.fullmatch(r"(\d+(?:\.\d+)?)\s*(?:y|yr|yrs|years?)?", s):
        v = float(m.group(1)); return v, v, "parsed_exact"
    if m := re.fullmatch(r"(\d+)\s*-\s*(\d+)\s*(?:y|years?)?", s):
        return float(m.group(1)), float(m.group(2)), "parsed_range"
    if m := re.fullmatch(r"p(\d+)y", s):              # ISO 8601 duration
        v = float(m.group(1)); return v, v, "parsed_iso8601"
    return None, None, "coarse_or_unparseable"         # e.g. "adult": keep raw, no number


class Mapper:
    def __init__(self, ontos: dict[str, Ontology], use_decisions: bool = True):
        self.ontos = ontos
        self.memory = load_yaml(CFG / "mapping_memory.yaml")
        self.overrides = load_yaml(CFG / "curator_decisions.yaml") if use_decisions else {}

    def check(self, field: str, tid: str) -> str | None:
        """Semantic validation of a chosen ID. Returns an error string or None."""
        onto_key, root = FIELD_RULES[field]
        if field == "disease" and tid == NORMAL:
            onto_key, root = "EFO", None
        term = self.ontos[onto_key].get(tid)
        if term is None:
            return f"{tid} not found in pinned {onto_key} ({self.ontos[onto_key].data_version})"
        if term.obsolete:
            repl = term.replaced_by or term.consider
            return f"{tid} is obsolete" + (f"; candidate replacement {repl}" if repl else "")
        if root and not self.ontos[onto_key].is_descendant_of(tid, root):
            return f"{tid} ({term.name}) is outside the allowed branch {root}"
        return None

    def map(self, field: str, raw) -> dict:
        key = norm(raw)
        onto_key, root = FIELD_RULES[field]
        tid, status = None, None
        if key in (self.overrides.get(field) or {}):
            tid, status = self.overrides[field][key], "curated_override"
        elif key in (self.memory.get(field) or {}):
            tid, status = self.memory[field][key], "curated_memory"
        else:
            hits = self.ontos[onto_key].exact_lookup(key, root)
            if len(hits) == 1:
                tid, status = hits[0], "auto_exact"
            elif len(hits) > 1:
                return dict(normalized=key, target_id=None, target_label="",
                            status="AMBIGUOUS", note=f"exact matches: {hits}")
            else:
                return dict(normalized=key, target_id=None, target_label="",
                            status="UNMAPPED", note="no memory entry, no exact label/synonym match")
        err = self.check(field, tid)
        onto = self.ontos["EFO" if tid == NORMAL else onto_key]
        label = onto.get(tid).name if onto.get(tid) else ""
        if err:
            return dict(normalized=key, target_id=None, target_label=label, status="INVALID_TERM", note=err)
        return dict(normalized=key, target_id=tid, target_label=label, status=status, note="")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--initial", action="store_true",
                    help="ignore curator_decisions.yaml; write only curator_queue_initial.tsv "
                         "(what a curator would have been sent before any decisions)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    lock = json.loads((CFG / "ontologies.lock.json").read_text())
    ontos = {k: Ontology(str(ROOT / "data" / "ontologies" / v["file"]), k) for k, v in lock.items()}
    for k, o in ontos.items():
        print(f"[onto] {k}: {len(o.terms):,} terms, {o.data_version}")

    schema = json.loads((CFG / "target_schema.json").read_text())
    validator = jsonschema.Draft202012Validator(schema)
    xw = load_yaml(CFG / "crosswalks.yaml")
    mapper = Mapper(ontos, use_decisions=not args.initial)

    samples, log, queue, cell_rows, report = [], [], [], [], {"samples": {}}
    for study, spec in xw["studies"].items():
        df = pd.read_csv(RAW / f"study_{study}" / spec["samples_file"], dtype=str)
        f = spec["fields"]
        missing = [c for c in f.values() if c not in df.columns]
        if missing:
            raise SystemExit(f"crosswalk for study {study} names columns not in the file: {missing}")

        for _, row in df.iterrows():
            rec = {"sample_id": row[f["sample_id"]], "study_id": study, "donor_id": row[f["donor_id"]],
                   **xw["target_constants"], "schema_version": schema["x-schema-version"]}
            labels = {}
            for field, target in [("sex", "sex_ontology_term_id"), ("tissue", "tissue_ontology_term_id"),
                                  ("disease", "disease_ontology_term_id"), ("assay", "assay_ontology_term_id")]:
                raw = row[f[field]]
                m = mapper.map(field, raw)
                rec[target] = m["target_id"]
                labels[field] = m["target_label"]
                log.append(dict(study=study, sample_id=rec["sample_id"], field=field, source_column=f[field],
                                source_value=raw, **m))
                if m["target_id"] is None:
                    queue.append(dict(study=study, sample_id=rec["sample_id"], field=field,
                                      source_value=raw, normalized=m["normalized"], status=m["status"], note=m["note"]))
            amin, amax, astatus = parse_age(row[f["age"]])
            rec.update(age_min_years=amin, age_max_years=amax, age_raw=str(row[f["age"]]))
            log.append(dict(study=study, sample_id=rec["sample_id"], field="age", source_column=f["age"],
                            source_value=row[f["age"]], normalized=norm(row[f["age"]]),
                            target_id=f"{amin}-{amax}", target_label="years", status=astatus, note=""))
            if astatus == "coarse_or_unparseable":
                queue.append(dict(study=study, sample_id=rec["sample_id"], field="age",
                                  source_value=row[f["age"]], normalized=norm(row[f["age"]]),
                                  status="COARSE", note="kept as raw text; map to a HsapDv stage, do not invent a number"))

            errors = [f"{'/'.join(map(str, e.path)) or '(root)'}: {e.message}"
                      for e in validator.iter_errors(rec)]
            if rec["age_min_years"] is not None and rec["age_max_years"] < rec["age_min_years"]:
                errors.append("age_max_years < age_min_years")
            report["samples"][rec["sample_id"]] = errors
            samples.append({**rec, "sex_label": labels["sex"], "tissue_label": labels["tissue"],
                            "disease_label": labels["disease"], "assay_label": labels["assay"],
                            "valid": not errors})

        # cell-type labels (from the count matrix obs, not the sample sheet)
        obs = ad.read_h5ad(RAW / f"study_{study}" / "counts.h5ad", backed="r").obs
        rollup = mapper.memory["cell_type_rollup"]
        for author, n in obs[spec["cells"]["label_col"]].value_counts().items():
            m = mapper.map("cell_type", author)
            harm = None
            if m["target_id"]:
                lineage = {m["target_id"]} | ontos["CL"].ancestors(m["target_id"])
                harm = next((r for r in rollup if r in lineage), None)
            cell_rows.append(dict(study=study, author_label=author, n_cells=int(n),
                                  cl_id=m["target_id"], cl_label=m["target_label"], evidence="author_label",
                                  harmonized_id=harm, harmonized_label=ontos["CL"].get(harm).name if harm else "",
                                  status=m["status"], note=m["note"]))
            if m["target_id"] is None:
                queue.append(dict(study=study, sample_id="*", field="cell_type", source_value=author,
                                  normalized=m["normalized"], status=m["status"], note=m["note"]))

    queue_cols = ["study", "sample_id", "field", "source_value", "normalized", "status", "note"]
    if args.initial:
        pd.DataFrame(queue, columns=queue_cols).to_csv(OUT / "curator_queue_initial.tsv", sep="\t", index=False)
        n_valid = sum(s["valid"] for s in samples)
        print(f"[meta] INITIAL pass (no curator decisions): {n_valid}/{len(samples)} samples valid, "
              f"{len(queue)} item(s) -> data/harmonized/curator_queue_initial.tsv")
        for q in queue:
            print(f"[meta]   queue: {q['sample_id']:>6} {q['field']:<8} {q['source_value']!r:<16} {q['status']}: {q['note']}")
        return
    pd.DataFrame(samples).to_csv(OUT / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(log).to_csv(OUT / "curation_log.tsv", sep="\t", index=False)
    pd.DataFrame(queue, columns=queue_cols).to_csv(OUT / "curator_queue.tsv", sep="\t", index=False)
    pd.DataFrame(cell_rows).to_csv(OUT / "cell_labels.tsv", sep="\t", index=False)

    statuses = pd.DataFrame(log)["status"].value_counts().to_dict()
    report["summary"] = {
        "schema_version": schema["x-schema-version"],
        "ontology_versions": {k: v.data_version for k, v in ontos.items()},
        "n_samples": len(samples), "n_valid": sum(s["valid"] for s in samples),
        "value_status_counts": statuses, "n_curator_queue": len(queue),
    }
    (OUT / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")

    s = report["summary"]
    print(f"[meta] samples valid: {s['n_valid']}/{s['n_samples']}")
    print(f"[meta] value decisions: {statuses}")
    print(f"[meta] curator queue: {len(queue)} item(s) -> data/harmonized/curator_queue.tsv")
    for sid, errs in report["samples"].items():
        for e in errs:
            print(f"[meta]   {sid}: {e}")


if __name__ == "__main__":
    main()
