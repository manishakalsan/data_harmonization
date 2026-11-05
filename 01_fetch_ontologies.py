"""Download pinned ontology releases and write a lockfile (data-version + SHA-256).

Run once; later runs reuse the cached files and fail loudly if a file's
checksum no longer matches the lockfile.
"""
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "ontologies"
LOCK = ROOT / "config" / "ontologies.lock.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    cfg = yaml.safe_load(open(ROOT / "config" / "ontologies.yaml"))["ontologies"]
    CACHE.mkdir(parents=True, exist_ok=True)
    lock = json.loads(LOCK.read_text()) if LOCK.exists() else {}

    for prefix, spec in cfg.items():
        dest = CACHE / Path(spec["url"]).name
        if not dest.exists():
            print(f"[fetch] {prefix}: {spec['url']}")
            urllib.request.urlretrieve(spec["url"], dest)
        digest = sha256(dest)

        with open(dest, encoding="utf-8") as fh:
            head = "".join(next(fh) for _ in range(5))
        if spec["expect_data_version"] not in head:
            sys.exit(f"[error] {prefix}: data-version does not contain {spec['expect_data_version']!r}")

        if prefix in lock and lock[prefix]["sha256"] != digest:
            sys.exit(f"[error] {prefix}: checksum changed since lockfile was written. "
                     "Delete the cached file or update the pin deliberately.")
        lock[prefix] = {"file": dest.name, "url": spec["url"],
                        "data_version": spec["expect_data_version"], "sha256": digest}
        print(f"[ok]    {prefix}: {spec['expect_data_version']}  sha256={digest[:12]}…")

    LOCK.write_text(json.dumps(lock, indent=2) + "\n")
    print(f"[lock]  {LOCK.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
