#!/usr/bin/env bash
# Reproduce everything from scratch. Ontology downloads (~150 MB) are cached.
set -euo pipefail
cd "$(dirname "$0")/scripts"
python3 01_fetch_ontologies.py
python3 00_simulate_studies.py
python3 02_harmonize_metadata.py --initial
python3 02_harmonize_metadata.py
python3 03_harmonize_data.py
