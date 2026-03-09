#!/usr/bin/env python3
"""
Metadata CSV Inspector -- exam-day helper.

Reads any SRA metadata CSV and prints:
  - Column names
  - Auto-detected organism, layout, assay type
  - Candidate condition columns (few unique values)
  - Ready-to-paste YAML snippet for column_mapping, subset_filters, comparisons

Usage:
    ./py scripts/inspect_metadata.py /data/GSEXXXXX/metadata.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Set

STANDARD_SRA_COLS: Set[str] = {
    "Run", "Assay Type", "AvgSpotLen", "Bases", "BioProject", "BioSample",
    "Bytes", "Center Name", "Consent", "DATASTORE filetype",
    "DATASTORE provider", "DATASTORE region", "Experiment",
    "GEO_Accession (exp)", "Instrument", "LibraryLayout", "LibrarySelection",
    "LibrarySource", "Organism", "Platform", "ReleaseDate", "create_date",
    "version", "Sample Name", "SRA Study",
}

SKIP_COLS: Set[str] = STANDARD_SRA_COLS | {
    "source_name", "cell_line", "tissue",
}

MAX_UNIQUE_FOR_CONDITION = 10
MIN_UNIQUE_FOR_CONDITION = 2


def safe_name(text: str) -> str:
    """Convert a long string to a short safe identifier."""
    text = text.strip()
    if len(text) <= 20 and re.match(r"^[A-Za-z0-9_\-]+$", text):
        return text
    words = re.findall(r"[A-Za-z0-9]+", text)
    short = "_".join(w.capitalize() for w in words[:4])
    if len(short) > 20:
        short = short[:20]
    return short or "condition"


def inspect(csv_path: str) -> None:
    path = Path(csv_path)
    if not path.exists():
        sys.exit(f"File not found: {csv_path}")

    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            sys.exit("Could not parse CSV headers")
        columns = list(reader.fieldnames)
        rows = list(reader)

    n_rows = len(rows)
    print(f"\n{'=' * 60}")
    print(f"Metadata Inspector: {path.name}")
    print(f"{'=' * 60}")
    print(f"Path:    {path}")
    print(f"Rows:    {n_rows}")
    print(f"Columns: {len(columns)}")

    # Auto-detect standard fields
    organisms = Counter(r.get("Organism", "") for r in rows)
    layouts = Counter(r.get("LibraryLayout", "") for r in rows)
    assays = Counter(r.get("Assay Type", "") for r in rows)

    print(f"\nOrganism:    {', '.join(f'{k} ({v})' for k, v in organisms.most_common())}")
    print(f"Layout:      {', '.join(f'{k} ({v})' for k, v in layouts.most_common())}")
    print(f"Assay Types: {', '.join(f'{k} ({v})' for k, v in assays.most_common())}")

    primary_organism = organisms.most_common(1)[0][0] if organisms else ""
    primary_layout = layouts.most_common(1)[0][0] if layouts else ""
    primary_assay = assays.most_common(1)[0][0] if assays else ""

    # Species suggestion
    org_lower = primary_organism.lower()
    if "mus musculus" in org_lower or "mouse" in org_lower:
        species_flag = "mouse"
    elif "homo sapiens" in org_lower or "human" in org_lower:
        species_flag = "human"
    else:
        species_flag = "UNKNOWN"
    print(f"Suggested --species: {species_flag}")

    # Find candidate condition columns
    print(f"\n--- Candidate condition columns ---")
    candidates: List[Dict[str, Any]] = []

    for col in columns:
        if col in SKIP_COLS:
            continue
        values = Counter(r.get(col, "").strip() for r in rows if r.get(col, "").strip())
        n_unique = len(values)
        if MIN_UNIQUE_FOR_CONDITION <= n_unique <= MAX_UNIQUE_FOR_CONDITION:
            candidates.append({"col": col, "values": values, "n_unique": n_unique})

    if not candidates:
        print("  (No obvious condition columns found -- check CSV manually)")
        print(f"\n  All columns: {', '.join(columns)}")
    else:
        for c in sorted(candidates, key=lambda x: x["n_unique"]):
            print(f"\n  {c['col']} ({c['n_unique']} unique values):")
            for val, count in c["values"].most_common():
                print(f"    - \"{val}\" ({count} samples)")

    # Pick best candidate (fewest unique values >= 2)
    best = None
    if candidates:
        best = min(candidates, key=lambda x: x["n_unique"])

    # Generate YAML snippet
    print(f"\n{'=' * 60}")
    print("Suggested YAML snippets")
    print(f"{'=' * 60}")

    if best:
        cond_col = best["col"]
        cond_values = sorted(best["values"].keys())
        print(f"\ncolumn_mapping:")
        print(f'  run_id_col: "Run"')
        print(f"  condition_cols:")
        print(f'    - "{cond_col}"')
        print(f"  condition_map:")
        for val in cond_values:
            short = safe_name(val)
            print(f'    "{val}": "{short}"')
        print(f'  replicate_strategy: "sorted_run_id"')

        # Suggest subset filter
        print(f"\nsubset_filters:")
        print(f"  default:")
        if primary_assay:
            print(f'    "Assay Type": "{primary_assay}"')

        # Suggest comparisons
        short_names = [safe_name(v) for v in cond_values]
        ref = short_names[0]
        print(f"\ndeseq2:")
        print(f'  reference_level: "{ref}"')
        print(f"\ncomparisons:")
        for sn in short_names[1:]:
            print(f'  - name: "{sn}_vs_{ref}"')
            print(f'    numerator: "{sn}"')
            print(f'    denominator: "{ref}"')

    # Run IDs and FASTQ discovery
    run_col = "Run"
    run_ids = [r.get(run_col, "").strip() for r in rows if r.get(run_col, "").strip()]
    print(f"\n--- Run IDs ({len(run_ids)}) ---")
    for rid in sorted(run_ids):
        print(f"  {rid}")

    dataset_dir = path.parent
    fastqs = sorted(dataset_dir.glob("*.fastq.gz"))
    fq_runs: Set[str] = set()
    for fq in fastqs:
        m = re.match(r"([SED]RR\d+)", fq.name)
        if m:
            fq_runs.add(m.group(1))
    matched = set(run_ids) & fq_runs
    print(f"\nFASTQ files in {dataset_dir}: {len(fastqs)} ({len(fq_runs)} unique runs)")
    print(f"Matched to metadata: {len(matched)}/{len(run_ids)}")

    # Suggested run command
    dataset_name = dataset_dir.name
    print(f"\n{'=' * 60}")
    print("Suggested pipeline command")
    print(f"{'=' * 60}")
    print(f"\n./py scripts/run_pipeline.py \\")
    print(f"  --config config/config_exam.yaml \\")
    print(f"  --dataset {dataset_name} \\")
    print(f"  --species {species_flag} \\")
    print(f"  run --methods none cutadapt")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect SRA metadata CSV for pipeline configuration",
    )
    parser.add_argument("csv_path", help="Path to metadata.csv")
    args = parser.parse_args()
    inspect(args.csv_path)


if __name__ == "__main__":
    main()
