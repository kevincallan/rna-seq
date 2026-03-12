#!/usr/bin/env python3
"""
Strandedness test for GSE104853: run steps 7–9 with -s 0, -s 1, -s 2 and compare.

Prerequisites:
  - results/gse104853_test exists (full pipeline already run with strandedness 0).
  - Run from pipeline directory: ./py scripts/run_strand_test.py

This script:
  1. Copies the main branch (cutadapt + STAR default) to results/gse104853_strand1
     and results/gse104853_strand2, updating paths in mapping_summary.
  2. Runs steps 7, 8, 9 with config_GSE104853_strand1.yaml and strand2.yaml.
  3. Collects metrics from baseline (gse104853_test) and from strand1, strand2.
  4. Prints a comparison table and recommends strandedness for config fix.

Main branch only: trim=cutadapt, mapper=star, mapper_option=default.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

# Pipeline root (parent of scripts/)
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
BASELINE_RUN_ID = "gse104853_test"
STRAND_RUN_IDS = ("gse104853_strand1", "gse104853_strand2")
STRAND_CONFIGS = (
    (STRAND_RUN_IDS[0], "config/config_GSE104853_strand1.yaml"),
    (STRAND_RUN_IDS[1], "config/config_GSE104853_strand2.yaml"),
)
MAIN_BRANCH = {"trim_method": "cutadapt", "mapper": "star", "mapper_option_set": "default"}
DE_CONTRAST_DIR = "cutadapt/deseq2/star/default/mirna125a_vs_ctrl"
GENES_OF_INTEREST = ("IL10RA", "CD163")


def results_dir(run_id: str) -> Path:
    return PIPELINE_ROOT / "results" / run_id


def copy_main_branch_from_baseline(target_run_id: str) -> None:
    """Create results/target_run_id with only main-branch data from baseline."""
    baseline = results_dir(BASELINE_RUN_ID)
    target = results_dir(target_run_id)
    if not baseline.exists():
        raise SystemExit(f"Baseline results not found: {baseline}")

    target.mkdir(parents=True, exist_ok=True)

    # samples.tsv
    samples_tsv = baseline / "samples.tsv"
    if samples_tsv.exists():
        (target / "samples.tsv").write_text(samples_tsv.read_text())

    # mapping_summary: only cutadapt + star + default, paths updated to target_run_id
    mapping_in = baseline / "mapping_summary.tsv"
    mapping_out = target / "mapping_summary.tsv"
    if mapping_in.exists():
        with open(mapping_in, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = [r for r in reader if (
                r.get("trim_method") == MAIN_BRANCH["trim_method"]
                and r.get("mapper") == MAIN_BRANCH["mapper"]
                and r.get("mapper_option_set") == MAIN_BRANCH["mapper_option_set"]
            )]
            fieldnames = reader.fieldnames
        if rows and fieldnames:
            for r in rows:
                for k in ("bam_path", "filtered_bam_path"):
                    if k in r and r[k]:
                        r[k] = r[k].replace(BASELINE_RUN_ID, target_run_id)
            with open(mapping_out, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)

    # Copy cutadapt/mapping/star/default (BAMs and indexes)
    src_map = baseline / "cutadapt" / "mapping" / "star" / "default"
    dst_map = target / "cutadapt" / "mapping" / "star" / "default"
    if src_map.exists():
        dst_map.mkdir(parents=True, exist_ok=True)
        for f in src_map.iterdir():
            if f.is_file():
                (dst_map / f.name).write_bytes(f.read_bytes())


def run_steps_7_8_9(config_path: str, run_id: str) -> bool:
    """Run pipeline steps 7, 8, 9 with given config and run_id. Returns True on success."""
    cfg = PIPELINE_ROOT / config_path
    if not cfg.exists():
        print(f"Config not found: {cfg}", file=sys.stderr)
        return False
    cmd = [
        str(PIPELINE_ROOT / "py"),
        "scripts/run_pipeline.py",
        "--config", str(cfg),
        "--dataset", "GSE104853",
        "--species", "human",
        "--run-id", run_id,
        "run", "--methods", "cutadapt", "--steps", "7", "8", "9",
    ]
    r = subprocess.run(cmd, cwd=str(PIPELINE_ROOT))
    return r.returncode == 0


def collect_metrics(run_id: str, strand_label: str) -> dict:
    """Collect assignment %, genes after filter, DEGs, IL10RA, CD163 for main branch."""
    out = {
        "setting": strand_label,
        "assigned_pct": "",
        "genes_after_filter": "",
        "DEGs_padj_lt_005": "",
        "IL10RA_log2FC_padj": "",
        "CD163_log2FC_padj": "",
        "comments": "",
    }
    root = results_dir(run_id)
    if not root.exists():
        out["comments"] = "missing"
        return out

    # featureCounts: % assigned for cutadapt/star/default, option_set default
    fc_summary = root / "featurecounts_summary.tsv"
    if fc_summary.exists():
        with open(fc_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = [r for r in reader if (
                r.get("trim_method") == "cutadapt"
                and r.get("mapper_option_set") == "default"
                and r.get("option_set") == "default"
            )]
        if rows:
            numeric_cols = [k for k in rows[0] if k not in (
                "trim_method", "mapper", "mapper_option_set", "option_set", "sample"
            )]
            total = 0
            assigned = 0
            for r in rows:
                for k in numeric_cols:
                    try:
                        v = int(r.get(k, 0))
                        total += v
                        if k == "Assigned":
                            assigned += v
                    except (ValueError, TypeError):
                        pass
            if total:
                out["assigned_pct"] = f"{100 * assigned / total:.1f}%"

    # filtering_summary: genes_out for cutadapt/star/default/default
    filt_summary = root / "filtering_summary.tsv"
    if filt_summary.exists():
        with open(filt_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                if (r.get("trim_method") == "cutadapt"
                        and r.get("mapper_option_set") == "default"
                        and r.get("option_set") == "default"):
                    out["genes_after_filter"] = r.get("genes_out", "")
                    break

    # DE: de_all.tsv in contrast dir
    de_all = root / DE_CONTRAST_DIR / "de_all.tsv"
    if de_all.exists():
        import pandas as pd
        df = pd.read_csv(de_all, sep="\t", index_col=0)
        if "padj" in df.columns:
            sig = df["padj"].notna() & (df["padj"] < 0.05)
            out["DEGs_padj_lt_005"] = str(int(sig.sum()))
        for gene in GENES_OF_INTEREST:
            if gene in df.index:
                row = df.loc[gene]
                lfc = row.get("log2FoldChange", "")
                padj = row.get("padj", "")
                if hasattr(lfc, "__round__"):
                    lfc = round(float(lfc), 3)
                if hasattr(padj, "__round__"):
                    padj = f"{float(padj):.2e}" if float(padj) < 0.01 else round(float(padj), 3)
                key = f"{gene}_log2FC_padj"
                if key in out:
                    out[key] = f"{lfc} / {padj}"
    return out


def main() -> None:
    print("Strandedness test for GSE104853 (main branch: cutadapt + STAR default)")
    print("=" * 70)

    # Prepare copies for strand1 and strand2
    for run_id in STRAND_RUN_IDS:
        d = results_dir(run_id)
        if d.exists():
            print(f"  {d} already exists; reusing (delete it to re-copy from baseline).")
        else:
            print(f"  Copying main branch from {BASELINE_RUN_ID} -> {run_id} ...")
            copy_main_branch_from_baseline(run_id)

    # Run steps 7, 8, 9 for strand1 and strand2
    for run_id, config_path in STRAND_CONFIGS:
        print(f"\n  Running steps 7,8,9 with {config_path} (run-id={run_id}) ...")
        ok = run_steps_7_8_9(config_path, run_id)
        if not ok:
            print(f"  WARNING: pipeline failed for {run_id}", file=sys.stderr)

    # Collect metrics (baseline = strand 0, then strand1, strand2)
    print("\nCollecting metrics ...")
    baseline_metrics = collect_metrics(BASELINE_RUN_ID, "strand 0 (baseline)")
    strand1_metrics = collect_metrics(STRAND_RUN_IDS[0], "strand 1")
    strand2_metrics = collect_metrics(STRAND_RUN_IDS[1], "strand 2")

    rows = [baseline_metrics, strand1_metrics, strand2_metrics]
    headers = ["setting", "assigned_pct", "genes_after_filter", "DEGs_padj_lt_005",
               "IL10RA_log2FC_padj", "CD163_log2FC_padj", "comments"]

    print("\n" + "=" * 70)
    print("Comparison table (main branch: cutadapt + STAR default)")
    print("=" * 70)
    col_widths = [max(len(str(r.get(h, ""))) for r in rows) + 2 for h in headers]
    col_widths = [max(c, len(h) + 2) for c, h in zip(col_widths, headers)]
    fmt = "  ".join(f"{{:<{w}}}" for w in col_widths)
    print(fmt.format(*headers))
    print("-" * 70)
    for r in rows:
        print(fmt.format(*[str(r.get(h, "")) for h in headers]))
    print("=" * 70)

    # Recommendation
    best = None
    best_degs = -1
    for m in rows:
        try:
            n = int(m["DEGs_padj_lt_005"])
            if n > best_degs:
                best_degs = n
                best = m["setting"]
        except (ValueError, TypeError):
            pass
    if best and best_degs >= 0:
        print(f"\nRecommendation: Use strandedness corresponding to '{best}' (DEGs = {best_degs}).")
        if "strand 1" in best:
            print("  -> Set featurecounts.strandedness: 1 in config_GSE104853.yaml")
        elif "strand 2" in best:
            print("  -> Set featurecounts.strandedness: 2 in config_GSE104853.yaml")
        else:
            print("  -> Keep featurecounts.strandedness: 0 or investigate further.")
    print("\nNext: Update config/config_GSE104853.yaml with chosen strandedness, then run")
    print("  ./py scripts/run_pipeline.py --config config/config_GSE104853.yaml \\")
    print("    --dataset GSE104853 --species human --run-id gse104853_fix_strand run \\")
    print("    --methods cutadapt --steps 7 8 9 10 11")
    print("to produce the fixed run and compare with baseline.")


if __name__ == "__main__":
    main()
