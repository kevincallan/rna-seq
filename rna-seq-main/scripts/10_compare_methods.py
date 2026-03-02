#!/usr/bin/env python3
"""
Step 10 -- Compare trimming methods.

Produces a consolidated comparison report across all enabled trimming
methods.  Covers:
  - Mapping rate comparison
  - Assigned reads in featureCounts
  - Correlation of normalised counts between methods
  - DEG overlap between methods (per contrast)
  - Genes with method-sensitive significance (padj flips)
"""

from __future__ import annotations

import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import (
    ensure_dirs,
    get_enabled_methods,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> List[Dict[str, str]]:
    """Read a TSV file into a list of dicts."""
    rows: List[Dict[str, str]] = []
    if not path.exists():
        return rows
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def load_de_genes(
    results_dir: Path,
    method: str,
    contrast: str,
    fdr: float,
) -> Set[str]:
    """Load set of significant gene IDs from a DESeq2 result."""
    # Try rscript-style path first
    de_path = results_dir / method / "deseq2" / contrast / "de_all.tsv"
    if not de_path.exists():
        # Try wrapper-style path
        de_path = results_dir / method / "deseq2" / "DESeq2.de_all.tsv"
    if not de_path.exists():
        return set()

    genes: Set[str] = set()
    with open(de_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        # Find padj column
        padj_idx = -1
        for i, h in enumerate(header):
            if "padj" in h.lower():
                padj_idx = i
                break
        if padj_idx < 0:
            return genes

        for line in fh:
            parts = line.strip().split("\t")
            try:
                padj = float(parts[padj_idx])
                if padj < fdr:
                    # Gene ID is first column (or row name)
                    genes.add(parts[0])
            except (ValueError, IndexError):
                continue
    return genes


def load_padj_dict(
    results_dir: Path,
    method: str,
    contrast: str,
) -> Dict[str, float]:
    """Load gene -> padj mapping from DE results."""
    de_path = results_dir / method / "deseq2" / contrast / "de_all.tsv"
    if not de_path.exists():
        de_path = results_dir / method / "deseq2" / "DESeq2.de_all.tsv"
    if not de_path.exists():
        return {}

    padj_map: Dict[str, float] = {}
    with open(de_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        padj_idx = -1
        for i, h in enumerate(header):
            if "padj" in h.lower():
                padj_idx = i
                break
        if padj_idx < 0:
            return padj_map

        for line in fh:
            parts = line.strip().split("\t")
            try:
                padj_map[parts[0]] = float(parts[padj_idx])
            except (ValueError, IndexError):
                continue
    return padj_map


def load_norm_counts(
    results_dir: Path,
    method: str,
    contrast: str,
) -> Dict[str, List[float]]:
    """Load normalised counts from DESeq2 output."""
    nc_path = results_dir / method / "deseq2" / contrast / "normalized_counts.tsv"
    if not nc_path.exists():
        return {}

    counts: Dict[str, List[float]] = {}
    with open(nc_path, encoding="utf-8") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            gene = parts[0]
            try:
                vals = [float(v) for v in parts[1:]]
                counts[gene] = vals
            except ValueError:
                continue
    return counts


def compute_correlation(
    counts_a: Dict[str, List[float]],
    counts_b: Dict[str, List[float]],
) -> Optional[float]:
    """Compute Pearson correlation of gene-level mean normalised counts."""
    common = sorted(set(counts_a.keys()) & set(counts_b.keys()))
    if len(common) < 10:
        return None

    import math

    means_a = [sum(counts_a[g]) / len(counts_a[g]) for g in common]
    means_b = [sum(counts_b[g]) / len(counts_b[g]) for g in common]

    n = len(common)
    mean_a = sum(means_a) / n
    mean_b = sum(means_b) / n

    cov = sum((a - mean_a) * (b - mean_b) for a, b in zip(means_a, means_b)) / n
    std_a = math.sqrt(sum((a - mean_a) ** 2 for a in means_a) / n)
    std_b = math.sqrt(sum((b - mean_b) ** 2 for b in means_b) / n)

    if std_a == 0 or std_b == 0:
        return None

    return cov / (std_a * std_b)


def summarize_featurecounts_options(
    fc_rows: List[Dict[str, str]],
) -> Dict[Tuple[str, str], Dict[str, float]]:
    """Aggregate assigned-read stats by (method, option_set)."""
    agg: Dict[Tuple[str, str], Dict[str, float]] = defaultdict(
        lambda: {"n": 0.0, "sum": 0.0, "min": float("inf"), "max": float("-inf")}
    )
    for row in fc_rows:
        method = row.get("method", "")
        opt = row.get("option_set", "")
        try:
            val = float(row.get("assigned_reads", "0"))
        except ValueError:
            continue
        key = (method, opt)
        a = agg[key]
        a["n"] += 1
        a["sum"] += val
        a["min"] = min(a["min"], val)
        a["max"] = max(a["max"], val)
    return agg


def count_degs_at_thresholds(
    de_all_path: Path,
    fdr_values: List[float],
    lfc_abs_threshold: float = 1.0,
) -> Dict[str, int]:
    """Count DEGs under multiple FDR / LFC thresholds from one de_all.tsv."""
    out: Dict[str, int] = {f"fdr_{f}": 0 for f in fdr_values}
    out["fdr_0.05_lfc1"] = 0

    if not de_all_path.exists():
        return out

    with open(de_all_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        padj_idx = -1
        lfc_idx = -1
        for i, h in enumerate(header):
            hl = h.lower()
            if "padj" in hl:
                padj_idx = i
            if "log2foldchange" in hl:
                lfc_idx = i

        if padj_idx < 0:
            return out

        for line in fh:
            parts = line.strip().split("\t")
            try:
                padj = float(parts[padj_idx])
            except (ValueError, IndexError):
                continue
            if padj != padj:  # NaN guard
                continue

            for f in fdr_values:
                if padj < f:
                    out[f"fdr_{f}"] += 1

            if lfc_idx >= 0:
                try:
                    lfc = float(parts[lfc_idx])
                    if padj < 0.05 and abs(lfc) >= lfc_abs_threshold:
                        out["fdr_0.05_lfc1"] += 1
                except (ValueError, IndexError):
                    continue

    return out


# ---------------------------------------------------------------------------
# Main comparison logic
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 10."""
    logger.info("=" * 60)
    logger.info("STEP 10: Compare trimming methods")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = methods_override or get_enabled_methods(cfg)
    comparisons = cfg.get("comparisons", [])
    fdr = cfg.get("deseq2", {}).get("fdr_threshold", 0.05)

    report_dir = results_dir / "reports"
    ensure_dirs(report_dir)

    comparison_lines: List[str] = []
    comparison_lines.append("# Trimming Method Comparison\n")

    # === 1. Mapping rate comparison =========================================
    comparison_lines.append("\n## Mapping Rates\n")
    mapping_data = read_tsv(results_dir / "mapping_summary.tsv")
    if mapping_data:
        comparison_lines.append("| Method | Sample | Uniquely Mapped % | Multi-Mapped % |")
        comparison_lines.append("|--------|--------|-------------------|----------------|")
        for row in mapping_data:
            comparison_lines.append(
                f"| {row.get('method', '')} | {row.get('sample', '')} "
                f"| {row.get('uniquely_mapped_pct', 'N/A')} "
                f"| {row.get('multi_mapped_pct', 'N/A')} |"
            )
    else:
        comparison_lines.append("*(Mapping summary not available)*\n")

    # === 2. Assigned reads comparison =======================================
    comparison_lines.append("\n## featureCounts Assigned Reads\n")
    fc_data = read_tsv(results_dir / "featurecounts_summary.tsv")
    if fc_data:
        comparison_lines.append("| Method | Option Set | Sample | Assigned Reads |")
        comparison_lines.append("|--------|-----------|--------|----------------|")
        for row in fc_data:
            comparison_lines.append(
                f"| {row.get('method', '')} | {row.get('option_set', '')} "
                f"| {row.get('sample', '')} | {row.get('assigned_reads', 'N/A')} |"
            )
    else:
        comparison_lines.append("*(featureCounts summary not available)*\n")

    # === 2b. featureCounts option-set impact ================================
    comparison_lines.append("\n## featureCounts Option-Set Impact\n")
    if fc_data:
        option_stats = summarize_featurecounts_options(fc_data)
        comparison_lines.append(
            "| Method | Option Set | N samples | Mean assigned reads | Min | Max |"
        )
        comparison_lines.append(
            "|--------|------------|-----------|---------------------|-----|-----|"
        )
        for (method, opt), stats in sorted(option_stats.items()):
            n = int(stats["n"])
            mean = stats["sum"] / stats["n"] if stats["n"] else 0.0
            min_v = int(stats["min"]) if n else 0
            max_v = int(stats["max"]) if n else 0
            comparison_lines.append(
                f"| {method} | {opt} | {n} | {mean:.1f} | {min_v} | {max_v} |"
            )

        # Delta versus default per method
        comparison_lines.append("\n### Delta vs default option set\n")
        comparison_lines.append(
            "| Method | Option Set | Mean assigned | Delta vs default |"
        )
        comparison_lines.append(
            "|--------|------------|---------------|------------------|"
        )
        by_method: Dict[str, Dict[str, float]] = defaultdict(dict)
        for (method, opt), stats in option_stats.items():
            mean = stats["sum"] / stats["n"] if stats["n"] else 0.0
            by_method[method][opt] = mean
        for method, opt_map in sorted(by_method.items()):
            base = opt_map.get("default")
            for opt, mean in sorted(opt_map.items()):
                if base is None:
                    delta = "N/A"
                else:
                    delta = f"{mean - base:+.1f}"
                comparison_lines.append(
                    f"| {method} | {opt} | {mean:.1f} | {delta} |"
                )
    else:
        comparison_lines.append("*(featureCounts option impact unavailable)*\n")

    # === 3. Normalised count correlation ====================================
    comparison_lines.append("\n## Normalised Count Correlation Between Methods\n")
    if len(methods) >= 2 and comparisons:
        contrast = comparisons[0]["name"]
        norm_counts_by_method: Dict[str, Dict[str, List[float]]] = {}
        for m in methods:
            nc = load_norm_counts(results_dir, m, contrast)
            if nc:
                norm_counts_by_method[m] = nc

        if len(norm_counts_by_method) >= 2:
            comparison_lines.append(f"*(Contrast: {contrast})*\n")
            comparison_lines.append("| Method A | Method B | Pearson r |")
            comparison_lines.append("|----------|----------|-----------|")
            mlist = sorted(norm_counts_by_method.keys())
            for i in range(len(mlist)):
                for j in range(i + 1, len(mlist)):
                    r = compute_correlation(
                        norm_counts_by_method[mlist[i]],
                        norm_counts_by_method[mlist[j]],
                    )
                    r_str = f"{r:.6f}" if r is not None else "N/A"
                    comparison_lines.append(
                        f"| {mlist[i]} | {mlist[j]} | {r_str} |"
                    )
        else:
            comparison_lines.append("*(Not enough normalised count data for correlation)*\n")
    else:
        comparison_lines.append("*(Need >= 2 methods and >= 1 contrast)*\n")

    # === 4. DEG overlap =====================================================
    comparison_lines.append("\n## DEG Overlap Between Methods\n")
    for contrast in comparisons:
        cname = contrast["name"]
        comparison_lines.append(f"\n### Contrast: {cname}\n")

        deg_sets: Dict[str, Set[str]] = {}
        for m in methods:
            degs = load_de_genes(results_dir, m, cname, fdr)
            deg_sets[m] = degs
            comparison_lines.append(f"- **{m}**: {len(degs)} DEGs (FDR < {fdr})")

        # Pairwise overlaps
        if len(methods) >= 2:
            comparison_lines.append("\n| Method A | Method B | A only | Shared | B only |")
            comparison_lines.append("|----------|----------|--------|--------|--------|")
            mlist = sorted(deg_sets.keys())
            for i in range(len(mlist)):
                for j in range(i + 1, len(mlist)):
                    a = deg_sets[mlist[i]]
                    b = deg_sets[mlist[j]]
                    shared = len(a & b)
                    a_only = len(a - b)
                    b_only = len(b - a)
                    comparison_lines.append(
                        f"| {mlist[i]} | {mlist[j]} | {a_only} | {shared} | {b_only} |"
                    )

    # === 5. Method-sensitive genes ==========================================
    comparison_lines.append("\n## Method-Sensitive Genes (padj flips)\n")
    comparison_lines.append(
        "Genes that are significant in one method but not another "
        f"(padj flip across {fdr}).\n"
    )

    for contrast in comparisons:
        cname = contrast["name"]
        comparison_lines.append(f"\n### {cname}\n")

        padj_maps: Dict[str, Dict[str, float]] = {}
        for m in methods:
            padj_maps[m] = load_padj_dict(results_dir, m, cname)

        if len(padj_maps) < 2:
            comparison_lines.append("*(Need >= 2 methods)*\n")
            continue

        # Find genes that flip significance between any pair
        all_genes = set()
        for pm in padj_maps.values():
            all_genes.update(pm.keys())

        flipped: List[Dict[str, str]] = []
        mlist = sorted(padj_maps.keys())
        for gene in sorted(all_genes):
            for i in range(len(mlist)):
                for j in range(i + 1, len(mlist)):
                    p_a = padj_maps[mlist[i]].get(gene)
                    p_b = padj_maps[mlist[j]].get(gene)
                    if p_a is not None and p_b is not None:
                        sig_a = p_a < fdr
                        sig_b = p_b < fdr
                        if sig_a != sig_b:
                            flipped.append({
                                "gene": gene,
                                "method_a": mlist[i],
                                "padj_a": f"{p_a:.2e}",
                                "method_b": mlist[j],
                                "padj_b": f"{p_b:.2e}",
                            })

        if flipped:
            comparison_lines.append(
                f"Found {len(flipped)} method-sensitive gene-pair instances:\n"
            )
            comparison_lines.append("| Gene | Method A | padj A | Method B | padj B |")
            comparison_lines.append("|------|----------|--------|----------|--------|")
            # Show first 50
            for entry in flipped[:50]:
                comparison_lines.append(
                    f"| {entry['gene']} | {entry['method_a']} | {entry['padj_a']} "
                    f"| {entry['method_b']} | {entry['padj_b']} |"
                )
            if len(flipped) > 50:
                comparison_lines.append(f"\n*(... and {len(flipped) - 50} more)*\n")
        else:
            comparison_lines.append("No method-sensitive genes detected.\n")

    # === 6. DE threshold sensitivity ========================================
    comparison_lines.append("\n## DE Threshold Sensitivity\n")
    comparison_lines.append(
        "Counts of significant genes under different thresholds "
        "(computed from existing `de_all.tsv`, no rerun required).\n"
    )
    fdr_values = [0.01, 0.05, 0.1]
    comparison_lines.append(
        "| Method | Contrast | FDR<0.01 | FDR<0.05 | FDR<0.1 | FDR<0.05 & |log2FC|>=1 |"
    )
    comparison_lines.append(
        "|--------|----------|----------|----------|---------|-----------------------|"
    )
    for method in methods:
        for contrast in comparisons:
            cname = contrast["name"]
            de_path = results_dir / method / "deseq2" / cname / "de_all.tsv"
            if not de_path.exists():
                de_path = results_dir / method / "deseq2" / "DESeq2.de_all.tsv"
            stats = count_degs_at_thresholds(de_path, fdr_values)
            comparison_lines.append(
                f"| {method} | {cname} | {stats['fdr_0.01']} | "
                f"{stats['fdr_0.05']} | {stats['fdr_0.1']} | "
                f"{stats['fdr_0.05_lfc1']} |"
            )

    # --- Write comparison report --------------------------------------------
    report_path = report_dir / "method_comparison.md"
    report_path.write_text("\n".join(comparison_lines), encoding="utf-8")
    logger.info("Method comparison report -> %s", report_path)

    logger.info("STEP 10 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 10: Compare methods")
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--methods", nargs="*", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg, methods_override=args.methods)
