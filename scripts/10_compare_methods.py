#!/usr/bin/env python3
"""
Step 10 -- Compare methods and parameter effects.

Produces a consolidated comparison report across all enabled trimming
methods and featureCounts option sets.  Covers:
  - Mapping rate comparison
  - Assigned reads in featureCounts
  - Correlation of normalised counts between analysis units
  - DEG overlap between analysis units (per contrast)
  - Genes with method-sensitive significance (padj flips)
  - DE threshold sensitivity
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
from src.analysis_unit import (
    AnalysisUnit,
    build_full_analysis_units,
    build_mapping_units,
    infer_analysis_units_from_de_summary,
    mapping_unit_label,
    resolve_de_dir,
    resolve_de_base_legacy,
    unit_label,
    write_selected_analysis,
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


def _resolve_de_all(results_dir: Path, au: AnalysisUnit, contrast: str) -> Optional[Path]:
    """Find de_all.tsv for an analysis unit + contrast.

    Checks the new count-option-aware path first, then legacy paths.
    """
    de_dir = resolve_de_dir(
        results_dir, au.method, au.mapper,
        au.mapper_option_set, au.count_option_set,
    )
    de_path = de_dir / contrast / "de_all.tsv"
    if de_path.exists():
        return de_path

    legacy_base = resolve_de_base_legacy(
        results_dir, au.method, au.mapper, au.mapper_option_set,
    )
    legacy_path = legacy_base / contrast / "de_all.tsv"
    if legacy_path.exists():
        return legacy_path
    wrapper_path = legacy_base / "DESeq2.de_all.tsv"
    if wrapper_path.exists():
        return wrapper_path

    return None


def load_de_genes(
    results_dir: Path,
    au: AnalysisUnit,
    contrast: str,
    fdr: float,
) -> Set[str]:
    """Load set of significant gene IDs from a DESeq2 result."""
    de_path = _resolve_de_all(results_dir, au, contrast)
    if de_path is None:
        return set()

    genes: Set[str] = set()
    with open(de_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
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
                    genes.add(parts[0])
            except (ValueError, IndexError):
                continue
    return genes


def load_padj_dict(
    results_dir: Path,
    au: AnalysisUnit,
    contrast: str,
) -> Dict[str, float]:
    """Load gene -> padj mapping from DE results."""
    de_path = _resolve_de_all(results_dir, au, contrast)
    if de_path is None:
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
    au: AnalysisUnit,
    contrast: str,
) -> Dict[str, List[float]]:
    """Load normalised counts from DESeq2 output."""
    de_dir = resolve_de_dir(
        results_dir, au.method, au.mapper,
        au.mapper_option_set, au.count_option_set,
    )
    nc_path = de_dir / contrast / "normalized_counts.tsv"
    if not nc_path.exists():
        legacy_base = resolve_de_base_legacy(
            results_dir, au.method, au.mapper, au.mapper_option_set,
        )
        nc_path = legacy_base / contrast / "normalized_counts.tsv"
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
) -> Dict[Tuple[str, str, str, str], Dict[str, float]]:
    """Aggregate assigned-read stats by (method, mapper, mapper_option_set, option_set)."""
    agg: Dict[Tuple[str, str, str, str], Dict[str, float]] = defaultdict(
        lambda: {"n": 0.0, "sum": 0.0, "min": float("inf"), "max": float("-inf"),
                 "total_sum": 0.0}
    )
    for row in fc_rows:
        method = row.get("trim_method", row.get("method", ""))
        mapper = row.get("mapper", "star")
        mapper_opt = row.get("mapper_option_set", "default")
        opt = row.get("option_set", "")
        try:
            val = float(row.get("Assigned", row.get("assigned_reads", "0")))
        except ValueError:
            continue
        # Compute total reads for assigned % (sum all numeric status cols)
        total = 0.0
        for k, v in row.items():
            if k in ("trim_method", "method", "mapper", "mapper_option_set",
                      "option_set", "sample", "Assigned_pct"):
                continue
            try:
                total += float(v)
            except (ValueError, TypeError):
                continue

        key = (method, mapper, mapper_opt, opt)
        a = agg[key]
        a["n"] += 1
        a["sum"] += val
        a["min"] = min(a["min"], val)
        a["max"] = max(a["max"], val)
        a["total_sum"] += total
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
# Selected-branch logic
# ---------------------------------------------------------------------------

def _select_primary_branch(
    cfg: Dict[str, Any],
    results_dir: Path,
    units: List[AnalysisUnit],
    fc_data: List[Dict[str, str]],
    comparisons: List[Dict[str, str]],
    fdr: float,
) -> Optional[AnalysisUnit]:
    """Select the primary analysis branch and write selected_analysis.tsv.

    Selection rules:
    1. If config has ``selected_branch``, use that.
    2. Otherwise, among units with valid DE outputs and non-multimapper
       count options, pick the one with the highest mean assigned-read %.
    3. Tie-break: prefer 'default' mapper_option_set, then 'strict'.
    """
    sb = cfg.get("selected_branch")
    if sb:
        manual = AnalysisUnit(
            method=sb.get("trim_method", ""),
            mapper=sb.get("mapper", "star"),
            mapper_option_set=sb.get("mapper_option_set", "default"),
            count_option_set=sb.get("count_option_set", "default"),
        )
        write_selected_analysis(results_dir, manual, reason="manual_config")
        return manual

    if not units:
        return None

    option_sets = cfg.get("featurecounts", {}).get("option_sets", {})
    multimapper_opts = {
        name for name, params in option_sets.items()
        if isinstance(params, dict) and params.get("M", False)
    }

    # Build assigned-% aggregates from fc_data
    assigned_pct: Dict[Tuple[str, str, str, str], float] = {}
    if fc_data:
        option_stats = summarize_featurecounts_options(fc_data)
        for (method, mapper, mapper_opt, opt), stats in option_stats.items():
            pct = (100.0 * stats["sum"] / stats["total_sum"]
                   if stats["total_sum"] > 0 else 0.0)
            assigned_pct[(method, mapper, mapper_opt, opt)] = pct

    # Filter to units with valid DE output and non-multimapper
    candidates: List[Tuple[float, int, AnalysisUnit]] = []
    for au in units:
        if au.count_option_set in multimapper_opts:
            continue

        has_de = False
        for contrast in comparisons:
            de_path = _resolve_de_all(results_dir, au, contrast["name"])
            if de_path is not None and de_path.exists():
                has_de = True
                break
        if not has_de:
            continue

        pct = assigned_pct.get(
            (au.method, au.mapper, au.mapper_option_set, au.count_option_set), 0.0
        )
        # Tie-break ordering: prefer "default" mapper_opt then "strict" count_opt
        tiebreak = 0
        if au.mapper_option_set == "default":
            tiebreak += 2
        if au.count_option_set == "strict":
            tiebreak += 1

        candidates.append((pct, tiebreak, au))

    if not candidates:
        # Fall back to first unit
        write_selected_analysis(results_dir, units[0], reason="fallback_first_unit")
        return units[0]

    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best = candidates[0][2]
    write_selected_analysis(
        results_dir, best,
        reason="highest_assigned_pct_with_valid_de",
    )
    return best


# ---------------------------------------------------------------------------
# Main comparison logic
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 10."""
    logger.info("=" * 60)
    logger.info("STEP 10: Compare methods and parameter effects")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = methods_override or get_enabled_methods(cfg)
    comparisons = cfg.get("comparisons", [])
    fdr = cfg.get("deseq2", {}).get("fdr_threshold", 0.05)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})

    units = infer_analysis_units_from_de_summary(results_dir, methods)
    if not units:
        units = build_full_analysis_units(results_dir, methods, option_sets)

    report_dir = results_dir / "reports"
    ensure_dirs(report_dir)

    comparison_lines: List[str] = []
    comparison_lines.append("# Method and Parameter Comparison\n")
    comparison_lines.append(
        f"Analysis units: {', '.join(au.label for au in units)}\n"
    )

    # === 1. Mapping rate comparison =========================================
    comparison_lines.append("\n## Mapping Rates\n")
    mapping_data = read_tsv(results_dir / "mapping_summary.tsv")
    if mapping_data:
        comparison_lines.append(
            "| Trim Method | Mapper | Mapper Option | Sample | Uniquely Mapped % | Multi-Mapped % |"
        )
        comparison_lines.append(
            "|-------------|--------|---------------|--------|-------------------|----------------|"
        )
        for row in mapping_data:
            method = row.get("trim_method", row.get("method", ""))
            comparison_lines.append(
                f"| {method} | {row.get('mapper', 'star')} | {row.get('mapper_option_set', 'default')} "
                f"| {row.get('sample', '')} "
                f"| {row.get('uniquely_mapped_pct', 'N/A')} "
                f"| {row.get('multi_mapped_pct', 'N/A')} |"
            )
    else:
        comparison_lines.append("*(Mapping summary not available)*\n")

    # === 1b. Mapper option impact ===========================================
    comparison_lines.append("\n## Mapper Option Impact\n")
    if mapping_data:
        agg: Dict[Tuple[str, str, str], Dict[str, float]] = defaultdict(
            lambda: {"n": 0.0, "uniq_sum": 0.0}
        )
        for row in mapping_data:
            method = row.get("trim_method", row.get("method", ""))
            mapper = row.get("mapper", "star")
            mapper_opt = row.get("mapper_option_set", "default")
            uniq_pct = row.get("uniquely_mapped_pct", "").replace("%", "")
            try:
                uniq_val = float(uniq_pct)
            except ValueError:
                continue
            key = (method, mapper, mapper_opt)
            agg[key]["n"] += 1
            agg[key]["uniq_sum"] += uniq_val
        comparison_lines.append(
            "| Trim Method | Mapper | Mapper Option | N samples | Mean uniquely mapped % |"
        )
        comparison_lines.append(
            "|-------------|--------|---------------|-----------|------------------------|"
        )
        for (method, mapper, mapper_opt), stats in sorted(agg.items()):
            mean = stats["uniq_sum"] / stats["n"] if stats["n"] else 0.0
            comparison_lines.append(
                f"| {method} | {mapper} | {mapper_opt} | {int(stats['n'])} | {mean:.2f}% |"
            )
    else:
        comparison_lines.append("*(Mapper option impact unavailable)*\n")

    # === 2. Assigned reads comparison =======================================
    comparison_lines.append("\n## featureCounts Assigned Reads\n")
    fc_data = read_tsv(results_dir / "featurecounts_summary.tsv")
    if fc_data:
        comparison_lines.append(
            "| Trim Method | Mapper | Mapper Option | Count Option Set | Sample | Assigned Reads |"
        )
        comparison_lines.append(
            "|-------------|--------|---------------|------------------|--------|----------------|"
        )
        for row in fc_data:
            method = row.get("trim_method", row.get("method", ""))
            comparison_lines.append(
                f"| {method} | {row.get('mapper', 'star')} | {row.get('mapper_option_set', 'default')} "
                f"| {row.get('option_set', '')} "
                f"| {row.get('sample', '')} | {row.get('Assigned', row.get('assigned_reads', 'N/A'))} |"
            )
    else:
        comparison_lines.append("*(featureCounts summary not available)*\n")

    # === 2b. featureCounts option-set impact ================================
    comparison_lines.append("\n## featureCounts Option-Set Impact\n")
    if fc_data:
        option_stats = summarize_featurecounts_options(fc_data)
        comparison_lines.append(
            "| Trim Method | Mapper | Mapper Option | Count Option Set | N samples "
            "| Mean assigned reads | Mean assigned % | Min | Max |"
        )
        comparison_lines.append(
            "|-------------|--------|---------------|------------------|----------"
            "|---------------------|-----------------|-----|-----|"
        )
        for (method, mapper, mapper_opt, opt), stats in sorted(option_stats.items()):
            n = int(stats["n"])
            mean = stats["sum"] / stats["n"] if stats["n"] else 0.0
            mean_pct = (100.0 * stats["sum"] / stats["total_sum"]
                        if stats["total_sum"] > 0 else 0.0)
            min_v = int(stats["min"]) if n else 0
            max_v = int(stats["max"]) if n else 0
            comparison_lines.append(
                f"| {method} | {mapper} | {mapper_opt} | {opt} | {n} "
                f"| {mean:.1f} | {mean_pct:.1f}% | {min_v} | {max_v} |"
            )

        # Delta versus default per trim+mapper+mapper_option
        comparison_lines.append("\n### Delta vs default option set\n")
        comparison_lines.append(
            "| Trim Method | Mapper | Mapper Option | Count Option Set | Mean assigned | Delta vs default |"
        )
        comparison_lines.append(
            "|-------------|--------|---------------|------------------|---------------|------------------|"
        )
        by_unit: Dict[Tuple[str, str, str], Dict[str, float]] = defaultdict(dict)
        for (method, mapper, mapper_opt, opt), stats in option_stats.items():
            mean = stats["sum"] / stats["n"] if stats["n"] else 0.0
            by_unit[(method, mapper, mapper_opt)][opt] = mean
        for (method, mapper, mapper_opt), opt_map in sorted(by_unit.items()):
            base = opt_map.get("default")
            for opt, mean in sorted(opt_map.items()):
                if base is None:
                    delta = "N/A"
                else:
                    delta = f"{mean - base:+.1f}"
                comparison_lines.append(
                    f"| {method} | {mapper} | {mapper_opt} | {opt} | {mean:.1f} | {delta} |"
                )
    else:
        comparison_lines.append("*(featureCounts option impact unavailable)*\n")

    # === 3. Normalised count correlation ====================================
    comparison_lines.append("\n## Normalised Count Correlation Between Analysis Units\n")
    if len(units) >= 2 and comparisons:
        contrast = comparisons[0]["name"]
        norm_counts_by_unit: Dict[str, Dict[str, List[float]]] = {}
        for au in units:
            nc = load_norm_counts(results_dir, au, contrast)
            if nc:
                norm_counts_by_unit[au.label] = nc

        if len(norm_counts_by_unit) >= 2:
            comparison_lines.append(f"*(Contrast: {contrast})*\n")
            comparison_lines.append("| Unit A | Unit B | Pearson r |")
            comparison_lines.append("|----------|----------|-----------|")
            mlist = sorted(norm_counts_by_unit.keys())
            for i in range(len(mlist)):
                for j in range(i + 1, len(mlist)):
                    r = compute_correlation(
                        norm_counts_by_unit[mlist[i]],
                        norm_counts_by_unit[mlist[j]],
                    )
                    r_str = f"{r:.6f}" if r is not None else "N/A"
                    comparison_lines.append(
                        f"| {mlist[i]} | {mlist[j]} | {r_str} |"
                    )
        else:
            comparison_lines.append("*(Not enough normalised count data for correlation)*\n")
    else:
        comparison_lines.append("*(Need >= 2 analysis units and >= 1 contrast)*\n")

    # === 4. DEG overlap =====================================================
    comparison_lines.append("\n## DEG Overlap Between Analysis Units\n")
    for contrast in comparisons:
        cname = contrast["name"]
        comparison_lines.append(f"\n### Contrast: {cname}\n")

        deg_sets: Dict[str, Set[str]] = {}
        for au in units:
            degs = load_de_genes(results_dir, au, cname, fdr)
            deg_sets[au.label] = degs
            comparison_lines.append(f"- **{au.label}**: {len(degs)} DEGs (FDR < {fdr})")

        if len(units) >= 2:
            comparison_lines.append("\n| Unit A | Unit B | A only | Shared | B only |")
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
        "Genes that are significant in one analysis unit but not another "
        f"(padj flip across {fdr}).\n"
    )

    for contrast in comparisons:
        cname = contrast["name"]
        comparison_lines.append(f"\n### {cname}\n")

        padj_maps: Dict[str, Dict[str, float]] = {}
        for au in units:
            padj_maps[au.label] = load_padj_dict(results_dir, au, cname)

        if len(units) < 2:
            comparison_lines.append("*(Need >= 2 analysis units)*\n")
            continue

        all_genes: Set[str] = set()
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
        "| Trim Method | Mapper | Mapper Opt | Count Opt | Contrast "
        "| FDR<0.01 | FDR<0.05 | FDR<0.1 | FDR<0.05 & |log2FC|>=1 |"
    )
    comparison_lines.append(
        "|-------------|--------|-----------|-----------|----------"
        "|----------|----------|---------|-----------------------|"
    )
    for au in units:
        for contrast in comparisons:
            cname = contrast["name"]
            de_path = _resolve_de_all(results_dir, au, cname)
            if de_path is None:
                de_path = Path("/dev/null")  # count_degs_at_thresholds handles missing
            stats = count_degs_at_thresholds(de_path, fdr_values)
            comparison_lines.append(
                f"| {au.method} | {au.mapper} | {au.mapper_option_set} "
                f"| {au.count_option_set} | {cname} | {stats['fdr_0.01']} | "
                f"{stats['fdr_0.05']} | {stats['fdr_0.1']} | "
                f"{stats['fdr_0.05_lfc1']} |"
            )

    # --- Write comparison report --------------------------------------------
    report_path = report_dir / "method_comparison.md"
    report_path.write_text("\n".join(comparison_lines), encoding="utf-8")
    logger.info("Method comparison report -> %s", report_path)

    # === 7. Select primary analysis branch ==================================
    selected = _select_primary_branch(cfg, results_dir, units, fc_data, comparisons, fdr)
    if selected:
        logger.info("Selected primary branch: %s", selected.label)
    else:
        logger.warning("Could not determine a selected primary branch.")

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
