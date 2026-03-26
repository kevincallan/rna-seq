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
    get_effective_trim_methods,
    get_trim_config_summary,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)
from src.analysis_unit import (
    AnalysisUnit,
    build_full_analysis_units,
    infer_analysis_units_from_de_summary,
    resolve_de_dir,
    resolve_de_base_legacy,
    write_selected_analysis,
    write_selected_count_comparison,
    write_selected_visualisation,
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

def _to_float(value: str) -> Optional[float]:
    try:
        return float(str(value).replace("%", "").strip())
    except (TypeError, ValueError):
        return None


def _build_unit_metrics(
    units: List[AnalysisUnit],
    mapping_data: List[Dict[str, str]],
    fc_data: List[Dict[str, str]],
    filtering_data: List[Dict[str, str]],
    de_summary_data: List[Dict[str, str]],
    redundancy_data: List[Dict[str, str]],
) -> Dict[AnalysisUnit, Dict[str, Any]]:
    metrics: Dict[AnalysisUnit, Dict[str, Any]] = {
        au: {
            "uniq_vals": [],
            "assigned_vals": [],
            "unassigned_no_features_vals": [],
            "unassigned_multimapping_vals": [],
            "genes_in_vals": [],
            "genes_out_vals": [],
            "de_tested": 0,
            "de_sig": 0,
            "de_completed_rows": 0,
            "de_any_valid": False,
            "is_redundant": False,
            "canonical_label": "",
        }
        for au in units
    }

    map_by_key = defaultdict(list)
    for row in mapping_data:
        key = (
            row.get("trim_method", row.get("method", "")),
            row.get("mapper", "star"),
            row.get("mapper_option_set", "default"),
        )
        map_by_key[key].append(row)
    for au in units:
        for row in map_by_key[(au.method, au.mapper, au.mapper_option_set)]:
            uniq = _to_float(row.get("uniquely_mapped_pct", ""))
            if uniq is not None:
                metrics[au]["uniq_vals"].append(uniq)

    for row in fc_data:
        key = (
            row.get("trim_method", row.get("method", "")),
            row.get("mapper", "star"),
            row.get("mapper_option_set", "default"),
            row.get("option_set", row.get("count_option_set", "default")),
        )
        for au in units:
            if (au.method, au.mapper, au.mapper_option_set, au.count_option_set) != key:
                continue
            assigned = _to_float(row.get("Assigned_pct", ""))
            if assigned is not None:
                metrics[au]["assigned_vals"].append(assigned)
            no_feat = _to_float(
                row.get("Unassigned_NoFeatures", row.get("Unassigned_No_Features", ""))
            )
            if no_feat is not None:
                metrics[au]["unassigned_no_features_vals"].append(no_feat)
            multi = _to_float(row.get("Unassigned_MultiMapping", ""))
            if multi is not None:
                metrics[au]["unassigned_multimapping_vals"].append(multi)

    for row in filtering_data:
        key = (
            row.get("trim_method", row.get("method", "")),
            row.get("mapper", "star"),
            row.get("mapper_option_set", "default"),
            row.get("option_set", row.get("count_option_set", "default")),
        )
        for au in units:
            if (au.method, au.mapper, au.mapper_option_set, au.count_option_set) != key:
                continue
            genes_in = _to_float(row.get("genes_in", ""))
            genes_out = _to_float(row.get("genes_out", ""))
            if genes_in is not None:
                metrics[au]["genes_in_vals"].append(genes_in)
            if genes_out is not None:
                metrics[au]["genes_out_vals"].append(genes_out)

    for row in de_summary_data:
        key = (
            row.get("trim_method", row.get("method", "")),
            row.get("mapper", "star"),
            row.get("mapper_option_set", "default"),
            row.get("count_option_set", "default"),
        )
        status = str(row.get("status", "completed")).strip().lower()
        for au in units:
            if (au.method, au.mapper, au.mapper_option_set, au.count_option_set) != key:
                continue
            if status == "completed":
                tested = int(float(row.get("total_tested", "0") or 0))
                sig = int(float(row.get("sig_total", "0") or 0))
                metrics[au]["de_tested"] += tested
                metrics[au]["de_sig"] += sig
                metrics[au]["de_completed_rows"] += 1
                metrics[au]["de_any_valid"] = True

    for row in redundancy_data:
        key = (
            row.get("trim_method", ""),
            row.get("mapper", "star"),
            row.get("mapper_option_set", "default"),
            row.get("count_option_set", "default"),
        )
        is_redundant = str(row.get("is_redundant", "")).strip().lower() in {"1", "true", "yes"}
        if not is_redundant:
            continue
        canonical = "/".join([
            row.get("canonical_trim_method", ""),
            row.get("canonical_mapper", ""),
            row.get("canonical_mapper_option_set", ""),
            row.get("canonical_count_option_set", ""),
        ])
        for au in units:
            if (au.method, au.mapper, au.mapper_option_set, au.count_option_set) == key:
                metrics[au]["is_redundant"] = True
                metrics[au]["canonical_label"] = canonical

    return metrics


def _mean(values: List[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _apply_selection_policy(
    cfg: Dict[str, Any],
    results_dir: Path,
    units: List[AnalysisUnit],
    metrics: Dict[AnalysisUnit, Dict[str, Any]],
) -> Tuple[Optional[AnalysisUnit], Optional[AnalysisUnit], str, str]:
    """Apply the explicit selection policy from config.

    Overrides all older manual selectors (selected_count_comparison,
    selected_branch, selected_visualisation) when selection: is present.

    Returns (primary_au, comparison_au, primary_reason, comparison_reason).
    Raises RuntimeError if the primary branch cannot be found or is redundant.
    """
    sel = cfg["selection"]
    trim = sel["preferred_trim_method"]
    mapper = sel["preferred_mapper"]
    mapper_opt = sel["preferred_mapper_option_set"]
    primary_cos = sel["primary_count_option_set"]
    comparison_cos_list = sel.get("comparison_count_option_sets", [])

    candidates = [
        u for u in units
        if u.method == trim
        and u.mapper == mapper
        and u.mapper_option_set == mapper_opt
    ]
    logger.info(
        "Selection policy: %d candidates matching %s/%s/%s (from %d total units)",
        len(candidates), trim, mapper, mapper_opt, len(units),
    )

    def _is_selectable(au: AnalysisUnit) -> bool:
        m = metrics.get(au)
        if m is None:
            return False
        return not m["is_redundant"]

    def _warn_if_no_de(au: AnalysisUnit, role: str) -> None:
        m = metrics.get(au)
        if m is not None and not m["de_any_valid"]:
            logger.warning(
                "Selection policy: %s branch %s has no valid DE results. "
                "Proceeding anyway -- check upstream steps.",
                role, au.label,
            )

    # --- Primary ---
    primary = next(
        (u for u in candidates if u.count_option_set == primary_cos), None
    )
    if primary is None or not _is_selectable(primary):
        raise RuntimeError(
            f"Selection policy: primary branch {trim}/{mapper}/{mapper_opt}/"
            f"{primary_cos} not found in inferred units or is redundant. "
            f"Candidates: {[u.label for u in candidates]}"
        )
    _warn_if_no_de(primary, "primary")

    primary_reason = f"selection_policy:primary={primary.label}"
    write_selected_analysis(results_dir, primary, reason=primary_reason)
    write_selected_visualisation(results_dir, primary, reason=primary_reason)
    logger.info(
        "Selection policy: primary/analysis branch = %s (%s)",
        primary.label, primary_reason,
    )
    logger.info(
        "Selection policy: visualisation branch    = %s (%s)",
        primary.label, primary_reason,
    )

    # --- Comparison ---
    comparison = next(
        (u for u in candidates
         if u.count_option_set != primary_cos
         and u.count_option_set in comparison_cos_list
         and _is_selectable(u)),
        None,
    )

    if comparison is not None:
        _warn_if_no_de(comparison, "comparison")
        comp_reason = f"selection_policy:comparison={comparison.label}"
        write_selected_count_comparison(results_dir, comparison, reason=comp_reason)
        logger.info(
            "Selection policy: comparison branch       = %s (%s)",
            comparison.label, comp_reason,
        )
    else:
        comp_reason = "selection_policy:comparison_not_available"
        logger.warning(
            "Selection policy: no valid comparison branch found among %s "
            "(candidates: %s). selected_count_comparison.tsv NOT written.",
            comparison_cos_list, [u.label for u in candidates],
        )

    return primary, comparison, primary_reason, comp_reason


def _select_count_comparison_branch(
    cfg: Dict[str, Any],
    results_dir: Path,
    units: List[AnalysisUnit],
    metrics: Dict[AnalysisUnit, Dict[str, Any]],
) -> Tuple[Optional[AnalysisUnit], str]:
    manual = cfg.get("selected_count_comparison") or cfg.get("selected_branch")
    if manual:
        au = AnalysisUnit(
            method=manual.get("trim_method", ""),
            mapper=manual.get("mapper", "star"),
            mapper_option_set=manual.get("mapper_option_set", "default"),
            count_option_set=manual.get("count_option_set", "default"),
        )
        write_selected_count_comparison(results_dir, au, reason="manual_config")
        return au, "manual_config"

    ranked: List[Tuple[Tuple[float, float, float, int, int], AnalysisUnit]] = []
    for au in units:
        m = metrics[au]
        if m["is_redundant"]:
            continue
        valid = 1 if m["de_any_valid"] else 0
        prefer = 0
        if au.method == "none" and au.mapper == "star" and au.mapper_option_set == "default":
            prefer = 4
        assigned = _mean(m["assigned_vals"])
        retained = _mean(m["genes_out_vals"])
        non_redundant = 1
        ranked.append(((valid, prefer, assigned, int(retained), non_redundant), au))

    if not ranked:
        return None, "no_candidates"

    ranked.sort(key=lambda x: x[0], reverse=True)
    chosen = ranked[0][1]
    write_selected_count_comparison(
        results_dir, chosen, reason="auto_prefer_none_star_default_with_valid_de"
    )
    return chosen, "auto_prefer_none_star_default_with_valid_de"


def _select_visualisation_branch(
    cfg: Dict[str, Any],
    results_dir: Path,
    units: List[AnalysisUnit],
) -> Tuple[Optional[AnalysisUnit], str]:
    manual = cfg.get("selected_visualisation")
    if manual:
        au = AnalysisUnit(
            method=manual.get("trim_method", ""),
            mapper=manual.get("mapper", "star"),
            mapper_option_set=manual.get("mapper_option_set", "strict_unique"),
            count_option_set=manual.get("count_option_set", "default"),
        )
        write_selected_visualisation(results_dir, au, reason="manual_config")
        return au, "manual_config"

    preferred = AnalysisUnit("none", "star", "strict_unique", "default")
    if preferred in units:
        write_selected_visualisation(results_dir, preferred, reason="auto_prefer_none_star_strict_unique_default")
        return preferred, "auto_prefer_none_star_strict_unique_default"

    fallback = next((u for u in units if u.method == "none" and u.mapper == "star"), None)
    if fallback is None and units:
        fallback = units[0]
    if fallback:
        write_selected_visualisation(results_dir, fallback, reason="fallback_available_star_branch")
    return fallback, "fallback_available_star_branch"


# ---------------------------------------------------------------------------
# Main comparison logic
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 10."""
    logger.info("=" * 60)
    logger.info("STEP 10: Compare methods and parameter effects")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = get_effective_trim_methods(cfg, methods_override)
    trim_cfg = get_trim_config_summary(cfg)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})

    units = infer_analysis_units_from_de_summary(results_dir, methods)
    if not units:
        units = build_full_analysis_units(results_dir, methods, option_sets)

    mapping_data = read_tsv(results_dir / "mapping_summary.tsv")
    fc_data = read_tsv(results_dir / "featurecounts_summary.tsv")
    filtering_data = read_tsv(results_dir / "filtering_summary.tsv")
    de_summary_data = read_tsv(results_dir / "de_summary.tsv")
    redundancy_data = read_tsv(results_dir / "redundancy_summary.tsv")

    metrics = _build_unit_metrics(
        units, mapping_data, fc_data, filtering_data, de_summary_data, redundancy_data
    )

    if cfg.get("selection"):
        selected_primary, selected_comparison, primary_reason, comp_reason = (
            _apply_selection_policy(cfg, results_dir, units, metrics)
        )
        selected_count = selected_comparison
        selected_vis = selected_primary
        count_reason = comp_reason
        vis_reason = primary_reason
    else:
        selected_count, count_reason = _select_count_comparison_branch(
            cfg, results_dir, units, metrics,
        )
        selected_vis, vis_reason = _select_visualisation_branch(
            cfg, results_dir, units,
        )
        if selected_count is not None:
            write_selected_analysis(
                results_dir, selected_count,
                reason=f"compat_from_count_comparison:{count_reason}",
            )

    report_dir = results_dir / "reports"
    ensure_dirs(report_dir)

    lines: List[str] = []
    lines.append("# Method and Parameter Comparison")
    lines.append("")
    lines.append("## Trimming Configuration")
    lines.append(f"- Primary trim method: `{trim_cfg['primary_method']}`")
    lines.append(f"- Trim comparison enabled: `{trim_cfg['compare_methods']}`")
    lines.append(f"- Effective trim methods: `{', '.join(trim_cfg['effective_methods'])}`")
    if trim_cfg["non_selected_methods"]:
        lines.append(
            "- Intentionally not run (trim comparison disabled/excluded): "
            f"`{', '.join(trim_cfg['non_selected_methods'])}`"
        )
    lines.append("")
    lines.append("## Selection Outcomes")
    if cfg.get("selection"):
        lines.append(
            f"- Selected primary analysis branch: "
            f"`{selected_vis.label if selected_vis else 'N/A'}` ({vis_reason})"
        )
    lines.append(
        f"- Selected count-comparison branch: "
        f"`{selected_count.label if selected_count else 'N/A'}` ({count_reason})"
    )
    lines.append(
        f"- Selected visualisation branch: "
        f"`{selected_vis.label if selected_vis else 'N/A'}` ({vis_reason})"
    )
    lines.append("- Count comparison focuses on featureCounts option interpretation.")
    lines.append("- Visualisation branch is used for selected-only BigWig generation.")
    lines.append("")

    lines.append("## Branch Summary")
    lines.append(
        "| Branch | Role | Mean uniquely mapped % | Mean assigned % | "
        "Mean unassigned_no_features % | Mean unassigned_multimapping % | "
        "Genes before filtering | Genes after filtering | Genes removed | "
        "DE tested genes | DE significant genes | Redundancy | Selection |"
    )
    lines.append(
        "|--------|------|------------------------|------------------|"
        "-----------------------------|-------------------------------|"
        "------------------------|-----------------------|--------------|"
        "-----------------|---------------------|-----------|-----------|"
    )

    redundant_lines: List[str] = []
    for au in sorted(units):
        m = metrics[au]
        uniq = _mean(m["uniq_vals"])
        assigned = _mean(m["assigned_vals"])
        no_feat = _mean(m["unassigned_no_features_vals"])
        multi = _mean(m["unassigned_multimapping_vals"])
        genes_in = _mean(m["genes_in_vals"])
        genes_out = _mean(m["genes_out_vals"])
        genes_removed = max(0.0, genes_in - genes_out)

        if m["is_redundant"]:
            role = "redundant"
        elif selected_vis is not None and au == selected_vis:
            role = "visualisation"
        elif au.method == "none" and au.mapper == "star" and au.mapper_option_set == "default":
            role = "count_comparison"
        else:
            role = "alternative"

        selection = []
        if selected_count is not None and au == selected_count:
            selection.append("selected_count")
        if selected_vis is not None and au == selected_vis:
            selection.append("selected_visualisation")
        selection_label = ",".join(selection) if selection else "-"

        redundancy_label = "no"
        if m["is_redundant"]:
            redundancy_label = f"duplicate_of:{m['canonical_label']}"
            redundant_lines.append(f"- `{au.label}` duplicates `{m['canonical_label']}`")

        lines.append(
            f"| `{au.label}` | {role} | {uniq:.2f} | {assigned:.2f} | {no_feat:.2f} | {multi:.2f} | "
            f"{int(genes_in)} | {int(genes_out)} | {int(genes_removed)} | "
            f"{int(m['de_tested'])} | {int(m['de_sig'])} | {redundancy_label} | {selection_label} |"
        )

    lines.append("")
    lines.append("## Redundant Branches")
    if redundant_lines:
        lines.extend(redundant_lines)
    else:
        lines.append("- No redundant branches detected.")

    lines.append("")
    lines.append("## Why These Branches Were Selected")
    lines.append(
        f"- Count-comparison selection uses valid DE output, assigned percentage, genes retained, and non-redundancy; chosen: `{selected_count.label if selected_count else 'N/A'}`."
    )
    lines.append(
        f"- Visualisation selection prioritises `none/star/strict_unique/default` for unique-mapper BigWig generation; chosen: `{selected_vis.label if selected_vis else 'N/A'}`."
    )

    lines.append("")
    lines.append("## Trimming Impact Summary")
    if trim_cfg["compare_methods"] and len(trim_cfg["effective_methods"]) > 1:
        by_method: Dict[str, Dict[str, float]] = defaultdict(
            lambda: {"n": 0.0, "uniq": 0.0, "assigned": 0.0, "de_sig": 0.0}
        )
        for au in units:
            m = metrics[au]
            by_method[au.method]["n"] += 1
            by_method[au.method]["uniq"] += _mean(m["uniq_vals"])
            by_method[au.method]["assigned"] += _mean(m["assigned_vals"])
            by_method[au.method]["de_sig"] += float(m["de_sig"])
        lines.append("| Trim Method | Mean uniquely mapped % | Mean assigned % | Mean DE significant genes |")
        lines.append("|-------------|------------------------|------------------|---------------------------|")
        for method in trim_cfg["effective_methods"]:
            stats = by_method.get(method, {"n": 0.0, "uniq": 0.0, "assigned": 0.0, "de_sig": 0.0})
            n = stats["n"] if stats["n"] > 0 else 1.0
            lines.append(
                f"| {method} | {stats['uniq']/n:.2f} | {stats['assigned']/n:.2f} | {stats['de_sig']/n:.2f} |"
            )
        lines.append(
            "Trim comparison was enabled; table above summarizes material differences in QC/mapping/counting/DE outcomes."
        )
    else:
        lines.append(
            "Trim comparison was disabled for this run; non-primary trim methods were intentionally not run."
        )

    report_path = report_dir / "method_comparison.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
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
