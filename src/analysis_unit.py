"""
Shared analysis-unit helpers for the RNA-seq pipeline.

Centralises mapping-unit discovery, path resolution, and the AnalysisUnit
abstraction so that ``count_option_set`` is a first-class dimension
everywhere: counting, filtering, DE, comparison, BigWig, and reporting.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# AnalysisUnit dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True, order=True)
class AnalysisUnit:
    """Immutable 4-tuple identifying one analysis branch."""

    method: str
    mapper: str
    mapper_option_set: str
    count_option_set: str

    @property
    def label(self) -> str:
        return f"{self.method}/{self.mapper}/{self.mapper_option_set}/{self.count_option_set}"

    @property
    def mapping_key(self) -> Tuple[str, str, str]:
        return (self.method, self.mapper, self.mapper_option_set)

    def as_dict(self) -> Dict[str, str]:
        return asdict(self)


def unit_label(unit: Dict[str, str]) -> str:
    """Human-readable label for a dict-style analysis unit."""
    method = unit.get("method", "")
    mapper = unit.get("mapper", "star")
    mapper_opt = unit.get("mapper_option_set", "default")
    count_opt = unit.get("count_option_set", "")
    parts = [method, mapper, mapper_opt]
    if count_opt:
        parts.append(count_opt)
    return "/".join(parts)


def mapping_unit_label(unit: Dict[str, str]) -> str:
    """Label for the 3-tuple mapping level (no count_option_set)."""
    method = unit.get("method", "")
    mapper = unit.get("mapper", "star")
    mapper_opt = unit.get("mapper_option_set", "default")
    return f"{method}/{mapper}/{mapper_opt}"


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def resolve_fc_dir(
    results_dir: Path,
    method: str,
    mapper: str,
    mapper_opt: str,
) -> Path:
    """Return the featureCounts output directory for a mapping unit."""
    return results_dir / method / "featurecounts" / mapper / mapper_opt


def resolve_count_matrix(
    results_dir: Path,
    method: str,
    mapper: str,
    mapper_opt: str,
    count_opt: str,
) -> Optional[Path]:
    """Return the filtered (clean) count matrix path, with legacy fallback.

    Returns None if no file can be found.
    """
    fc_dir = resolve_fc_dir(results_dir, method, mapper, mapper_opt)
    primary = fc_dir / f"clean_matrix_{count_opt}.tsv"
    if primary.exists():
        return primary

    legacy = results_dir / method / "featurecounts" / f"clean_matrix_{count_opt}.tsv"
    if legacy.exists():
        return legacy

    return None


def resolve_de_dir(
    results_dir: Path,
    method: str,
    mapper: str,
    mapper_opt: str,
    count_opt: str,
) -> Path:
    """Return the DE output directory for a full analysis unit."""
    return results_dir / method / "deseq2" / mapper / mapper_opt / count_opt


def resolve_de_base_legacy(
    results_dir: Path,
    method: str,
    mapper: str,
    mapper_opt: str,
) -> Path:
    """Legacy DE base (no count_option_set). Used for backwards compatibility."""
    nested = results_dir / method / "deseq2" / mapper / mapper_opt
    if nested.exists():
        return nested
    return results_dir / method / "deseq2"


# ---------------------------------------------------------------------------
# Mapping-unit builders (with BAM references)
# ---------------------------------------------------------------------------

def build_mapping_units_with_bams(
    results_dir: Path,
    samples: list,
    methods: List[str],
    *,
    use_filtered_bam: bool = False,
    prefer_filtered: bool = False,
) -> List[Dict[str, Any]]:
    """Collect mapping units with per-sample BAM paths.

    Unified version used by steps that need BAM files (counting, BigWig).

    Parameters
    ----------
    use_filtered_bam : bool
        When True and ``filtered_bam_path`` exists in the mapping summary,
        substitute it for ``bam_path``.  Used by featureCounts.
    prefer_filtered : bool
        When True and ``filtered_bam_path`` exists on disk, prefer it.
        Used by BigWig generation.
    """
    mapping_summary = results_dir / "mapping_summary.tsv"
    units: Dict[Tuple[str, str, str], Dict[str, Any]] = {}

    if mapping_summary.exists():
        with open(mapping_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("method") or row.get("trim_method") or ""
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                sample = row.get("sample", "")
                bam_path = row.get("bam_path", "")

                if use_filtered_bam and row.get("filtered_bam_path"):
                    bam_path = row.get("filtered_bam_path", "")
                elif prefer_filtered:
                    filtered = row.get("filtered_bam_path", "")
                    if filtered and Path(filtered).exists():
                        bam_path = filtered

                if method not in methods or not sample or not bam_path:
                    continue

                key = (method, mapper, mapper_opt)
                if key not in units:
                    units[key] = {
                        "method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "sample_to_bam": {},
                    }
                units[key]["sample_to_bam"][sample] = Path(bam_path)

    if not units:
        for method in methods:
            sample_to_bam: Dict[str, Path] = {}
            star_dir = results_dir / method / "star"
            for s in samples:
                bam = star_dir / f"{s.sample_name}_Aligned.sortedByCoord.out.bam"
                if bam.exists():
                    sample_to_bam[s.sample_name] = bam
            if sample_to_bam:
                units[(method, "star", "default")] = {
                    "method": method,
                    "mapper": "star",
                    "mapper_option_set": "default",
                    "sample_to_bam": sample_to_bam,
                }

    return [units[k] for k in sorted(units.keys())]


# ---------------------------------------------------------------------------
# Mapping-unit builder (lightweight, no BAM references)
# ---------------------------------------------------------------------------

def build_mapping_units(
    results_dir: Path,
    methods: List[str],
) -> List[Dict[str, str]]:
    """Infer 3-tuple mapping units from summary files.

    Used by steps that do not need BAM paths (DE, filter, compare, report).
    """
    units: Set[Tuple[str, str, str]] = set()
    mapping_summary = results_dir / "mapping_summary.tsv"

    if mapping_summary.exists():
        with open(mapping_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("method") or row.get("trim_method") or ""
                if method not in methods:
                    continue
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                units.add((method, mapper, mapper_opt))

    if not units:
        for method in methods:
            units.add((method, "star", "default"))

    return [
        {"method": m, "mapper": p, "mapper_option_set": o}
        for (m, p, o) in sorted(units)
    ]


def build_full_analysis_units(
    results_dir: Path,
    methods: List[str],
    option_sets: Dict[str, Any],
) -> List[AnalysisUnit]:
    """Build the complete 4-tuple analysis units.

    Cross-product of mapping units x featureCounts option sets.
    """
    mapping_units = build_mapping_units(results_dir, methods)
    units: List[AnalysisUnit] = []
    for mu in mapping_units:
        for opt_name in option_sets:
            units.append(AnalysisUnit(
                method=mu["method"],
                mapper=mu["mapper"],
                mapper_option_set=mu["mapper_option_set"],
                count_option_set=opt_name,
            ))
    return sorted(units)


def infer_analysis_units_from_de_summary(
    results_dir: Path,
    methods: List[str],
) -> List[AnalysisUnit]:
    """Infer 4-tuple analysis units from de_summary.tsv."""
    units: Set[Tuple[str, str, str, str]] = set()
    de_summary = results_dir / "de_summary.tsv"
    if de_summary.exists():
        with open(de_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("method") or row.get("trim_method") or ""
                if method not in methods:
                    continue
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                count_opt = row.get("count_option_set", "default")
                units.add((method, mapper, mapper_opt, count_opt))

    if not units:
        mapping = build_mapping_units(results_dir, methods)
        for mu in mapping:
            units.add((mu["method"], mu["mapper"], mu["mapper_option_set"], "default"))

    return sorted(
        AnalysisUnit(m, p, o, c)
        for (m, p, o, c) in units
    )


# ---------------------------------------------------------------------------
# Selected-branch helpers
# ---------------------------------------------------------------------------

def read_selected_analysis(results_dir: Path) -> Optional[AnalysisUnit]:
    """Read the selected analysis branch from selected_analysis.tsv."""
    path = results_dir / "selected_analysis.tsv"
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return AnalysisUnit(
                method=row.get("trim_method", ""),
                mapper=row.get("mapper", "star"),
                mapper_option_set=row.get("mapper_option_set", "default"),
                count_option_set=row.get("count_option_set", "default"),
            )
    return None


def write_selected_analysis(
    results_dir: Path,
    unit: AnalysisUnit,
    reason: str = "",
) -> Path:
    """Write the selected analysis branch to selected_analysis.tsv."""
    path = results_dir / "selected_analysis.tsv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["trim_method", "mapper", "mapper_option_set",
                         "count_option_set", "selection_reason"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({
            "trim_method": unit.method,
            "mapper": unit.mapper,
            "mapper_option_set": unit.mapper_option_set,
            "count_option_set": unit.count_option_set,
            "selection_reason": reason,
        })
    logger.info("Selected analysis branch -> %s", path)
    return path
