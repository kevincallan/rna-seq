"""
Metadata parsing and sample design table generation.

Reads an SRA-style metadata CSV, applies subset filters, builds a
deterministic sample table with replicate numbering, and creates
FASTQ symlinks using the convention ``condition_rep[_mate].fastq.gz``.
"""

from __future__ import annotations

import csv
import logging
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.utils import safe_name

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

class Sample:
    """Represents a single biological sample (one Run ID)."""

    __slots__ = (
        "run_id", "condition", "replicate", "sample_name",
        "layout", "r1", "r2",
    )

    def __init__(
        self,
        run_id: str,
        condition: str,
        replicate: int,
        layout: str = "paired",
        r1: str = "",
        r2: str = "",
    ):
        self.run_id = run_id
        self.condition = condition
        self.replicate = replicate
        self.layout = layout
        self.r1 = r1
        self.r2 = r2
        self.sample_name = f"{condition}_{replicate}"

    def __repr__(self) -> str:
        return (
            f"Sample({self.sample_name!r}, run={self.run_id}, "
            f"cond={self.condition}, rep={self.replicate}, layout={self.layout})"
        )


# ---------------------------------------------------------------------------
# Metadata parsing
# ---------------------------------------------------------------------------

def parse_metadata(cfg: Dict[str, Any]) -> List[Dict[str, str]]:
    """Read metadata CSV and return list of row dicts.

    Parameters
    ----------
    cfg : dict
        Full pipeline config.

    Returns
    -------
    list of dict
        One dict per CSV row.
    """
    csv_path = cfg["data"]["metadata_csv"]
    logger.info("Reading metadata from %s", csv_path)

    rows: List[Dict[str, str]] = []
    with open(csv_path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(dict(row))

    logger.info("Read %d rows from metadata CSV", len(rows))
    return rows


def validate_metadata_columns(
    rows: List[Dict[str, str]],
    cfg: Dict[str, Any],
) -> None:
    """Validate that configured column names actually exist in the CSV.

    Raises ``ValueError`` for missing required columns (run_id_col,
    condition_cols).  Logs a warning for missing subset-filter columns.
    """
    if not rows:
        raise ValueError("Metadata CSV is empty (0 rows)")

    csv_cols = set(rows[0].keys())
    col_map = cfg.get("column_mapping", {})

    # --- run_id_col ---
    run_col = col_map.get("run_id_col", "Run")
    if run_col not in csv_cols:
        raise ValueError(
            f"run_id_col '{run_col}' not found in metadata CSV. "
            f"Available columns: {sorted(csv_cols)}"
        )

    # --- condition_cols ---
    for cc in col_map.get("condition_cols", []):
        if cc not in csv_cols:
            raise ValueError(
                f"condition_col '{cc}' not found in metadata CSV. "
                f"Available columns: {sorted(csv_cols)}"
            )

    # --- subset_filters (warn only) ---
    subset_key = cfg.get("active_subset", "default")
    all_filters = cfg.get("subset_filters", {})

    if subset_key not in all_filters and all_filters:
        logger.warning(
            "active_subset '%s' not found in subset_filters (available: %s)",
            subset_key, sorted(all_filters.keys()),
        )

    filters = all_filters.get(subset_key, {})
    for fc in filters:
        if fc not in csv_cols:
            logger.warning(
                "Subset filter column '%s' not found in metadata CSV. "
                "Available columns: %s", fc, sorted(csv_cols),
            )


def apply_subset_filters(
    rows: List[Dict[str, str]],
    cfg: Dict[str, Any],
    subset_name: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Filter metadata rows based on the active (or named) subset.

    Parameters
    ----------
    rows : list of dict
        Raw metadata rows.
    cfg : dict
        Full config.
    subset_name : str, optional
        Override the subset key from config.

    Returns
    -------
    list of dict
        Filtered rows.
    """
    subset_key = subset_name or cfg.get("active_subset", "default")
    filters = cfg.get("subset_filters", {}).get(subset_key)

    if not filters:
        logger.warning("No subset filters found for '%s'; using all rows", subset_key)
        return rows

    logger.info("Applying subset filter '%s': %s", subset_key, filters)

    filtered: List[Dict[str, str]] = []
    for row in rows:
        match = True
        for col, expected in filters.items():
            actual = row.get(col, "").strip()
            if actual != expected:
                match = False
                break
        if match:
            filtered.append(row)

    logger.info(
        "Subset '%s': %d / %d rows passed filter",
        subset_key, len(filtered), len(rows),
    )

    if not filtered:
        raise ValueError(
            f"Subset filter '{subset_key}' matched 0 rows. "
            "Check column names and filter values in config.yaml."
        )

    return filtered


# ---------------------------------------------------------------------------
# Condition string building
# ---------------------------------------------------------------------------

def build_condition(row: Dict[str, str], cfg: Dict[str, Any]) -> str:
    """Build a condition string from one metadata row.

    Uses ``column_mapping.condition_cols`` (list) and optionally maps
    values through ``column_mapping.condition_map``.

    Raises ``ValueError`` if the resulting condition is empty.
    """
    col_map = cfg["column_mapping"]
    cond_cols: List[str] = col_map["condition_cols"]
    value_map: Dict[str, str] = col_map.get("condition_map", {})

    parts: List[str] = []
    for col in cond_cols:
        raw = row.get(col, "").strip()
        mapped = value_map.get(raw, raw)
        parts.append(safe_name(mapped))

    result = "_".join(p for p in parts if p)
    if not result:
        run_id = row.get(col_map.get("run_id_col", "Run"), "?")
        raise ValueError(
            f"Empty condition for run '{run_id}'. "
            "Check condition_cols and condition_map in config."
        )
    return result


# ---------------------------------------------------------------------------
# Design table (samples list)
# ---------------------------------------------------------------------------

def _resolve_fastq_r1(run_id: str, fastq_dir: str) -> str:
    """Find the R1 FASTQ, supporting both ``{run_id}_1.fastq.gz`` and ``{run_id}.fastq.gz``."""
    candidate = Path(fastq_dir) / f"{run_id}_1.fastq.gz"
    if candidate.exists():
        return str(candidate)
    bare = Path(fastq_dir) / f"{run_id}.fastq.gz"
    if bare.exists():
        return str(bare)
    return str(candidate)


def build_design_table(
    rows: List[Dict[str, str]],
    cfg: Dict[str, Any],
) -> List[Sample]:
    """Build a list of Sample objects with deterministic replicate numbering.

    Parameters
    ----------
    rows : list of dict
        Filtered metadata rows.
    cfg : dict
        Full config.

    Returns
    -------
    list of Sample
        Sorted by condition then replicate.
    """
    col_map = cfg["column_mapping"]
    run_col = col_map["run_id_col"]
    layout = cfg["data"].get("layout", "paired")
    fastq_dir = cfg["data"]["fastq_dir"]
    strategy = col_map.get("replicate_strategy", "sorted_run_id")
    strict = cfg.get("_strict", False)

    # Check for duplicate run IDs
    run_ids = [row[run_col].strip() for row in rows]
    seen_runs: Dict[str, int] = {}
    for rid in run_ids:
        seen_runs[rid] = seen_runs.get(rid, 0) + 1
    duplicates = {k: v for k, v in seen_runs.items() if v > 1}
    if duplicates:
        raise ValueError(
            f"Duplicate run IDs in filtered metadata: {duplicates}. "
            "Each run ID must appear only once."
        )

    # Group runs by condition
    cond_runs: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        cond = build_condition(row, cfg)
        cond_runs[cond].append(row)

    samples: List[Sample] = []

    for cond in sorted(cond_runs.keys()):
        group = cond_runs[cond]

        # Deterministic ordering
        if strategy == "sorted_run_id":
            group.sort(key=lambda r: r[run_col])
        else:
            group.sort(key=lambda r: r[run_col])

        for rep_idx, row in enumerate(group, start=1):
            run_id = row[run_col].strip()

            sample_layout = _detect_layout(run_id, fastq_dir, layout)

            r1 = _resolve_fastq_r1(run_id, fastq_dir)
            r2 = str(Path(fastq_dir) / f"{run_id}_2.fastq.gz") if sample_layout == "paired" else ""

            # FASTQ existence check
            missing: List[str] = []
            if not Path(r1).exists():
                missing.append(r1)
            if sample_layout == "paired" and r2 and not Path(r2).exists():
                missing.append(r2)
            if missing:
                msg = f"FASTQ not found for {run_id}: {missing}"
                if strict:
                    raise FileNotFoundError(msg)
                logger.warning(msg)

            s = Sample(
                run_id=run_id,
                condition=cond,
                replicate=rep_idx,
                layout=sample_layout,
                r1=r1,
                r2=r2,
            )
            samples.append(s)
            logger.info("  %s -> %s", run_id, s.sample_name)

    # Check for duplicate derived sample names
    name_counts: Dict[str, int] = {}
    for s in samples:
        name_counts[s.sample_name] = name_counts.get(s.sample_name, 0) + 1
    dup_names = {k: v for k, v in name_counts.items() if v > 1}
    if dup_names:
        raise ValueError(
            f"Duplicate derived sample names: {dup_names}. "
            "This usually means condition_map produces the same label for "
            "different raw values. Check condition_cols and condition_map."
        )

    logger.info("Design table: %d samples across %d conditions",
                len(samples), len(set(s.condition for s in samples)))
    return samples


def _detect_layout(
    run_id: str, fastq_dir: str, default: str
) -> str:
    """Detect paired vs single from file existence."""
    default = default.strip().lower()
    if default in ("paired", "single"):
        return default

    # Auto-detect
    r2 = Path(fastq_dir) / f"{run_id}_2.fastq.gz"
    return "paired" if r2.exists() else "single"


# ---------------------------------------------------------------------------
# Symlink creation
# ---------------------------------------------------------------------------

def create_symlinks(
    samples: List[Sample],
    work_dir: Path,
) -> None:
    """Create FASTQ symlinks with standardised names.

    Creates symlinks in ``work_dir/fastq_links/`` using the convention:
    - paired: ``condition_rep_1.fastq.gz``, ``condition_rep_2.fastq.gz``
    - single: ``condition_rep.fastq.gz``
    """
    link_dir = work_dir / "fastq_links"
    link_dir.mkdir(parents=True, exist_ok=True)

    for s in samples:
        if s.layout == "paired":
            targets = [
                (s.r1, link_dir / f"{s.sample_name}_1.fastq.gz"),
                (s.r2, link_dir / f"{s.sample_name}_2.fastq.gz"),
            ]
        else:
            targets = [
                (s.r1, link_dir / f"{s.sample_name}.fastq.gz"),
            ]

        for src, dst in targets:
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(os.path.abspath(src), dst)
            logger.debug("Symlink: %s -> %s", dst.name, src)

    logger.info("Created %d symlinks in %s",
                sum(2 if s.layout == "paired" else 1 for s in samples),
                link_dir)


# ---------------------------------------------------------------------------
# Write samples.tsv
# ---------------------------------------------------------------------------

def write_samples_tsv(samples: List[Sample], path: Path) -> None:
    """Write the design table to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("sample\tcondition\trun_id\treplicate\tlayout\n")
        for s in samples:
            fh.write(
                f"{s.sample_name}\t{s.condition}\t{s.run_id}"
                f"\t{s.replicate}\t{s.layout}\n"
            )
    logger.info("Wrote samples.tsv -> %s (%d samples)", path, len(samples))


# ---------------------------------------------------------------------------
# Write DESeq2 sample description
# ---------------------------------------------------------------------------

def write_sample_description(
    samples: List[Sample],
    path: Path,
) -> None:
    """Write sample_description.txt for DESeq2.

    Format:
        name<TAB>condition
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("name\tcondition\n")
        for s in samples:
            fh.write(f"{s.sample_name}\t{s.condition}\n")
    logger.info("Wrote sample_description.txt -> %s", path)


# ---------------------------------------------------------------------------
# Read back samples.tsv
# ---------------------------------------------------------------------------

def read_samples_tsv(path: Path) -> List[Sample]:
    """Read a samples.tsv back into Sample objects."""
    samples: List[Sample] = []
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            s = Sample(
                run_id=row["run_id"],
                condition=row["condition"],
                replicate=int(row["replicate"]),
                layout=row.get("layout", "paired"),
            )
            samples.append(s)
    return samples
