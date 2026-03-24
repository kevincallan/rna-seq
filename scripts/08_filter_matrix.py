#!/usr/bin/env python3
"""
Step 08 -- Filter count matrix.

Implements the lecture-style filter: drop genes where the maximum count
across all samples is below a configurable threshold (default: number of
samples -- rule of thumb ~2 reads per sample).
"""

from __future__ import annotations

import hashlib
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv
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

logger = logging.getLogger(__name__)


from src.analysis_unit import build_mapping_units


def _matrix_fingerprint(path: Path) -> str:
    """Stable sha256 fingerprint for a filtered matrix file."""
    hasher = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def filter_matrix(
    input_path: Path,
    output_path: Path,
    threshold: int,
) -> Dict[str, int]:
    """Filter a gene count matrix, keeping genes with max count >= threshold.

    Parameters
    ----------
    input_path : Path
        Input TSV with Geneid + sample columns.
    output_path : Path
        Output filtered TSV.
    threshold : int
        Minimum max-count across samples to keep a gene.

    Returns
    -------
    dict
        Statistics: genes_in, genes_out, genes_removed.
    """
    genes_in = 0
    genes_out = 0

    with open(input_path, "r", encoding="utf-8") as fin, \
         open(output_path, "w", encoding="utf-8") as fout:

        header = fin.readline()
        fout.write(header)

        for line in fin:
            genes_in += 1
            parts = line.strip().split("\t")
            gene_id = parts[0]
            try:
                counts = [int(float(v)) for v in parts[1:]]
            except ValueError:
                # Skip lines with non-numeric values
                continue

            max_count = max(counts) if counts else 0

            if max_count >= threshold:
                fout.write(line)
                genes_out += 1

    stats = {
        "genes_in": genes_in,
        "genes_out": genes_out,
        "genes_removed": genes_in - genes_out,
    }
    logger.info(
        "Filtered %s: %d -> %d genes (removed %d, threshold=%d)",
        input_path.name, genes_in, genes_out, stats["genes_removed"], threshold,
    )
    return stats


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 08."""
    logger.info("=" * 60)
    logger.info("STEP 08: Filter count matrices")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = get_effective_trim_methods(cfg, methods_override)
    trim_cfg = get_trim_config_summary(cfg)
    logger.info(
        "Trim config: primary=%s compare_methods=%s effective=%s",
        trim_cfg["primary_method"],
        trim_cfg["compare_methods"],
        ",".join(methods),
    )

    # Determine threshold
    filter_cfg = cfg.get("filtering", {})
    min_max = filter_cfg.get("min_max_count", "auto")
    if min_max == "auto":
        threshold = len(samples)
    else:
        threshold = int(min_max)

    logger.info("Filter threshold: max count >= %d", threshold)

    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})
    all_stats: List[Dict[str, Any]] = []
    filtered_matrices: List[Dict[str, Any]] = []

    mapping_units = build_mapping_units(results_dir, methods)
    for unit in mapping_units:
        method = unit["method"]
        mapper = unit["mapper"]
        mapper_opt = unit["mapper_option_set"]
        fc_dir = results_dir / method / "featurecounts" / mapper / mapper_opt

        for opt_name in option_sets:
            input_matrix = fc_dir / f"count_matrix_{opt_name}.tsv"
            if not input_matrix.exists():
                logger.warning("Count matrix not found: %s", input_matrix)
                continue

            output_matrix = fc_dir / f"clean_matrix_{opt_name}.tsv"
            stats = filter_matrix(input_matrix, output_matrix, threshold)
            stats["trim_method"] = method
            stats["mapper"] = mapper
            stats["mapper_option_set"] = mapper_opt
            stats["option_set"] = opt_name
            all_stats.append(stats)
            filtered_matrices.append({
                "trim_method": method,
                "mapper": mapper,
                "mapper_option_set": mapper_opt,
                "count_option_set": opt_name,
                "matrix_path": str(output_matrix),
                "fingerprint": _matrix_fingerprint(output_matrix),
            })

    # Write filtering summary
    summary_path = results_dir / "filtering_summary.tsv"
    if all_stats:
        with open(summary_path, "w", encoding="utf-8") as fh:
            fh.write(
                "trim_method\tmapper\tmapper_option_set\t"
                "option_set\tgenes_in\tgenes_out\tgenes_removed\n"
            )
            for s in all_stats:
                fh.write(
                    f"{s['trim_method']}\t{s['mapper']}\t"
                    f"{s['mapper_option_set']}\t{s['option_set']}\t"
                    f"{s['genes_in']}\t{s['genes_out']}\t{s['genes_removed']}\n"
                )
        logger.info("Filtering summary -> %s", summary_path)

    # Write redundancy summary from clean-matrix fingerprints
    redundancy_path = results_dir / "redundancy_summary.tsv"
    canonical_by_fp: Dict[str, Tuple[str, str, str, str]] = {}
    for row in filtered_matrices:
        key = (
            row["trim_method"],
            row["mapper"],
            row["mapper_option_set"],
            row["count_option_set"],
        )
        fp = row["fingerprint"]
        if fp not in canonical_by_fp or key < canonical_by_fp[fp]:
            canonical_by_fp[fp] = key

    with open(redundancy_path, "w", encoding="utf-8") as fh:
        fh.write(
            "trim_method\tmapper\tmapper_option_set\tcount_option_set\t"
            "fingerprint\tis_redundant\tcanonical_trim_method\tcanonical_mapper\t"
            "canonical_mapper_option_set\tcanonical_count_option_set\tredundancy_reason\n"
        )
        for row in sorted(
            filtered_matrices,
            key=lambda r: (
                r["trim_method"], r["mapper"], r["mapper_option_set"], r["count_option_set"]
            ),
        ):
            key = (
                row["trim_method"],
                row["mapper"],
                row["mapper_option_set"],
                row["count_option_set"],
            )
            canonical = canonical_by_fp[row["fingerprint"]]
            is_redundant = key != canonical
            reason = "identical_clean_matrix" if is_redundant else "canonical"
            fh.write(
                f"{row['trim_method']}\t{row['mapper']}\t{row['mapper_option_set']}\t"
                f"{row['count_option_set']}\t{row['fingerprint']}\t"
                f"{'true' if is_redundant else 'false'}\t"
                f"{canonical[0]}\t{canonical[1]}\t{canonical[2]}\t{canonical[3]}\t{reason}\n"
            )
    logger.info("Redundancy summary -> %s", redundancy_path)

    logger.info("STEP 08 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 08: Filter count matrix")
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
