#!/usr/bin/env python3
"""
Step 08 -- Filter count matrix.

Implements the lecture-style filter: drop genes where the maximum count
across all samples is below a configurable threshold (default: number of
samples -- rule of thumb ~2 reads per sample).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv
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


from src.analysis_unit import build_mapping_units


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
    methods = methods_override or get_enabled_methods(cfg)

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
