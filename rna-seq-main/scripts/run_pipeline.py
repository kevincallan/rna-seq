#!/usr/bin/env python3
"""
RNA-seq Pipeline Orchestrator
==============================

Runs the full pipeline or individual steps.  Configuration comes from
``config/config.yaml``; dataset and species can be overridden on the
command line so only one argument changes on exam day.

Usage examples (JupyterHub)
---------------------------
# Full pipeline -- exam day command:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run

# Explicit paths:
./py scripts/run_pipeline.py --data-root /data --dataset GSE48519 --species mouse run

# Specific steps only:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --steps 0 1 2 5

# Override trimming methods:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --methods none cutadapt

# Use a specific subset filter:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --subset day3

# Single step:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --steps 5
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

# Allow running from repository root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import (
    get_enabled_methods,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger("pipeline")


# ---------------------------------------------------------------------------
# Species -> reference path mapping
# ---------------------------------------------------------------------------

SPECIES_MAP: Dict[str, Dict[str, str]] = {
    "mouse": {
        "index_dir": "mm39",
        "genome_dir": "/data/genomes/mouse/GRCm39",
    },
    "human": {
        "index_dir": "hg38",
        "genome_dir": "/data/genomes/human",
    },
}


# ---------------------------------------------------------------------------
# Runtime & data validation (runs before heavy compute)
# ---------------------------------------------------------------------------

def _check_runtime() -> None:
    """Print interpreter info, verify packages, fail fast if anything missing."""
    print("=" * 60)
    print("Runtime check")
    print("=" * 60)
    print(f"Interpreter: {sys.executable}")
    print(f"Python:      {sys.version.split(chr(10))[0]}")

    if not sys.executable.startswith("/opt/jupyterhub/"):
        print(f"WARNING: Not running under /opt/jupyterhub/bin/python3")
        print(f"         Actual: {sys.executable}")
        print(f"         Use:    ./py scripts/run_pipeline.py ...")

    missing: List[str] = []
    # HTSeq only required when featurecounts.backend is "htseq"; default is featurecounts
    for pkg, import_name in [
        ("pydeseq2", "pydeseq2"),
        ("pandas", "pandas"),
        ("scikit-learn", "sklearn"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("matplotlib", "matplotlib"),
        ("pyyaml", "yaml"),
        ("pysam", "pysam"),
    ]:
        try:
            mod = __import__(import_name)
            ver = getattr(mod, "__version__", "?")
            print(f"  {pkg:15s} {ver}")
        except ImportError:
            missing.append(pkg)
            print(f"  {pkg:15s} MISSING")

    if missing:
        sys.exit(f"FATAL: Missing packages: {', '.join(missing)}")

    print("=" * 60)


def _find_star_index(genome_dir: Path) -> str:
    """Auto-discover STAR index inside a genome directory.

    Looks for a directory containing Genome, SA, SAindex (STAR index files).
    Tries: genome_dir/STAR, genome_dir/star, genome_dir itself.
    """
    for candidate in [genome_dir / "STAR", genome_dir / "star", genome_dir]:
        if (candidate / "Genome").exists() or (candidate / "SA").exists():
            return str(candidate)
    return str(genome_dir / "STAR")


def _find_gtf(genome_dir: Path, build: str) -> str:
    """Auto-discover GTF file inside a genome directory.

    Preference order: gene_names GTF > chr GTF > any GTF > gzipped versions.
    """
    if build != "*":
        for pattern in [f"{build}.gtf", f"{build}*.gtf"]:
            matches = sorted(genome_dir.glob(pattern))
            if matches:
                return str(matches[0])

    # Prefer gene_names variants, then chr variants, then any
    for pattern in ["*gene_names*.gtf.gz", "*gene_names*.gtf",
                    "*chr*.gtf.gz", "*chr*.gtf",
                    "*.gtf", "*.gtf.gz"]:
        matches = sorted(genome_dir.glob(pattern))
        if matches:
            return str(matches[-1])  # latest version (highest number)
    return str(genome_dir / f"{build}.gtf")


def _find_fasta(genome_dir: Path, build: str) -> str:
    """Auto-discover genome FASTA inside a genome directory.

    Prefers uncompressed, then gzipped. Prefers primary_assembly.
    """
    if build != "*":
        for pattern in [f"{build}.fa", f"{build}.fasta",
                        f"{build}*.fa", f"{build}*.fasta"]:
            matches = sorted(genome_dir.glob(pattern))
            if matches:
                return str(matches[0])

    for pattern in ["*primary_assembly*.fa", "*primary_assembly*.fa.gz",
                    "*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz"]:
        matches = sorted(genome_dir.glob(pattern))
        if matches:
            return str(matches[0])
    return str(genome_dir / f"{build}.fa")


def _apply_data_pointer(args: argparse.Namespace, cfg: Dict[str, Any]) -> None:
    """Override config paths based on --dataset / --species CLI flags."""
    if args.dataset:
        dataset_dir = Path(args.data_root) / args.dataset
        cfg["data"]["fastq_dir"] = str(dataset_dir)
        cfg["data"]["metadata_csv"] = str(dataset_dir / "metadata.csv")
        cfg["_data_root"] = args.data_root
        cfg["_dataset"] = args.dataset

    if args.metadata:
        cfg["data"]["metadata_csv"] = args.metadata

    if args.species:
        if args.species not in SPECIES_MAP:
            sys.exit(
                f"FATAL: Unknown species '{args.species}'. "
                f"Valid: {', '.join(sorted(SPECIES_MAP))}"
            )
        sp = SPECIES_MAP[args.species]
        idx_base = Path(args.data_root) / "indices" / sp["index_dir"]
        genome_dir = Path(sp["genome_dir"])
        cfg["references"]["genome_index"] = _find_star_index(idx_base)
        cfg["references"]["gtf"] = _find_gtf(genome_dir, "*")
        cfg["references"]["genome_fasta"] = _find_fasta(genome_dir, "*")

    if args.outdir:
        cfg["project"]["results_dir"] = args.outdir


def _validate_data(
    cfg: Dict[str, Any],
    strict: bool = False,
    steps: Optional[List[int]] = None,
) -> None:
    """Confirm data directories, metadata, FASTQs, and references exist.

    References (STAR index, GTF) are only required when running steps 5+.
    QC-only runs (steps 0-4) can proceed without them.
    """
    print("=" * 60)
    print("Data validation")
    print("=" * 60)

    dataset_dir = cfg["data"]["fastq_dir"]
    metadata_path = cfg["data"]["metadata_csv"]
    data_root = cfg.get("_data_root") or str(Path(dataset_dir).parent)

    # List available data under data root (helps when path is wrong)
    data_root_path = Path(data_root)
    print(f"Data root:    {data_root}")
    if data_root_path.exists():
        subdirs = sorted(d.name for d in data_root_path.iterdir() if d.is_dir())
        print(f"  Available datasets: {', '.join(subdirs) if subdirs else '(none)'}")
        indices_path = data_root_path / "indices"
        if indices_path.exists():
            genomes = sorted(d.name for d in indices_path.iterdir() if d.is_dir())
            print(f"  Reference genomes:  {', '.join(genomes) if genomes else '(none)'}")
    else:
        print(f"  WARNING: Data root not found: {data_root}")

    print(f"Dataset dir:  {dataset_dir}")
    print(f"Metadata:     {metadata_path}")
    print(f"STAR index:   {cfg['references']['genome_index']}")
    print(f"GTF:          {cfg['references']['gtf']}")

    errors: List[str] = []

    if not Path(dataset_dir).is_dir():
        errors.append(f"Dataset directory not found: {dataset_dir}")
    if not Path(metadata_path).is_file():
        errors.append(f"Metadata CSV not found: {metadata_path}")

    # References only needed for mapping/counting/DE (steps 5-9)
    steps_needing_refs = {5, 6, 7, 8, 9}
    steps_to_run = set(steps) if steps else steps_needing_refs
    if steps_to_run & steps_needing_refs:
        ref_idx = cfg["references"]["genome_index"]
        ref_gtf = cfg["references"]["gtf"]
        if not Path(ref_idx).exists():
            errors.append(f"STAR genome index not found: {ref_idx}")
            parent = Path(ref_idx).parent
            if parent.exists():
                contents = sorted(p.name for p in parent.iterdir())[:15]
                print(f"  Contents of {parent}: {', '.join(contents)}")
        if not Path(ref_gtf).exists():
            errors.append(f"GTF annotation not found: {ref_gtf}")
    else:
        print("  (Skipping reference check -- steps 0-4 only, no mapping/counting)")

    if errors:
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(f"FATAL: {len(errors)} path(s) not found. Check --dataset / --species.")

    run_col = cfg["column_mapping"]["run_id_col"]
    with open(metadata_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    run_ids = {row[run_col].strip() for row in rows if row.get(run_col)}

    fastqs = sorted(Path(dataset_dir).glob("*.fastq.gz"))
    fq_run_ids: set = set()
    for fq in fastqs:
        m = re.match(r"([SED]RR\d+|[SED]RA\d+)", fq.name)
        if m:
            fq_run_ids.add(m.group(1))

    matched = run_ids & fq_run_ids
    missing_fq = run_ids - fq_run_ids

    print(f"Runs in metadata:  {len(run_ids)}")
    print(f"FASTQ files found: {len(fastqs)} ({len(fq_run_ids)} unique run IDs)")
    print(f"Matched runs:      {len(matched)}")

    if len(matched) == 0:
        sys.exit("FATAL: No metadata run IDs match any FASTQ filenames in dataset dir")

    if missing_fq:
        msg = f"{len(missing_fq)} run(s) in metadata have no FASTQ: {sorted(missing_fq)[:5]}"
        if strict:
            sys.exit(f"FATAL (--strict): {msg}")
        print(f"  WARNING: {msg}")

    print("=" * 60)


# ---------------------------------------------------------------------------
# Step registry
# ---------------------------------------------------------------------------

STEPS = {
    0:  ("00_validate_env",     "Validate environment"),
    1:  ("01_prepare_samples",  "Prepare samples"),
    2:  ("02_trim_reads",       "Trim reads"),
    3:  ("03_qc_fastqc",       "FastQC"),
    4:  ("04_multiqc",         "MultiQC"),
    5:  ("05_map_star",        "STAR mapping"),
    6:  ("06_bigwig",          "BigWig generation"),
    7:  ("07_featurecounts",   "featureCounts"),
    8:  ("08_filter_matrix",   "Filter count matrix"),
    9:  ("09_deseq2",         "DESeq2"),
    10: ("10_compare_methods", "Compare methods"),
    11: ("11_make_report",     "Generate report"),
}

ALL_STEPS = sorted(STEPS.keys())


def import_step(step_num: int):
    """Dynamically import a step module."""
    module_name, _ = STEPS[step_num]
    # Import from scripts package
    import importlib
    mod = importlib.import_module(f"scripts.{module_name}")
    return mod


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def run_pipeline(
    cfg: Dict[str, Any],
    run_id: str,
    steps: Optional[List[int]] = None,
    methods_override: Optional[List[str]] = None,
    subset_name: Optional[str] = None,
) -> None:
    """Run the pipeline steps in order.

    Parameters
    ----------
    cfg : dict
        Loaded config.
    run_id : str
        Unique run identifier.
    steps : list of int, optional
        Which steps to run.  Default: all.
    methods_override : list of str, optional
        Override trimming methods from config.
    subset_name : str, optional
        Override subset filter from config.
    """
    steps = steps if steps is not None else ALL_STEPS

    logger.info("=" * 70)
    logger.info("RNA-seq Pipeline   |   Run ID: %s", run_id)
    logger.info("Steps to run: %s", steps)
    logger.info("=" * 70)

    t0 = time.time()

    for step_num in steps:
        if step_num not in STEPS:
            logger.error("Unknown step: %d (valid: %s)", step_num, ALL_STEPS)
            continue

        _, description = STEPS[step_num]
        logger.info("\n>>> STEP %d: %s", step_num, description)
        step_start = time.time()

        mod = import_step(step_num)

        try:
            if step_num == 0:
                cfg = mod.main(cfg, run_id)
            elif step_num == 1:
                mod.main(cfg, subset_name=subset_name)
            elif step_num in (2, 3, 4, 5, 6, 7, 8):
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 9:
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 10:
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 11:
                mod.main(cfg)
            else:
                mod.main(cfg)

        except Exception as exc:
            logger.error(
                "STEP %d FAILED: %s: %s", step_num, type(exc).__name__, exc,
                exc_info=True,
            )
            raise

        elapsed = time.time() - step_start
        logger.info("<<< STEP %d completed in %.1f s", step_num, elapsed)

    total = time.time() - t0
    logger.info("\n" + "=" * 70)
    logger.info("Pipeline complete.  Total time: %.1f s", total)
    logger.info("Results directory: %s", cfg.get("_results_dir", "N/A"))
    logger.info("=" * 70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description="RNA-seq Analysis Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="Path to config.yaml (default: config/config.yaml)",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Override run ID (default: auto-generated from project name + timestamp)",
    )

    # --- Data pointer args (override config paths) --------------------------
    parser.add_argument(
        "--data-root",
        default="/data",
        help="Root data directory on server (default: /data)",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="Dataset accession (e.g. GSE48519). Derives fastq_dir and metadata_csv.",
    )
    parser.add_argument(
        "--species",
        default=None,
        choices=sorted(SPECIES_MAP.keys()),
        help="Species name -- derives genome index, GTF, FASTA paths.",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="Explicit metadata CSV path (overrides --dataset default).",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Override results directory (default: results/<run_id>).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if any metadata run IDs lack matching FASTQ files.",
    )

    sub = parser.add_subparsers(dest="command", help="Pipeline commands")

    # --- run subcommand -----------------------------------------------------
    run_p = sub.add_parser("run", help="Run the pipeline (all steps or selected)")
    run_p.add_argument(
        "--steps",
        nargs="*",
        type=int,
        default=None,
        help=f"Step numbers to run (default: all). Valid: {ALL_STEPS}",
    )
    run_p.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Override trimming methods (e.g., none cutadapt fastp)",
    )
    run_p.add_argument(
        "--subset",
        default=None,
        help="Subset filter name from config (e.g., day3)",
    )
    run_p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Override project thread count (e.g. 1 to reduce memory on Colab)",
    )

    # --- list subcommand ----------------------------------------------------
    sub.add_parser("list", help="List available pipeline steps")

    return parser


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "list":
        print("\nAvailable pipeline steps:")
        print("-" * 50)
        for num in ALL_STEPS:
            module, desc = STEPS[num]
            print(f"  {num:2d}  {desc:30s}  ({module}.py)")
        print()
        sys.exit(0)

    # --- Runtime check (interpreter + packages) -----------------------------
    _check_runtime()

    # --- Load config --------------------------------------------------------
    cfg = load_config(args.config)

    # --- Apply data pointer overrides from CLI ------------------------------
    _apply_data_pointer(args, cfg)

    run_id = args.run_id or get_run_id(cfg)

    # Set up logging
    setup_logging(cfg["project"]["logs_dir"], run_id)

    # --- Validate data paths before heavy compute ---------------------------
    steps_to_run = getattr(args, "steps", None) or ALL_STEPS
    _validate_data(cfg, strict=args.strict, steps=steps_to_run)

    # Pre-enrich config with runtime paths (step 00 will enrich further)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    if args.command == "run" and getattr(args, "threads", None) is not None:
        cfg["project"]["threads"] = args.threads

    # --- Run ----------------------------------------------------------------
    if args.command == "run":
        run_pipeline(
            cfg,
            run_id,
            steps=args.steps,
            methods_override=args.methods,
            subset_name=args.subset,
        )


if __name__ == "__main__":
    main()
