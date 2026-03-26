#!/usr/bin/env python3
"""
Metadata CSV Inspector.

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
from typing import Any, Dict, List, Optional, Set, Tuple

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


# ---------------------------------------------------------------------------
# Smart biological short-name generation (used by --write-config only)
# ---------------------------------------------------------------------------

_BIO_TERMS: List[Tuple[re.Pattern, str]] = [
    (re.compile(r"^wild[\s-]?type$", re.I), "wt"),
    (re.compile(r"^(wt|wildtype)$", re.I), "wt"),
    (re.compile(r"^control$", re.I), "ctrl"),
    (re.compile(r"negative\s+control", re.I), "ctrl"),
    (re.compile(r"^untreated$", re.I), "untreated"),
    (re.compile(r"^treated$", re.I), "treated"),
    (re.compile(r"^knockout$", re.I), "ko"),
    (re.compile(r"^knock[\s-]?out$", re.I), "ko"),
    (re.compile(r"over[\s-]?express\w*", re.I), "oe"),
    (re.compile(r"^heterozygous$", re.I), "het"),
    (re.compile(r"^homozygous$", re.I), "hom"),
]

_CONTROL_KEYWORDS = {"wt", "ctrl", "control", "untreated", "wild-type", "wildtype"}

_SPECIES_REFS: Dict[str, Dict[str, str]] = {
    "mouse": {
        "genome_index": "/data/indices/mm39/star",
        "gtf": "/data/genomes/mouse/GRCm39/Mus_musculus.GRCm39.115.chr_gene_names.gtf.gz",
        "genome_fasta": "/data/genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
        "index_dir": "mm39",
    },
    "human": {
        "genome_index": "/data/indices/hg38/star",
        "gtf": "/data/genomes/human/hg38/Homo_sapiens.GRCh38.115.chr.gtf.gz",
        "genome_fasta": "/data/genomes/human/hg38/hg38.fa",
        "index_dir": "hg38",
    },
}


def _smart_short_name(text: str) -> str:
    """Generate a concise, biologically meaningful short name for DE."""
    text = text.strip()
    if not text:
        return "condition"

    for pattern, replacement in _BIO_TERMS:
        if pattern.search(text):
            return replacement

    # Knockout notation: geneX-/- or geneX -/-
    m = re.match(r"^(\w+)\s*-/-", text)
    if m:
        return m.group(1).lower()

    # Short alphanumeric already (e.g. "tet1")
    if len(text) <= 15 and re.match(r"^[A-Za-z0-9_\-]+$", text):
        return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()

    # Fallback: lowercase first 3 meaningful words, cap at 15
    words = re.findall(r"[A-Za-z0-9]+", text)
    short = "_".join(w.lower() for w in words[:3])
    if len(short) > 15:
        short = short[:15].rstrip("_")
    return short or "condition"


def _strip_shared_affixes(values: List[str]) -> List[str]:
    """Remove leading/trailing text shared by ALL condition values."""
    if len(values) < 2:
        return values

    originals = list(values)

    # Shared suffix (split on comma/semicolon boundaries)
    def _common_suffix(strings: List[str]) -> str:
        rev = [s[::-1] for s in strings]
        prefix = []
        for chars in zip(*rev):
            if len(set(chars)) == 1:
                prefix.append(chars[0])
            else:
                break
        return "".join(prefix)[::-1]

    suffix = _common_suffix(values)
    # Only strip at word/punctuation boundaries
    m = re.match(r"^[\s,;]+", suffix)
    if m:
        strip_len = len(suffix)
        values = [v[: len(v) - strip_len].strip() for v in values]
    elif len(suffix) > 3:
        strip_len = len(suffix)
        values = [v[: len(v) - strip_len].strip().rstrip(",;: ") for v in values]

    # Shared prefix
    def _common_prefix(strings: List[str]) -> str:
        prefix = []
        for chars in zip(*strings):
            if len(set(chars)) == 1:
                prefix.append(chars[0])
            else:
                break
        return "".join(prefix)

    prefix = _common_prefix(values)
    if len(prefix) > 3:
        m = re.search(r"[\s,;]+$", prefix)
        if m:
            strip_len = len(prefix)
            values = [v[strip_len:].strip() for v in values]

    return [v if v else orig for v, orig in zip(values, originals)]


def _dedup_labels(label_map: Dict[str, str]) -> Dict[str, str]:
    """Resolve collisions where multiple raw values map to the same short name."""
    seen: Dict[str, List[str]] = defaultdict(list)
    for raw, short in label_map.items():
        seen[short].append(raw)

    result: Dict[str, str] = {}
    for raw, short in label_map.items():
        if len(seen[short]) > 1:
            idx = seen[short].index(raw)
            if idx > 0:
                result[raw] = f"{short}{idx + 1}"
            else:
                result[raw] = short
        else:
            result[raw] = short
    return result


def _safe_key_name(text: str) -> str:
    """Generate a clean, concise YAML-safe key from a value string.

    Prefers the trailing (distinguishing) portion of compound labels.
    Examples: "differentiation day 3" -> "day_3", "day 0" -> "day_0".
    """
    text = text.strip().lower()
    key = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    key = re.sub(r"_+", "_", key)
    if len(key) <= 15:
        return key
    parts = key.split("_")
    for i in range(len(parts)):
        candidate = "_".join(parts[i:])
        if len(candidate) <= 15:
            return candidate
    return key[:15].rstrip("_")


# ---------------------------------------------------------------------------
# FASTQ matching and layout inference
# ---------------------------------------------------------------------------

def _extract_fq_run_ids(dataset_dir: Path) -> Tuple[Set[str], List[Path]]:
    """Glob FASTQs and extract SRA run IDs."""
    fastqs = sorted(dataset_dir.glob("*.fastq.gz"))
    fq_runs: Set[str] = set()
    for fq in fastqs:
        m = re.match(r"([SED]RR\d+)", fq.name)
        if m:
            fq_runs.add(m.group(1))
    return fq_runs, fastqs


def _match_to_fastqs(
    rows: List[Dict[str, str]],
    dataset_dir: Path,
    assay: str,
) -> Tuple[List[Dict[str, str]], Set[str], Set[str]]:
    """Filter metadata rows to those with FASTQs on disk AND matching assay.

    Returns (matched_rows, matched_run_ids, all_fq_run_ids).
    """
    fq_runs, _ = _extract_fq_run_ids(dataset_dir)
    matched_rows: List[Dict[str, str]] = []
    matched_ids: Set[str] = set()
    for row in rows:
        run_id = row.get("Run", "").strip()
        row_assay = row.get("Assay Type", "").strip()
        if run_id and run_id in fq_runs and row_assay == assay:
            matched_rows.append(row)
            matched_ids.add(run_id)
    return matched_rows, matched_ids, fq_runs


def _infer_layout(
    matched_rows: List[Dict[str, str]],
    dataset_dir: Path,
) -> Tuple[str, bool]:
    """Infer paired/single from matched FASTQ files on disk.

    Returns (layout_string, is_confident).
    Cross-references LibraryLayout column when available.
    """
    if not matched_rows:
        return "paired", False

    paired_count = 0
    single_count = 0
    for row in matched_rows:
        run_id = row.get("Run", "").strip()
        if not run_id:
            continue
        has_r1 = any(dataset_dir.glob(f"{run_id}*_1.fastq.gz"))
        has_r2 = any(dataset_dir.glob(f"{run_id}*_2.fastq.gz"))
        if has_r1 and has_r2:
            paired_count += 1
        elif has_r1 or any(dataset_dir.glob(f"{run_id}*.fastq.gz")):
            single_count += 1

    total = paired_count + single_count
    if total == 0:
        lib_layouts = Counter(
            r.get("LibraryLayout", "").strip().upper() for r in matched_rows
        )
        if lib_layouts.get("PAIRED", 0) > lib_layouts.get("SINGLE", 0):
            return "paired", False
        elif lib_layouts.get("SINGLE", 0) > 0:
            return "single", False
        return "paired", False

    if paired_count / total >= 0.8:
        return "paired", True
    if single_count / total >= 0.8:
        return "single", True
    return "paired", False


# ---------------------------------------------------------------------------
# Condition column selection (matched assay subset only)
# ---------------------------------------------------------------------------

def _pick_condition_column(
    matched_rows: List[Dict[str, str]],
    columns: List[str],
) -> Optional[Dict[str, Any]]:
    """Pick the best condition column from the matched (assay-filtered) subset.

    Scoring favours columns with 2-6 balanced, non-ID-like unique values.
    Only rows that passed the assay + FASTQ filter are considered.
    """
    candidates: List[Dict[str, Any]] = []

    for col in columns:
        if col in SKIP_COLS:
            continue
        values = Counter(
            r.get(col, "").strip()
            for r in matched_rows
            if r.get(col, "").strip()
        )
        n_unique = len(values)
        if n_unique < MIN_UNIQUE_FOR_CONDITION or n_unique > MAX_UNIQUE_FOR_CONDITION:
            continue

        counts = list(values.values())
        max_c, min_c = max(counts), min(counts)
        balance = max_c / min_c if min_c > 0 else 999.0

        # Penalize columns whose values look like IDs
        id_like = sum(
            1 for v in values if re.match(r"^[\d]+$", v) or len(v) > 60
        )

        score = 0.0
        if 2 <= n_unique <= 6:
            score += 10
        elif n_unique <= 10:
            score += 5
        if balance < 3:
            score += 5
        elif balance < 5:
            score += 2
        score -= id_like * 3
        # Slight preference for fewer groups
        score -= (n_unique - 2) * 0.5

        candidates.append({
            "col": col,
            "values": values,
            "n_unique": n_unique,
            "score": score,
        })

    if not candidates:
        return None
    return max(candidates, key=lambda c: c["score"])


def _find_secondary_grouping_columns(
    matched_rows: List[Dict[str, str]],
    columns: List[str],
    chosen_condition_col: str,
) -> List[Dict[str, Any]]:
    """Return plausible secondary biological grouping columns.

    Searches the matched (assay-filtered) subset for columns that also vary,
    excluding SKIP_COLS and the chosen condition column.  Accepts candidates
    in two modes:

    A. 2-6 balanced non-empty unique values (standard case).
    B. Exactly 1 non-empty unique value with at least 1 missing/blank row
       (partially-annotated column -- still signals a biological subgroup).

    No column names are hardcoded.
    """
    total = len(matched_rows)
    results: List[Dict[str, Any]] = []
    for col in columns:
        if col in SKIP_COLS or col == chosen_condition_col:
            continue

        all_vals = [r.get(col, "").strip() for r in matched_rows]
        non_empty = Counter(v for v in all_vals if v)
        n_missing = sum(1 for v in all_vals if not v)
        n_unique = len(non_empty)

        if n_unique == 0:
            continue

        id_like = sum(
            1 for v in non_empty if re.match(r"^\d+$", v) or len(v) > 60
        )
        if id_like > 0:
            continue

        accept = False
        score = 0.0

        if 2 <= n_unique <= 6:
            counts = list(non_empty.values())
            max_c, min_c = max(counts), min(counts)
            balance = max_c / min_c if min_c > 0 else 999.0
            if balance < 5:
                accept = True
                if 2 <= n_unique <= 4:
                    score += 10
                else:
                    score += 5
                if balance < 3:
                    score += 5
                score -= (n_unique - 2) * 0.5
        elif n_unique == 1 and n_missing > 0:
            accept = True
            score += 6

        if not accept:
            continue

        results.append({
            "col": col,
            "values": non_empty,
            "n_unique": n_unique,
            "n_missing": n_missing,
            "score": score,
        })

    results.sort(key=lambda c: c["score"], reverse=True)
    return results


def _detect_species(organism: str) -> str:
    """Map organism string to species flag."""
    org = organism.lower()
    if "mus musculus" in org or "mouse" in org:
        return "mouse"
    if "homo sapiens" in org or "human" in org:
        return "human"
    return "UNKNOWN"


# ---------------------------------------------------------------------------
# Full draft YAML config generation
# ---------------------------------------------------------------------------

def generate_config(
    csv_path: str,
    output_path: str,
    *,
    assay: str = "RNA-Seq",
    dataset_name: Optional[str] = None,
    condition_col_override: Optional[str] = None,
    reference_level_override: Optional[str] = None,
    layout_override: Optional[str] = None,
    minimal: bool = False,
) -> None:
    """Generate a full draft YAML config for the RNA-seq pipeline.

    Bases all recommendations on the matched subset: rows whose Run ID has
    FASTQs on disk AND whose Assay Type matches ``assay``.
    """
    path = Path(csv_path)
    if not path.exists():
        sys.exit(f"File not found: {csv_path}")

    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            sys.exit("Could not parse CSV headers")
        columns = list(reader.fieldnames)
        all_rows = list(reader)

    dataset_dir = path.parent
    if dataset_name is None:
        dataset_name = dataset_dir.name

    # --- Match to FASTQs + assay filter ---
    matched_rows, matched_ids, fq_runs = _match_to_fastqs(
        all_rows, dataset_dir, assay,
    )
    if not matched_rows:
        sys.exit(
            f"No rows matched assay={assay!r} with FASTQs in {dataset_dir}.\n"
            f"  Total metadata rows: {len(all_rows)}\n"
            f"  FASTQ run IDs found: {len(fq_runs)}\n"
            f"  Try a different --assay value."
        )

    # --- Layout ---
    todos: List[str] = []
    if layout_override:
        layout = layout_override
        layout_confident = True
    else:
        layout, layout_confident = _infer_layout(matched_rows, dataset_dir)
    layout_todo = ""
    if not layout_confident:
        layout_todo = "  # TODO: verify layout"
        todos.append("Verify layout")

    # --- Species ---
    organisms = Counter(r.get("Organism", "") for r in matched_rows)
    primary_organism = organisms.most_common(1)[0][0] if organisms else ""
    species = _detect_species(primary_organism)

    if species in _SPECIES_REFS:
        refs = _SPECIES_REFS[species]
        ref_index = refs["genome_index"]
        ref_gtf = refs["gtf"]
        ref_fasta = refs["genome_fasta"]
        ref_comment = f"# Auto-detected: {species}"
        idx_dir = refs["index_dir"]
    else:
        ref_index = "/data/indices/UNKNOWN/star"
        ref_gtf = "/data/genomes/UNKNOWN/UNKNOWN.gtf.gz"
        ref_fasta = "/data/genomes/UNKNOWN/UNKNOWN.fa.gz"
        ref_comment = "# TODO: set reference genome paths (unknown species)"
        idx_dir = "UNKNOWN"
        todos.append("Set reference genome paths")

    # --- Condition column (from matched assay subset only) ---
    if condition_col_override:
        cond_values = Counter(
            r.get(condition_col_override, "").strip()
            for r in matched_rows
            if r.get(condition_col_override, "").strip()
        )
        best_col = condition_col_override
        cond_todo = ""
        if len(cond_values) < 2:
            cond_todo = "  # TODO: column has <2 values in matched subset"
            todos.append("Verify condition column (few values)")
    else:
        picked = _pick_condition_column(matched_rows, columns)
        if picked:
            best_col = picked["col"]
            cond_values = picked["values"]
            cond_todo = "  # TODO: verify condition column choice"
            todos.append("Verify condition column choice")
        else:
            best_col = "CHANGE_ME"
            cond_values = Counter()
            cond_todo = "  # TODO: no condition column auto-detected"
            todos.append("Set condition column manually")

    # --- Short names ---
    raw_values = sorted(cond_values.keys())
    stripped = _strip_shared_affixes(raw_values)
    label_map: Dict[str, str] = {}
    for raw, cleaned in zip(raw_values, stripped):
        label_map[raw] = _smart_short_name(cleaned)
    label_map = _dedup_labels(label_map)
    short_names = [label_map[v] for v in raw_values]

    # --- Reference level ---
    if reference_level_override:
        ref_level = reference_level_override
        ref_level_todo = ""
    else:
        ref_level = short_names[0] if short_names else "ctrl"
        for raw, short in label_map.items():
            if short.lower() in _CONTROL_KEYWORDS or raw.lower() in _CONTROL_KEYWORDS:
                ref_level = short
                break
        ref_level_todo = "  # TODO: confirm reference level"
        todos.append("Confirm reference level")

    # --- Comparisons (cautious for 3+ groups) ---
    treatments = [sn for sn in short_names if sn != ref_level]
    comparisons_lines: List[str] = []
    if not treatments:
        comparisons_lines.append(
            '  - name: "treatment_vs_ctrl"  # TODO: set comparison\n'
            '    numerator: "treatment"\n'
            '    denominator: "ctrl"'
        )
        todos.append("Set comparison manually")
    elif len(treatments) == 1:
        t = treatments[0]
        comparisons_lines.append(
            f'  - name: "{t}_vs_{ref_level}"\n'
            f'    numerator: "{t}"\n'
            f'    denominator: "{ref_level}"'
        )
    else:
        comparisons_lines.append(
            "  # NOTE: multiple treatment groups detected.\n"
            "  # Only pairwise-vs-reference comparisons are generated.\n"
            "  # Add or remove contrasts as needed."
        )
        for t in treatments:
            comparisons_lines.append(
                f'  - name: "{t}_vs_{ref_level}"\n'
                f'    numerator: "{t}"\n'
                f'    denominator: "{ref_level}"'
            )
        todos.append("Review comparisons (multiple treatment groups)")

    # --- Condition map YAML fragment ---
    cmap_lines = []
    for raw in raw_values:
        cmap_lines.append(f'    "{raw}": "{label_map[raw]}"')
    cmap_block = "\n".join(cmap_lines)

    comparisons_block = "\n".join(comparisons_lines)

    # --- Condition counts for summary ---
    cond_summary_parts = []
    for raw in raw_values:
        cond_summary_parts.append(f"{label_map[raw]} ({cond_values[raw]})")

    # --- Secondary grouping columns (generic detection) ---
    secondary_cols = _find_secondary_grouping_columns(
        matched_rows, columns, best_col,
    )
    secondary_filter_block = ""
    secondary_summary_lines: List[str] = []
    if secondary_cols:
        top_sec = secondary_cols[0]
        sec_col = top_sec["col"]
        sec_vals = sorted(top_sec["values"].keys())
        n_missing = top_sec.get("n_missing", 0)

        sec_val_summary = ", ".join(
            f"{v} ({top_sec['values'][v]})" for v in sec_vals
        )
        secondary_summary_lines.append(
            f"  Secondary grouping: {sec_col}"
        )
        secondary_summary_lines.append(
            f"    Non-empty values: {sec_val_summary}"
        )
        if n_missing > 0:
            secondary_summary_lines.append(
                f"    Missing values: {n_missing}"
            )
            secondary_summary_lines.append(
                "    WARNING: matched subset mixes annotated and "
                "unannotated samples for this column."
            )
        else:
            secondary_summary_lines.append(
                "    WARNING: another biological axis varies in the "
                "matched subset."
            )
        secondary_summary_lines.append(
            "    Subgroup subset_filters have been added to the draft config."
        )
        todos.append(
            f"Review secondary grouping column: {sec_col} "
            f"({top_sec['n_unique']} non-empty values"
            f"{f', {n_missing} missing' if n_missing else ''})"
        )

        sec_lines = []
        sec_lines.append("")
        sec_lines.append(
            f'  # WARNING: column "{sec_col}" also varies in the '
            f"matched subset."
        )
        for sv in sec_vals:
            sec_lines.append(
                f"  # Non-empty values: {sv} ({top_sec['values'][sv]})"
            )
        if n_missing > 0:
            sec_lines.append(f"  # Missing/blank rows: {n_missing}")
        sec_lines.append(
            "  # Consider switching active_subset to one subgroup "
            "for a cleaner DE comparison."
        )
        for sv in sec_vals:
            filter_name = _safe_key_name(sv)
            sec_lines.append(f"  {filter_name}:")
            sec_lines.append(f'    "Assay Type": "{assay}"')
            sec_lines.append(f'    "{sec_col}": "{sv}"')
        secondary_filter_block = "\n".join(sec_lines)

    # --- Build conditional YAML sections ---

    # Mapping section
    if minimal:
        mapping_block = f"""\
mapping:
  backends:
    star:
      enabled: true
      genome_index: "{ref_index}"
      option_sets:
        default:
          outFilterMultimapNmax: 20
          extra_args: ""
    hisat2:
      enabled: false"""
    else:
        mapping_block = f"""\
mapping:
  backends:
    star:
      enabled: true
      genome_index: "{ref_index}"
      option_sets:
        default:
          outFilterMultimapNmax: 20
          extra_args: ""
        strict_unique:
          outFilterMultimapNmax: 1
          extra_args: ""
    hisat2:
      enabled: false
      index_prefix: ""
      option_sets:
        default:
          extra_args: \"\""""

    # Trimming section
    if minimal:
        trimming_block = """\
trimming:
  primary_method: "cutadapt"
  compare_methods: false
  comparison_methods: []

  none:
    enabled: false

  cutadapt:
    enabled: true
    adapter_fwd: ""
    adapter_rev: ""
    quality: 20
    min_length: 25
    extra_args: \"\""""
    else:
        trimming_block = """\
trimming:
  primary_method: "cutadapt"
  compare_methods: false
  comparison_methods: []

  none:
    enabled: true

  cutadapt:
    enabled: true
    adapter_fwd: ""
    adapter_rev: ""
    quality: 20
    min_length: 25
    extra_args: ""

  fastp:
    enabled: false
    quality: 20
    min_length: 25
    detect_adapter: true
    extra_args: ""

  trimmomatic:
    enabled: false
    adapter_file: "TruSeq3-PE-2.fa"
    leading: 3
    trailing: 3
    slidingwindow: "4:20"
    min_length: 25
    extra_args: \"\""""

    # featureCounts section
    if minimal:
        fc_block = """\
featurecounts:
  backend: "featurecounts"
  strandedness: 0  # TODO: check dataset protocol (0=unstranded, 1=stranded, 2=reverse-stranded)
  use_filtered_bam: false
  feature_type: "exon"
  attribute: "gene_name"  # TODO: change to "gene_id" if Ensembl IDs are preferred
  gtf_filter: "gene_name"
  option_sets:
    default:
      label: "Default"
      B: false
      P: false
      C: false
      Q: 0
      extra_args: ""
    # multimapper:
    #   label: "Multi-mapper aware (M+fraction)"
    #   M: true
    #   fraction: true
    #   B: false
    #   P: false
    #   C: false
    #   Q: 0
    #   extra_args: \"\""""
    else:
        fc_block = """\
featurecounts:
  backend: "featurecounts"
  strandedness: 0  # TODO: check dataset protocol (0=unstranded, 1=stranded, 2=reverse-stranded)
  use_filtered_bam: false
  feature_type: "exon"
  attribute: "gene_name"  # TODO: change to "gene_id" if Ensembl IDs are preferred
  gtf_filter: "gene_name"
  option_sets:
    default:
      label: "Default"
      B: false
      P: false
      C: false
      Q: 0
      extra_args: ""
    strict:
      label: "Strict (B+P+C+Q10)"
      B: true
      P: true
      C: true
      Q: 10
      extra_args: ""
    stringent:
      label: "Stringent (B+P+C+Q255)"
      B: true
      P: true
      C: true
      Q: 255
      extra_args: ""
    multimapper:
      label: "Multi-mapper aware (M+fraction)"
      M: true
      fraction: true
      B: false
      P: false
      C: false
      Q: 0
      extra_args: \"\""""

    # Selection policy block
    if minimal:
        selection_block = """\
# -- Selection policy ---------------------------------------------------------
# Single-branch config: only cutadapt / star / default / default is generated.
# Uncomment and adjust if you add more branches later.
#
# selection:
#   preferred_trim_method: "cutadapt"
#   preferred_mapper: "star"
#   preferred_mapper_option_set: "default"
#   primary_count_option_set: "default"
#   comparison_count_option_sets:
#     - "default\""""
    else:
        selection_block = """\
# -- Selection policy (optional) ----------------------------------------------
# Uncomment and fill in to override automatic branch selection in step 10.
# selected_analysis.tsv and selected_visualisation.tsv -> primary branch.
# selected_count_comparison.tsv -> comparison branch.
#
# selection:
#   preferred_trim_method: "cutadapt"
#   preferred_mapper: "star"
#   preferred_mapper_option_set: "default"
#   primary_count_option_set: "default"
#   comparison_count_option_sets:
#     - "default"
#     - "multimapper\""""

    # --- Build YAML ---
    yaml_text = f"""\
# =============================================================================
# RNA-seq Pipeline Configuration -- {dataset_name}
# =============================================================================
# Auto-generated by inspect_metadata.py --write-config
# Source: {path.name}
# Assay filter: {assay}
# Matched runs: {len(matched_rows)} / {len(all_rows)}
{"# Mode: minimal (single-branch)" if minimal else "# Mode: full"}
# =============================================================================

# -- Project metadata ---------------------------------------------------------
project:
  name: "{dataset_name}"
  results_dir: "results"
  work_dir: "work"
  cache_dir: "cache"
  logs_dir: "logs"
  threads: 6

# -- Input data ---------------------------------------------------------------
data:
  metadata_csv: "/data/{dataset_name}/metadata.csv"
  fastq_dir: "/data/{dataset_name}"
  layout: "{layout}"{layout_todo}

# -- Column mapping -----------------------------------------------------------
{cond_todo}
column_mapping:
  run_id_col: "Run"
  condition_cols:
    - "{best_col}"
  condition_map:
{cmap_block}
  replicate_strategy: "sorted_run_id"

# -- Subset filters -----------------------------------------------------------
# Filter to {assay} rows only.
subset_filters:
  default:
    "Assay Type": "{assay}"
{secondary_filter_block}
active_subset: "default"{f'  # TODO: consider switching to a subgroup filter' if secondary_filter_block else ''}

# -- Reference genome ---------------------------------------------------------
{ref_comment}
references:
  genome_index: "{ref_index}"
  gtf: "{ref_gtf}"
  genome_fasta: "{ref_fasta}"

# -- Tool paths ---------------------------------------------------------------
tools:
  fastqc: "fastqc"
  multiqc: "multiqc"
  cutadapt: "cutadapt"
  fastp: "fastp"
  trimmomatic: "trimmomatic"
  star: "STAR"
  hisat2: "hisat2"
  samtools: "samtools"
  featurecounts: "featureCounts"
  bamcoverage: "bamCoverage"

# -- Mapping approaches -------------------------------------------------------
{mapping_block}

# -- STAR mapping parameters --------------------------------------------------
star_params:
  outSAMtype: "BAM SortedByCoordinate"
  readFilesCommand: "zcat"
  outFilterMultimapNmax: 20
  extra_args: ""

# -- Trimming methods ---------------------------------------------------------
{trimming_block}

# -- featureCounts options ----------------------------------------------------
{fc_block}

# -- Matrix filtering ---------------------------------------------------------
filtering:
  min_max_count: "auto"

# -- DESeq2 parameters --------------------------------------------------------
deseq2:
  fdr_threshold: 0.05
  lfc_threshold: 0.0
  method: "pydeseq2"
  reference_level: "{ref_level}"{ref_level_todo}

# -- Contrasts / comparisons --------------------------------------------------
comparisons:
{comparisons_block}

{selection_block}

# -- BigWig parameters --------------------------------------------------------
bigwig:
  enabled: false
  mode: "selected_only"
  normalization: "CPM"
  mapq_filter: 255
  bin_size: 50
  use_deseq2_sizefactors: true
"""

    if "strandedness" not in [t.lower() for t in todos]:
        todos.append("Verify strandedness")

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml_text, encoding="utf-8")

    # --- Summary ---
    comparison_names = [f"{t}_vs_{ref_level}" for t in treatments]

    print(f"\nDraft config written: {out}")
    print(f"  Mode:            {'minimal' if minimal else 'full'}")
    print(f"  Assay filter:    {assay}")
    print(f"  Matched runs:    {len(matched_rows)} / {len(all_rows)}")
    print(f"  Layout:          {layout}{'  (confident)' if layout_confident else '  (uncertain)'}")
    print(f"  Condition col:   {best_col} ({len(cond_values)} values)")
    print(f"  Conditions:      {', '.join(cond_summary_parts)}")
    print(f"  Reference level: {ref_level}")
    print(f"  Comparisons:     {', '.join(comparison_names) if comparison_names else '(none)'}")
    print(f"  Species:         {species}")
    for line in secondary_summary_lines:
        print(line)
    if todos:
        print(f"\n  TODOs in generated config:")
        for td in todos:
            print(f"    - {td}")
    print()


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
    print(f"  --config config/config_{dataset_name}.yaml \\")
    print(f"  --dataset {dataset_name} \\")
    print(f"  --species {species_flag} \\")
    print(f"  run --profile primary --methods none cutadapt")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect SRA metadata CSV for pipeline configuration",
    )
    parser.add_argument("csv_path", help="Path to metadata.csv")
    parser.add_argument(
        "--write-config", default=None, metavar="PATH",
        help="Write a full draft YAML config to this path",
    )
    parser.add_argument(
        "--assay", default="RNA-Seq",
        help="Assay type filter for config generation (default: RNA-Seq)",
    )
    parser.add_argument(
        "--dataset-name", default=None,
        help="Short project name (default: inferred from CSV parent directory)",
    )
    parser.add_argument(
        "--condition-col", default=None,
        help="Override auto-detected condition column name",
    )
    parser.add_argument(
        "--reference-level", default=None,
        help="Override auto-detected reference level for DE",
    )
    parser.add_argument(
        "--layout", default=None, choices=["paired", "single"],
        help="Override auto-detected library layout",
    )
    parser.add_argument(
        "--minimal", action="store_true",
        help="Generate a lean single-branch config (cutadapt-only, STAR-default-only)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="When --write-config is active, also print the full legacy inspector dump",
    )
    args = parser.parse_args()

    if args.write_config:
        if args.verbose:
            inspect(args.csv_path)
        generate_config(
            args.csv_path,
            args.write_config,
            assay=args.assay,
            dataset_name=args.dataset_name,
            condition_col_override=args.condition_col,
            reference_level_override=args.reference_level,
            layout_override=args.layout,
            minimal=args.minimal,
        )
    else:
        inspect(args.csv_path)


if __name__ == "__main__":
    main()
