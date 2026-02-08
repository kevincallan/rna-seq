# RNA-seq Analysis Pipeline

A complete, reproducible RNA-seq analysis pipeline designed for exam-day use. Runs as standalone Python scripts (one per stage), driven entirely by `config/config.yaml` + a metadata CSV. Supports comparing multiple trimming methodologies side-by-side.

---

## Quick Start

```bash
cd pipeline/

# 1. Edit config/config.yaml with your paths
# 2. Run the full pipeline:
python scripts/run_pipeline.py --config config/config.yaml run

# Or use Make:
make run
```

---

## Repository Structure

```
pipeline/
  README.md                  # This file
  Makefile                   # Convenience targets
  .gitignore
  config/
    config.yaml              # ALL configuration lives here
  env/
    environment.yml          # Conda environment
    requirements.txt         # pip requirements
  src/
    __init__.py
    utils.py                 # Logging, subprocess, hashing, config loader
    metadata.py              # CSV parsing, design table, symlinks
    reporting.py             # Markdown/HTML report builder
  scripts/
    run_pipeline.py          # Orchestrator (run full or individual steps)
    00_validate_env.py       # Check tools, paths, create directories
    01_prepare_samples.py    # Parse metadata, create symlinks, samples.tsv
    02_trim_reads.py         # Trim with none/cutadapt/fastp/trimmomatic
    03_qc_fastqc.py          # FastQC on raw + trimmed reads
    04_multiqc.py            # MultiQC reports per method
    05_map_star.py           # STAR mapping per method
    06_bigwig.py             # BigWig (CPM + DESeq2-scaled) per method
    07_featurecounts.py      # featureCounts with multiple option sets
    08_filter_matrix.py      # Filter low-expression genes
    09_deseq2.py             # DESeq2 per method/contrast
    10_compare_methods.py    # Cross-method comparison report
    11_make_report.py        # Final Markdown + HTML report
    deseq2_run.R             # Standalone R script for DESeq2
    download_sra.py          # Optional: download FASTQs from SRA
  tests/
    test_metadata.py         # Tests for metadata parsing
    test_naming.py           # Tests for naming conventions
```

---

## Exam-Day Checklist

### What to change in `config/config.yaml`:

1. **`data.metadata_csv`** -- path to the new metadata CSV
2. **`data.fastq_dir`** -- path to the directory with FASTQs
3. **`column_mapping.run_id_col`** -- column name with Run IDs (e.g., `Run`)
4. **`column_mapping.condition_cols`** -- columns to build condition strings
5. **`column_mapping.condition_map`** -- map raw values to short names
6. **`subset_filters`** -- define your subset (e.g., "day 3 WT vs Tet1")
7. **`active_subset`** -- select which subset to use
8. **`references`** -- genome index, GTF, FASTA paths
9. **`comparisons`** -- which contrasts to test (numerator vs denominator)
10. **`deseq2.reference_level`** -- the control condition

### What NOT to change:

- Tool paths (unless they differ on your machine)
- Trimming method parameters (unless you want to experiment)
- Script code (it's all config-driven)

### Steps to run:

```bash
# Navigate to pipeline directory
cd pipeline/

# Full pipeline:
python scripts/run_pipeline.py --config config/config.yaml run

# With a specific subset:
python scripts/run_pipeline.py --config config/config.yaml run --subset day3_wt_vs_tet1

# Only specific trimming methods:
python scripts/run_pipeline.py --config config/config.yaml run --methods none cutadapt fastp

# Only specific steps (e.g., just QC):
python scripts/run_pipeline.py --config config/config.yaml run --steps 0 1 2 3 4
```

---

## Pipeline Stages

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_validate_env.py` | Validate config, check tools, create directories, write manifest |
| 1 | `01_prepare_samples.py` | Parse metadata, filter subset, create symlinks + samples.tsv |
| 2 | `02_trim_reads.py` | Trim reads with each enabled method (none/cutadapt/fastp) |
| 3 | `03_qc_fastqc.py` | Run FastQC on raw and trimmed reads |
| 4 | `04_multiqc.py` | Generate MultiQC reports per method |
| 5 | `05_map_star.py` | Map reads with STAR per method, index BAMs |
| 6 | `06_bigwig.py` | Generate BigWigs (CPM + optional DESeq2-scaled) |
| 7 | `07_featurecounts.py` | Count reads with featureCounts (multiple option sets) |
| 8 | `08_filter_matrix.py` | Filter low-expression genes |
| 9 | `09_deseq2.py` | Run DESeq2 for each method and contrast |
| 10 | `10_compare_methods.py` | Compare trimming methods (mapping, counts, DEGs) |
| 11 | `11_make_report.py` | Generate final Markdown + HTML report |

---

## Trimming Method Comparison

The pipeline runs all enabled trimming methods in parallel and tracks differences:

- **none**: No trimming (baseline)
- **cutadapt**: Adapter + quality trimming with cutadapt
- **fastp**: All-in-one trimming with fastp (auto-detects adapters)
- **trimmomatic**: Traditional sliding-window trimmer (disabled by default)

Each downstream step (mapping, counting, DE) runs independently per method. Step 10 produces a comparison report showing:

- Mapping rate differences
- Assigned reads differences
- Normalised count correlations
- DEG overlap (Venn-style counts)
- Method-sensitive genes (padj flips)

---

## featureCounts Option Comparison

The pipeline supports multiple featureCounts option sets (configured in YAML):

- **default**: No special flags
- **strict**: `-B -P -C -Q 10` (both ends mapped, proper pairs, no chimeras, MAPQ>=10)
- **stringent**: `-B -P -C -Q 255` (unique mappers only)

This follows the lecture homework: "Try different options in featureCounts (B, P, C, Q) and compare results."

---

## BigWig Generation

For each sample and method, three types of BigWig files can be generated:

1. **CPM-normalised** with MAPQ=255 (unique mappers only)
2. **DESeq2 size-factor scaled** with MAPQ=255 (if size factors available)

This follows the lecture: `bamCoverage --scaleFactor 1/FACTOR --minMappingQuality 255`

---

## Output Structure

```
results/<run_id>/
  samples.tsv                          # Design table
  sample_description.txt               # DESeq2 sample file
  mapping_summary.tsv                  # STAR stats across methods
  featurecounts_summary.tsv            # Assigned reads across methods
  filtering_summary.tsv                # Filtering stats
  de_summary.tsv                       # DEG counts across methods
  trimming_summary.tsv                 # Trimming metrics
  run_manifest.json                    # Reproducibility manifest
  reports/
    report.md                          # Final report (Markdown)
    report.html                        # Final report (HTML)
    method_comparison.md               # Detailed method comparison
    multiqc_raw/                       # Raw reads QC
    multiqc_none/                      # Per-method MultiQC
    multiqc_cutadapt/
    multiqc_fastp/
  none/                                # Per-method results
    star/                              # BAMs + STAR logs
    featurecounts/                     # Count matrices
    deseq2/                            # DE results + plots
    bigwig/                            # BigWig tracks
    qc/fastqc/                         # FastQC results
  cutadapt/
    ...
  fastp/
    ...
```

---

## Worked Example: GSE48519 (5mC Oxidation Paper)

### Undifferentiated WT vs Tet1 vs Tet2

This is the default configuration in `config.yaml`. The metadata filters for:
- `Assay Type == "RNA-Seq"`
- `differentiation == ""` (undifferentiated mESCs)

```bash
python scripts/run_pipeline.py --config config/config.yaml run
```

This maps 6 samples (wt_1, wt_2, tet1_1, tet1_2, tet2_1, tet2_2) and runs two contrasts: tet1 vs wt and tet2 vs wt.

### Day 3 WT vs Tet1 (Exam Homework)

To analyse the differentiation day 3 subset:

1. In `config.yaml`, set `active_subset: "day3_wt_vs_tet1"` or use the CLI:

```bash
python scripts/run_pipeline.py --config config/config.yaml run --subset day3_wt_vs_tet1
```

2. Update the `comparisons` section if needed for the new conditions.

### Using a Completely Different Dataset

1. Place your FASTQs in a directory (e.g., `/data/NEW_STUDY/`)
2. Place your metadata CSV alongside them
3. Edit `config.yaml`:

```yaml
data:
  metadata_csv: "/data/NEW_STUDY/metadata.csv"
  fastq_dir: "/data/NEW_STUDY"
  layout: "paired"   # or "single" or "auto"

column_mapping:
  run_id_col: "Run"                    # whatever your CSV calls it
  condition_cols: ["treatment"]        # your condition column(s)
  condition_map:
    "control": "ctrl"
    "drug_treated": "drug"

references:
  genome_index: "/data/indices/hg38/STAR"   # for human, etc.
  gtf: "/data/indices/hg38/hg38.gtf"

comparisons:
  - name: "drug_vs_ctrl"
    numerator: "drug"
    denominator: "ctrl"
```

4. Run: `python scripts/run_pipeline.py --config config/config.yaml run`

No code changes needed.

---

## Running Tests

```bash
cd pipeline/
python -m pytest tests/ -v
```

---

## Dependencies

Install via conda (recommended):

```bash
conda env create -f env/environment.yml
conda activate rnaseq_pipeline
```

Or pip (Python packages only -- bioinformatics tools must be installed separately):

```bash
pip install -r env/requirements.txt
```

Required external tools: STAR, samtools, featureCounts (subread), FastQC, MultiQC, cutadapt, fastp, bamCoverage (deeptools), R + DESeq2.

---

## Design Principles

- **One file per step**: Each script is independently runnable
- **Config-driven**: All parameters in `config/config.yaml`
- **Deterministic**: Same inputs + config = same outputs (sorted Run IDs for replicate numbering)
- **No hard-coded IDs**: Works with any dataset via config
- **Safe subprocess calls**: `subprocess.run([...], check=True)`, no `shell=True`
- **Full provenance**: `run_manifest.json` captures config hash, tool versions, timestamps
