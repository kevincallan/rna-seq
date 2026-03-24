# RNA-seq Analysis Pipeline

A config-driven, end-to-end RNA-seq pipeline for differential expression analysis. Produces QC reports, count matrices, DE results with statistical tables, publication-quality figures, and gene lists ready for functional enrichment analysis.

---

## Quick Start

```bash
cd ~/rna-seq-pipeline

# Primary analysis -- one command to a complete report:
./py scripts/run_pipeline.py --dataset GSEXXXXX --species mouse run --profile primary

# Full pipeline including BigWig coverage tracks:
./py scripts/run_pipeline.py --dataset GSEXXXXX --species mouse run --profile full
```

`./py` is a wrapper that ensures the correct interpreter is used. The `primary` profile runs the fastest path to a complete, report-ready result:

1. Verify the interpreter and all required packages
2. Derive data paths from `--dataset` (FASTQs and metadata in `/data/GSEXXXXX/`)
3. Validate that files exist and FASTQ filenames match metadata run IDs
4. Run QC, trimming, mapping, counting (all option sets), filtering, DE analysis (all option sets), comparison, branch selection, and report generation

---

## Execution Profiles

The pipeline supports three execution profiles:

| Profile | Command | What it does |
|---------|---------|-------------|
| `primary` | `--profile primary` | Fastest path to a complete report. Runs QC through DE, comparison, and report. **No BigWig.** This is the default recommended workflow. |
| `primary_bw` | `--profile primary_bw` | Same as primary, plus BigWig generation for the selected analysis branch. |
| `full` | `--profile full` | All steps including BigWig. |

```bash
# Report-ready run (recommended):
./py scripts/run_pipeline.py --config config/config.yaml \
    --dataset GSEXXXXX --species mouse run --profile primary

# Add BigWig tracks for the selected branch:
./py scripts/run_pipeline.py --config config/config.yaml \
    --dataset GSEXXXXX --species mouse run --profile primary_bw

# Everything:
./py scripts/run_pipeline.py --config config/config.yaml \
    --dataset GSEXXXXX --species mouse run --profile full
```

You can also select individual steps with `--steps`:

```bash
./py scripts/run_pipeline.py --dataset GSEXXXXX run --steps 0 1 2 3 4 5
```

---

## What the Pipeline Does

```
Raw FASTQs
  --> Quality control (FastQC / MultiQC)
  --> Adapter trimming (cutadapt, with untrimmed baseline)
  --> Read mapping to reference genome (STAR and/or HISAT2)
  --> Read counting per gene (featureCounts, all configured option sets)
  --> Low-expression gene filtering (all option sets)
  --> Differential expression (PyDESeq2, for every count option set)
  --> Method and parameter comparison (mapping, counting, DE overlap)
  --> Primary branch selection (automatic or manual)
  --> [Optional] BigWig coverage tracks for IGV (selected branch only)
  --> Final summary report with pipeline health checks
```

All steps are automated and config-driven. The pipeline compares multiple trimming approaches and featureCounts parameter sets side-by-side, runs DE for each, and selects a primary analysis branch.

---

## Pipeline Steps

| Step | Module | Description | Key outputs |
|------|--------|-------------|-------------|
| 0 | `00_validate_env` | Validate environment, tools, paths | `run_manifest.json` |
| 1 | `01_prepare_samples` | Parse metadata, apply subset filters, create symlinks | `samples.tsv` |
| 2 | `02_trim_reads` | Run each enabled trimming method | Trimmed FASTQs, `trimming_summary.tsv` |
| 3 | `03_qc_fastqc` | FastQC on raw and trimmed reads | FastQC HTML reports |
| 4 | `04_multiqc` | Aggregate FastQC into MultiQC | MultiQC HTML reports |
| 5 | `05_map_star` | Align reads (STAR/HISAT2), validate BAM integrity, then index BAMs | Sorted/indexed BAMs, `mapping_summary.tsv` |
| 6 | `07_featurecounts` | Count reads per gene (all option sets) | Count matrices, `featurecounts_summary.tsv` |
| 7 | `08_filter_matrix` | Remove lowly-expressed genes | Filtered matrices, `filtering_summary.tsv` |
| 8 | `09_deseq2` | DE analysis for every count option set | DE tables, PCA/MA/volcano plots, gene lists |
| 9 | `10_compare_methods` | Compare methods, select count-comparison + visualisation branches | `method_comparison.md`, `selected_count_comparison.tsv`, `selected_visualisation.tsv`, `selected_analysis.tsv` (legacy) |
| 10 | `06_bigwig` | Generate coverage tracks (CPM + DESeq2-scaled) | BigWig files for IGV |
| 11 | `11_make_report` | Generate combined report with health checks | `report.md`, `report.html` |

BigWig generation (step 10) runs **after** DESeq2 so that size factors are available for scaled tracks. The `primary` profile skips it by default for speed.

---

## Key Design: count_option_set as a First-Class Dimension

The pipeline treats the featureCounts option set as a first-class analysis dimension alongside trimming method, mapper, and mapper option set. This means:

- **Step 6 (featureCounts)** runs all configured option sets (default, strict, stringent, multimapper)
- **Step 7 (filter)** filters each option set independently
- **Step 8 (DESeq2)** runs DE for **every** option set, not just the first
- **Step 9 (compare)** compares across all 4-tuple analysis units
- **DE output paths** include the count option: `<trim>/deseq2/<mapper>/<mapper_opt>/<count_opt>/<contrast>/`
- **The report** shows results for every analysis unit with clear labelling

After comparison, the pipeline writes two explicit selected branches:

- `selected_count_comparison.tsv` -- used for featureCounts/DE interpretation
- `selected_visualisation.tsv` -- used for selected-only BigWig generation

For backward compatibility, `selected_analysis.tsv` is still written and mirrors the count-comparison selection.

---

## Command-Line Options

### Data pointer

| Flag | Default | Purpose |
|------|---------|---------|
| `--dataset GSEXXXXX` | *(none)* | Dataset accession. Derives `fastq_dir` and `metadata_csv` |
| `--data-root /path` | `/data` | Root directory for datasets and indices |
| `--metadata /path/file.csv` | *(derived)* | Override metadata CSV path |
| `--species mouse\|human` | *(from config)* | Derives reference genome paths |
| `--outdir path` | `results/<run_id>` | Override output directory |
| `--strict` | off | Fail if metadata run IDs have no matching FASTQs |

### Run control

| Flag | Purpose |
|------|---------|
| `--config path` | Path to config YAML (default: `config/config.yaml`) |
| `--run-id ID` | Override auto-generated run ID |
| `--profile primary\|primary_bw\|full` | Execution profile (sets step list) |
| `--steps 0 1 2 5` | Run only specific steps (overrides `--profile`) |
| `--methods none cutadapt` | Override trimming methods |
| `--subset day3` | Use a specific subset filter from config |
| `--threads N` | Override thread count |
| `--force` | Regenerate existing outputs |

### Examples

```bash
# Primary analysis (recommended):
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --profile primary

# Full pipeline:
./py scripts/run_pipeline.py --dataset GSE48519 --species mouse run --profile full

# Just QC steps:
./py scripts/run_pipeline.py --dataset GSE48519 run --steps 0 1 2 3 4

# Only the 'none' trimming method:
./py scripts/run_pipeline.py --dataset GSE48519 run --methods none

# Strandedness test (before main analysis):
./py scripts/run_strand_test.py --config config/config.yaml \
    --dataset GSE48519 --species mouse
```

---

## Strandedness Verification

Wrong strandedness is the most common cause of low assigned-read rates. The pipeline:

1. **Warns automatically** if the mean assigned rate drops below 30% during featureCounts
2. Provides a dedicated **strandedness test** command that runs featureCounts with -s 0, -s 1, and -s 2 on a single mapping unit and recommends the best setting:

```bash
./py scripts/run_strand_test.py --config config/config.yaml \
    --dataset GSE48519 --species mouse --run-id my_run
```

This requires step 5 (mapping) to have completed first.

---

## Configuration

All parameters live in `config/config.yaml`. Key sections:

### featureCounts option sets

Multiple counting stringency levels run in parallel. DE is run for each:

```yaml
featurecounts:
  backend: "featurecounts"
  strandedness: 2            # 0=unstranded, 1=stranded, 2=reverse-stranded
  option_sets:
    default:
      B: false
      P: false
      C: false
      Q: 0
    strict:
      B: true
      P: true
      C: true
      Q: 10
    stringent:
      B: true
      P: true
      C: true
      Q: 255
    multimapper:
      M: true
      fraction: true
```

### Selected analysis branches

The pipeline auto-selects separate branches for interpretation and visualisation. Optional manual overrides:

```yaml
selected_count_comparison:
  trim_method: "none"
  mapper: "star"
  mapper_option_set: "default"
  count_option_set: "strict"

selected_visualisation:
  trim_method: "none"
  mapper: "star"
  mapper_option_set: "strict_unique"
  count_option_set: "default"
```

### BigWig

BigWig generation is off by default in the `primary` profile. Enable and configure:

```yaml
bigwig:
  enabled: false
  mode: "selected_only"   # "selected_only" | "all_units"
  normalization: "CPM"
  mapq_filter: 255
  use_deseq2_sizefactors: true
```

### Comparisons (contrasts)

```yaml
comparisons:
  - name: "tet1_vs_wt"
    numerator: "tet1"
    denominator: "wt"
deseq2:
  reference_level: "wt"
  fdr_threshold: 0.05
```

---

## Output Structure

```
results/<run_id>/
  samples.tsv                                # Design table
  mapping_summary.tsv                        # Mapping rates
  featurecounts_summary.tsv                  # Assigned reads + Assigned_pct
  filtering_summary.tsv                      # Genes filtered per option set
  redundancy_summary.tsv                     # Duplicate clean-matrix relationships
  de_summary.tsv                             # DEG summary + status per branch/contrast
  selected_count_comparison.tsv              # Selected interpretation branch
  selected_visualisation.tsv                 # Selected BigWig branch
  selected_analysis.tsv                      # Legacy compatibility selection file
  reports/
    report.md / report.html                  # Combined report with health checks
    method_comparison.md                     # Full comparison report
  <trim_method>/
    featurecounts/<mapper>/<mapper_opt>/
      counts_<opt>.tsv                       # Raw featureCounts output
      count_matrix_<opt>.tsv                 # Cleaned gene x sample matrix
      clean_matrix_<opt>.tsv                 # Filtered matrix (DE input)
    deseq2/<mapper>/<mapper_opt>/<count_opt>/<contrast>/
      de_all.tsv                             # All genes with log2FC, padj
      de_significant.tsv                     # FDR-significant genes
      normalized_counts.tsv                  # DESeq2-normalised counts
      size_factors.tsv                       # DESeq2 size factors
      pca.pdf, ma_plot.pdf, volcano.pdf      # Plots
      top_genes_for_enrichment.txt           # Gene lists for enrichment tools
    bigwig/<mapper>/<mapper_opt>/
      <sample>.CPM.MAPQ255.bw                # CPM-normalised track
      <sample>.DESeq2scaled.MAPQ255.bw       # Size-factor-scaled track
```

---

## Reliability Model

- **Run-scoped work paths:** intermediates are isolated under `work/<run_id>/...` (for example `fastq_links`, `trimmed`) to prevent cross-run interference.
- **Shared reusable cache:** reference-derived filtered GTF files are cached under `cache/filtered_gtf/` and reused across runs.
- **Concurrent-safe cache build:** filtered GTF generation uses lock + atomic rename so two runs do not corrupt shared cache entries.
- **BAM safety checks:** mapping validates BAM integrity before indexing, and indexing failures are wrapped with actionable context (run/sample/branch/path).
- **Downstream BAM guards:** BAM-consuming steps fail early when BAMs are missing/unindexed rather than producing ambiguous later errors.

---

## Interpreting the Results

### DE results (`de_all.tsv`)

| Column | Meaning |
|--------|---------|
| `baseMean` | Average normalised expression across all samples |
| `log2FoldChange` | Effect size (positive = up in numerator) |
| `padj` | FDR-adjusted p-value. **Use this for significance.** |

### Pipeline health

The report includes a **Pipeline Health** section showing:
- Whether each expected output file exists
- Branch-level consistency between counting/filtering and DE summaries
- Clear distinction between `MISSING`, `REDUNDANT_SKIPPED`, `COMPLETED_ZERO_DEG`, `MALFORMED`, and passing rows

### Selected branches

The report highlights both selected roles:

- **Count-comparison branch** for featureCounts/DE interpretation
- **Visualisation branch** for selected-only BigWig generation

The legacy `selected_analysis.tsv` is retained for compatibility and mirrors the count-comparison selection.

---

## Functional Enrichment Analysis

The pipeline exports enrichment-ready gene lists per contrast:
- `top_genes_for_enrichment.txt` -- all significant DEGs
- `top_up_genes_for_enrichment.txt` -- upregulated only
- `top_down_genes_for_enrichment.txt` -- downregulated only
- `top200_*` variants for tools with upload limits

Recommended tools: **g:Profiler**, **Enrichr**, **iDEP**.

---

## Repository Structure

```
py                               # Interpreter wrapper
Makefile                         # Convenience targets
README.md                        # This file
config/
  config.yaml                    # Mouse GSE48519 configuration
  config_template.yaml           # Blank template for new datasets
src/
  utils.py                       # Logging, subprocess, config loading
  metadata.py                    # CSV parsing, design table, symlinks
  reporting.py                   # Markdown/HTML report builder
  analysis_unit.py               # Shared analysis-unit helpers, path resolution
scripts/
  run_pipeline.py                # Orchestrator (profiles, step ordering)
  run_strand_test.py             # Strandedness verification helper
  inspect_metadata.py            # Metadata CSV inspector
  00_validate_env.py             # Environment validation
  01_prepare_samples.py          # Metadata parsing + sample prep
  02_trim_reads.py               # Read trimming
  03_qc_fastqc.py                # FastQC
  04_multiqc.py                  # MultiQC
  05_map_star.py                 # Read mapping (STAR + optional HISAT2)
  06_bigwig.py                   # Coverage tracks (runs after DE)
  07_featurecounts.py            # Read counting (all option sets)
  08_filter_matrix.py            # Low-expression filtering
  09_deseq2.py                   # DE analysis (all option sets)
  10_compare_methods.py          # Method comparison + branch selection
  11_make_report.py              # Report generation + health checks
tests/
  test_metadata.py               # Metadata/design tests
  test_naming.py                 # Naming/symlink tests
  test_analysis_unit.py          # Path resolution, analysis unit tests
  test_path_integrity.py         # Cross-step path chain tests
  test_profiles.py               # Execution profile and step-order tests
  test_reliability_paths.py      # Run isolation/cache/BAM reliability tests
```

---

## Dependencies

**Python packages**: pydeseq2, pandas, numpy, scipy, matplotlib, scikit-learn, pysam, HTSeq, pyyaml, markdown

**System tools**: STAR, FastQC, MultiQC, cutadapt, featureCounts, samtools, deeptools (bamCoverage)

```bash
# Conda (includes all tools):
conda env create -f env/environment.yml

# Or pip (Python packages only):
pip install -r env/requirements.txt
```

---

## Troubleshooting

- **Low assigned rate (<30%)**: Almost certainly wrong strandedness. Run `run_strand_test.py`.
- **`FATAL: Missing packages`**: Check your Python environment. Use `./py -m pip list`.
- **`FATAL: path(s) not found`**: Check `--dataset` and `--species` arguments.
- **featureCounts paired-end error**: Verify `data.layout` in config matches actual library type.
- **BAM integrity/indexing failure**: check storage/I/O health, remove the failed run directory, and rerun the same command. Mapping now fails early with run/sample/branch context.
- **BigWig step fails**: BigWig now runs after DE. Use `--profile primary_bw` or `--profile full`.
