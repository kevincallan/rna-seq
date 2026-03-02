# RNA-seq Analysis Pipeline

A config-driven, end-to-end RNA-seq pipeline for differential expression analysis. Designed for exam-day use on the university JupyterHub server. Produces QC reports, count matrices, DE results with statistical tables, publication-quality figures, and gene lists ready for functional enrichment analysis.

---

## Exam Day Quick Start

```bash
cd ~/rna-seq-pipeline

# One command -- change only the dataset ID:
./py scripts/run_pipeline.py --dataset GSEXXXXX run
```

`./py` is a wrapper that ensures the correct JupyterHub interpreter (`/opt/jupyterhub/bin/python3`) is used. The pipeline will:

1. Verify the interpreter and all required packages
2. Derive data paths from `--dataset` (FASTQs and metadata in `/data/GSEXXXXX/`)
3. Validate that files exist and FASTQ filenames match metadata run IDs
4. Run all 12 steps: QC, trimming, mapping, counting, filtering, DE analysis, plots, and report generation

On exam day, only the dataset accession changes. Everything else (references, parameters, tool paths) is already configured.

---

## What the Pipeline Does

This pipeline takes raw RNA-seq FASTQ files and produces a complete differential expression analysis:

```
Raw FASTQs
  --> Quality control (FastQC / MultiQC)
  --> Adapter trimming (cutadapt, with untrimmed baseline)
  --> Read mapping to reference genome (STAR and/or HISAT2 option sets)
  --> Read counting per gene (featureCounts / HTSeq)
  --> Low-expression gene filtering
  --> Differential expression analysis (PyDESeq2)
  --> Figures: PCA, MA plot, volcano plot
  --> Gene lists for functional enrichment
  --> Method comparison report
  --> Final summary report
```

All steps are automated and config-driven. The pipeline compares multiple trimming approaches and featureCounts parameter sets side-by-side, which demonstrates how preprocessing choices affect downstream results.

---

## Pipeline Stages

| Step | Script | What it does | Key outputs |
|------|--------|-------------|-------------|
| 0 | `00_validate_env.py` | Checks all tools and paths exist, creates output directories | `run_manifest.json` |
| 1 | `01_prepare_samples.py` | Parses metadata CSV, applies subset filters, creates FASTQ symlinks with standardised names | `samples.tsv`, `sample_description.txt` |
| 2 | `02_trim_reads.py` | Runs each enabled trimming method (none, cutadapt) on all samples | Trimmed FASTQs, `trimming_summary.tsv` |
| 3 | `03_qc_fastqc.py` | Runs FastQC on raw and trimmed reads | FastQC HTML reports per sample |
| 4 | `04_multiqc.py` | Aggregates FastQC reports into MultiQC summaries | MultiQC HTML reports |
| 5 | `05_map_star.py` | Aligns reads to the reference genome using STAR, indexes BAMs | Sorted BAMs, `mapping_summary.tsv` |
| 6 | `06_bigwig.py` | Generates coverage tracks (CPM and DESeq2-scaled) | BigWig files for genome browser |
| 7 | `07_featurecounts.py` | Counts reads per gene using multiple parameter sets | Count matrices per option set |
| 8 | `08_filter_matrix.py` | Removes lowly-expressed genes | Filtered count matrices, `filtering_summary.tsv` |
| 9 | `09_deseq2.py` | Runs PyDESeq2 for each contrast: normalisation, Wald test, plots | DE tables, PCA, MA, volcano plots, gene lists |
| 10 | `10_compare_methods.py` | Compares trimming methods: mapping rates, DEG overlap, correlations | `method_comparison.md` |
| 11 | `11_make_report.py` | Generates a combined Markdown + HTML report | `report.md`, `report.html` |

---

## Command-Line Options

### Data pointer (change on exam day)

| Flag | Default | Purpose |
|------|---------|---------|
| `--dataset GSEXXXXX` | *(none)* | Dataset accession. Derives `fastq_dir = /data/GSEXXXXX` and `metadata_csv = /data/GSEXXXXX/metadata.csv` |
| `--data-root /path` | `/data` | Root directory where datasets and indices live |
| `--metadata /path/file.csv` | *(derived)* | Override metadata CSV path if not in the standard location |
| `--species mouse\|human` | *(from config)* | Optional. Derives reference genome paths (`/data/indices/mm39/` or `/data/indices/hg38/`). Only needed if switching organism. |
| `--outdir path` | `results/<run_id>` | Override output directory |
| `--strict` | off | Fail if any metadata run IDs have no matching FASTQ files |

### Run control

| Flag | Purpose |
|------|---------|
| `--config path` | Path to config YAML (default: `config/config.yaml`) |
| `--run-id ID` | Override the auto-generated run ID |
| `run --steps 0 1 2 5` | Run only specific steps |
| `run --methods none cutadapt` | Override which trimming methods to run |
| `run --subset day3` | Use a specific subset filter from config |
| `run --threads N` | Override thread count |

### Examples

```bash
# Full pipeline with dataset pointer:
./py scripts/run_pipeline.py --dataset GSE48519 run

# Just QC steps:
./py scripts/run_pipeline.py --dataset GSE48519 run --steps 0 1 2 3 4

# Run through mapping only:
./py scripts/run_pipeline.py --dataset GSE48519 run --steps 0 1 2 3 4 5

# Only the 'none' trimming method (fastest):
./py scripts/run_pipeline.py --dataset GSE48519 run --methods none

# Use Make shortcuts:
make run
make qc       # steps 0-4
make map      # steps 0-5
make count    # steps 0-8
make de       # steps 0-9
```

---

## Configuration Parameters

All parameters live in `config/config.yaml`. Key sections you may want to adjust:

### Column mapping

Maps your metadata CSV columns to pipeline concepts:

```yaml
column_mapping:
  run_id_col: "Run"                    # Column with SRR/ERR IDs
  condition_cols: ["Genotype"]         # Column(s) defining experimental conditions
  condition_map:                       # Map raw values to short names
    "wild-type": "wt"
    "tet1-/-": "tet1"
```

### Subset filters

Select which samples to analyse (e.g., only day 3 differentiation):

```yaml
subset_filters:
  default:
    "Assay Type": "RNA-Seq"
    "differentiation": ""
  day3:
    "Assay Type": "RNA-Seq"
    "differentiation": "day 3"
active_subset: "default"
```

### Comparisons (contrasts)

Define which conditions to compare:

```yaml
comparisons:
  - name: "tet1_vs_wt"
    numerator: "tet1"
    denominator: "wt"
deseq2:
  reference_level: "wt"       # Control condition (denominator)
  fdr_threshold: 0.05          # Significance threshold
  lfc_threshold: 0.0           # log2FC cutoff (0 = no filter)
```

### featureCounts option sets

Three counting stringency levels run in parallel to show parameter effects:

```yaml
featurecounts:
  backend: "featurecounts"     # or "htseq"
  option_sets:
    default:                   # No special flags
      B: false
      P: false
      C: false
      Q: 0
    strict:                    # Both ends mapped, proper pairs, no chimeras
      B: true
      P: true
      C: true
      Q: 10
    stringent:                 # Unique mappers only (MAPQ=255)
      B: true
      P: true
      C: true
      Q: 255
```

### Trimming methods

Toggle trimming approaches. Each enabled method produces independent downstream results:

```yaml
trimming:
  none:
    enabled: true              # Untrimmed baseline
  cutadapt:
    enabled: true              # Adapter + quality trimming
    quality: 20
    min_length: 25
```

---

## Output Structure and Key Files

```
results/<run_id>/
  run_manifest.json                          # Reproducibility: config hash, tool versions
  samples.tsv                                # Design table (sample -> condition -> replicate)
  mapping_summary.tsv                        # Mapping rates per trim/mapping approach/sample
  featurecounts_summary.tsv                  # Assigned reads per trim/mapping/count option set
  filtering_summary.tsv                      # Genes before/after filtering (all dimensions)
  de_summary.tsv                             # DEG counts per trim/mapping/contrast
  reports/
    report.md / report.html                  # Combined final report
    method_comparison.md                     # Trimming + mapper + parameter comparison
    multiqc_raw/multiqc_report.html          # Raw reads QC
    multiqc_<method>/multiqc_report.html     # Per-method QC
  <method>/                                  # e.g. none/, cutadapt/
    mapping/<mapper>/<mapper_option>/<sample>...
    featurecounts/<mapper>/<mapper_option>/count_matrix_<set>.tsv
    featurecounts/<mapper>/<mapper_option>/clean_matrix_<set>.tsv
    deseq2/<mapper>/<mapper_option>/<contrast>/
      de_all.tsv                             # All genes: log2FC, pvalue, padj
      de_significant.tsv                     # FDR-significant genes only
      normalized_counts.tsv                  # DESeq2-normalised count matrix
      size_factors.tsv                       # DESeq2 size factors
      top_genes_for_enrichment.txt           # All significant genes
      top_up_genes_for_enrichment.txt        # Significant upregulated genes only
      top_down_genes_for_enrichment.txt      # Significant downregulated genes only
      top200_genes_for_enrichment.txt        # Capped list for web tools
      top200_up_genes_for_enrichment.txt     # Capped up list
      top200_down_genes_for_enrichment.txt   # Capped down list
      pca.pdf                                # PCA plot
      ma_plot.pdf                            # MA plot
      volcano.pdf                            # Volcano plot
      session_info.txt                       # Package versions
```

---

## Interpreting the Results

### DE results table (`de_all.tsv`, `de_significant.tsv`)

Each row is a gene. Key columns:

| Column | Meaning |
|--------|---------|
| `baseMean` | Average normalised expression across all samples |
| `log2FoldChange` | Effect size: positive = upregulated in numerator, negative = downregulated |
| `lfcSE` | Standard error of the log2FC estimate |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value from the Wald test |
| `padj` | Benjamini-Hochberg adjusted p-value (FDR). **Use this for significance.** |

A gene is significantly DE when `padj < 0.05` (configurable via `deseq2.fdr_threshold`).

### PCA plot (`pca.pdf`)

Shows the first two principal components of variance-stabilised counts. Samples should cluster by condition. If they don't, check for batch effects, mislabelling, or outlier samples. The percentage on each axis indicates how much variance that component explains.

### MA plot (`ma_plot.pdf`)

Plots log2 fold change (y) against mean expression (x). Red points are significant genes. Most genes cluster near zero fold change. The spread of red points shows the magnitude and distribution of differential expression.

### Volcano plot (`volcano.pdf`)

Plots log2 fold change (x) against statistical significance as -log10(padj) (y). Points in the upper corners are both highly significant and have large effect sizes -- these are the most biologically interesting candidates. The dashed horizontal line marks the FDR threshold.

### Method comparison (`method_comparison.md`)

Compares analysis units (trim method + mapper + mapper option) on:
- **Mapping rates**: Higher is generally better, but very high rates may indicate reference contamination
- **Mapper option impact**: Mean unique mapping percentage across mapper settings
- **Assigned reads**: Fraction of mapped reads that overlap annotated features
- **featureCounts option impact**: Delta vs default counting stringency
- **Count correlations**: Pearson r between normalised counts from different analysis units
- **DEG overlap**: Jaccard-style overlap of significant gene sets
- **Method-sensitive genes**: Genes whose significance flips between units -- these findings are fragile and should be reported cautiously

### What the results mean biologically

- **Large number of DEGs** (hundreds-thousands): Strong transcriptomic response to the condition
- **Few DEGs**: Subtle or no effect; or insufficient replicates/sequencing depth
- **Asymmetric up/down**: Condition may activate or repress specific pathways
- **PCA separation**: Conditions have distinct expression profiles
- **High method agreement**: Results are robust to preprocessing choices

---

## Functional Enrichment Analysis

The pipeline exports multiple enrichment-ready lists (one gene ID per line):
- `top_genes_for_enrichment.txt` (all significant DEGs)
- `top_up_genes_for_enrichment.txt` (upregulated only)
- `top_down_genes_for_enrichment.txt` (downregulated only)
- `top200_*` variants for tools with practical upload-size limits

### Recommended tools

1. **g:Profiler** (https://biit.cs.ut.ee/gprofiler/gost)
   - Paste gene list, select organism
   - Returns enriched GO terms (Biological Process, Molecular Function, Cellular Component), KEGG pathways, Reactome
   - Produces publication-quality plots

2. **Enrichr** (https://maayanlab.cloud/Enrichr/)
   - Paste gene list
   - Searches many gene-set libraries (GO, KEGG, WikiPathways, transcription factor targets, etc.)
   - Combined score ranks enrichment by significance and effect size

3. **iDEP** (https://bioinformatics.sdstate.edu/idep/)
   - Upload the full count matrix (`normalized_counts.tsv`) or DE table (`de_all.tsv`)
   - Provides PCA, clustering, enrichment, and pathway analysis in one interface
   - Useful for exploratory analysis beyond what this pipeline produces

### Enrichment workflow (recommended and reproducible)

1. Choose the exact analysis unit path:
   `results/<run_id>/<trim_method>/deseq2/<mapper>/<mapper_option>/<contrast>/`
2. Use **upregulated genes** (`top_up_genes_for_enrichment.txt`) for activated pathways.
3. Use **downregulated genes** (`top_down_genes_for_enrichment.txt`) for repressed pathways.
4. If the tool is slow or rejects long lists, use `top200_*.txt` files.
5. In g:Profiler/Enrichr, set organism explicitly (for example *Mus musculus*).
6. Save top terms + adjusted p-values and include 1-2 pathway findings in the report.
7. For iDEP, upload `normalized_counts.tsv` from the same analysis unit to keep comparisons consistent.

### Interpreting enrichment results

- **GO Biological Process**: What biological processes are affected? (e.g., "cell differentiation", "immune response")
- **KEGG/Reactome pathways**: Which specific signalling or metabolic pathways are perturbed?
- **Significance**: Use adjusted p-values. Many enriched terms may be redundant (parent-child GO relationships)
- **Compare with published results**: Do your enriched pathways match the original paper's findings? Concordance strengthens your analysis; discrepancies should be discussed

---

## Parameter Effects to Explore (Advanced)

The pipeline is designed to demonstrate how parameter choices affect RNA-seq results. These are good discussion points for the report:

### 1. Trimming vs no trimming

Compare `none/` and `cutadapt/` results:
- Does trimming change the mapping rate? (Check `mapping_summary.tsv`)
- Does it change which genes are called DE? (Check `method_comparison.md`)
- Implication: Adapter contamination inflates multi-mapping; trimming improves specificity

### 2. Mapping approach / mapper option

Compare mapper settings (for example STAR default vs STAR strict_unique):
- Stricter mapping often raises confidence but can reduce mapped read counts.
- Check downstream impact on assigned reads and DEG overlap.
- Report whether biological conclusions are stable across mapper settings.

### 3. featureCounts stringency

Compare count option sets (default, strict, stringent):
- `strict` requires both mates mapped and proper pairing -- removes ambiguous fragments
- `stringent` (Q=255) keeps only uniquely mapped reads -- most conservative
- Implication: Stricter counting reduces noise but may lose signal for repetitive loci

### 4. FDR threshold sensitivity

The default threshold is 0.05. Discuss what happens at 0.01 or 0.1:
- More stringent = fewer DEGs but higher confidence
- Less stringent = more DEGs but more false positives
- This can be explored post-hoc from `de_all.tsv` without re-running the pipeline

### 5. Replicate count

Using more than the minimum 2 replicates increases statistical power. The pipeline automatically uses all replicates matching the subset filter. More replicates = narrower confidence intervals = more DEGs detected.

---

## Exam-Day Workflow

1. **Receive** your dataset ID (e.g., `GSE12345`), subset assignment, and metadata CSV
2. **Edit** `config/config.yaml`:
   - Update `column_mapping` if the metadata CSV has different column names
   - Update `condition_map` to map raw condition values to short names
   - Update `subset_filters` to define your assigned subset
   - Update `comparisons` to define your contrast (numerator vs denominator)
   - Update `deseq2.reference_level` to set the control condition
3. **Run**:

```bash
cd ~/rna-seq-pipeline
./py scripts/run_pipeline.py --dataset GSE12345 run
```

4. **Collect outputs** from `results/<run_id>/`:
   - Figures for report: `pca.pdf`, `volcano.pdf`, `ma_plot.pdf` (choose best 4)
   - Tables: `de_significant.tsv` (top genes), `mapping_summary.tsv` (QC stats)
   - Gene list: `top_genes_for_enrichment.txt` --> paste into g:Profiler
5. **Run enrichment** analysis in browser (g:Profiler / Enrichr)
6. **Write** 800-1200 word report with 4 figures

### Report structure suggestion

1. **Introduction**: Dataset description, biological question, experimental design
2. **Methods**: Pipeline overview (QC, mapping, counting, DE method, parameters used)
3. **Results**: Key figures (PCA showing sample grouping, volcano plot showing DE landscape), top DE genes table, enrichment findings
4. **Discussion**: Compare with published findings, discuss parameter effects, note limitations

---

## Repository Structure

```
pipeline/
  py                           # Interpreter wrapper (./py ensures correct Python)
  Makefile                     # Convenience targets (make run, make qc, etc.)
  README.md                    # This file
  config/
    config.yaml                # ALL configuration
  env/
    environment.yml            # Conda environment spec
    requirements.txt           # pip requirements
  src/
    utils.py                   # Logging, subprocess, config loading
    metadata.py                # CSV parsing, design table, symlinks
    reporting.py               # Markdown/HTML report builder
  scripts/
    run_pipeline.py            # Orchestrator
    00_validate_env.py         # Environment validation
    01_prepare_samples.py      # Metadata parsing + sample prep
    02_trim_reads.py           # Read trimming
    03_qc_fastqc.py            # FastQC
    04_multiqc.py              # MultiQC
    05_map_star.py             # STAR alignment
    06_bigwig.py               # Coverage tracks
    07_featurecounts.py        # Read counting
    08_filter_matrix.py        # Low-expression filtering
    09_deseq2.py               # Differential expression + plots
    10_compare_methods.py      # Method comparison
    11_make_report.py          # Report generation
    deseq2_run.R               # R DESeq2 backend (optional)
  tests/
    test_metadata.py
    test_naming.py
```

---

## Dependencies

On the university JupyterHub server, all dependencies are pre-installed. The `./py` wrapper ensures the correct environment is used.

**Python packages**: pydeseq2, pandas, numpy, scipy, matplotlib, scikit-learn, pysam, HTSeq, pyyaml, markdown

**System tools**: STAR, FastQC, MultiQC, cutadapt, featureCounts, samtools, deeptools (bamCoverage)

For other environments, install via:

```bash
# Conda (includes all tools):
conda env create -f env/environment.yml

# Or pip (Python packages only; system tools must be installed separately):
./py -m pip install -r env/requirements.txt
```

---

## Running on Google Colab

See [docs/COLAB_COMMANDS.md](docs/COLAB_COMMANDS.md) for Colab-specific instructions, chunked runs, and memory tips. Colab uses `scripts/setup_colab.sh` to install dependencies and configure a chr19-only reference for fast testing.
