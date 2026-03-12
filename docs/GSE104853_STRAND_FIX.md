# GSE104853 investigation: 0-DEG result

## Summary

The pipeline produced 0 significant DEGs (FDR < 0.05) across all branches for GSE104853
(miR-125a-5p overexpression vs control in THP-1 cells). Two hypotheses were tested:

1. **Wrong strandedness** -- ruled out. Strand 0 (unstranded) is correct.
2. **Counting from unfiltered vs filtered BAMs** -- ruled out. Counts are identical.

The 0-DEG result is genuine for this dataset with n=3 per group. The pipeline is
analytically correct. Key genes from the paper (IL10RA, CD163) show the expected
direction of effect but do not reach significance.

## What was tested

### 1. Strandedness (tested 2026-03-12)

| setting            | assigned % | genes after filter | DEGs padj<0.05 | IL10RA log2FC/padj | CD163 log2FC/padj |
|--------------------|------------|---------------------|----------------|---------------------|-------------------|
| strand 0 (baseline)| 42.2%      | 6694                | 0              | -0.798 / 1.0        | -1.188 / 1.0      |
| strand 1           | 22.3%      | 4660                | 0              | -0.601 / 1.0        | (filtered out)     |
| strand 2           | 21.8%      | 4693                | 0              | -1.057 / 1.0        | (filtered out)     |

**Conclusion:** Strand 0 has ~2x the assignment rate of 1 or 2. The library is unstranded.
Strandedness was not the cause of 0 DEGs.

### 2. Filtered BAMs (tested 2026-03-12)

| BAM type             | genes after filter | DEGs padj<0.05 | baseMean VPS13D | baseMean RLF  |
|----------------------|---------------------|----------------|-----------------|---------------|
| unfiltered (default) | 6694                | 0              | 14.503910       | 130.129394    |
| filtered (.filtered) | 6694                | 0              | 14.503910       | 130.129394    |

**Conclusion:** Counts are **identical**. featureCounts in default mode (no `-M` flag) already
ignores multi-mapped reads, which is exactly what the step-5 BAM filter removes. The filtered
BAMs contain only uniquely-mapped reads, but featureCounts was already only counting those.
No effect on DE results.

## Why 0 DEGs is the correct pipeline output

- **Sample size:** n=3 per group. With high within-group variance, FDR correction
  pushes all padj toward 1.0. The minimum padj across 6694 genes is ~0.10 (MT-CYB).
- **Effect direction is correct:** IL10RA log2FC = -0.80 (paper: decreased), CD163
  log2FC = -1.19 (paper: reduced). Both are in the expected direction but not
  individually significant.
- **The paper used edgeR + hg19 + Tophat**, not DESeq2 + hg38 + STAR. Differences in
  reference genome, aligner, and DE method may produce different sensitivity. The paper
  also reported pathway-level enrichment (cytokine-cytokine receptor interactions), not
  necessarily large numbers of individual DEGs.
- **The pipeline is not wrong.** It correctly identifies no individually significant
  genes at FDR < 0.05 with the current data and methods. The biological signal is real
  but below the detection threshold for this sample size.

## What was changed in the codebase

1. **Config `config_GSE104853.yaml`:**
   - `featurecounts.strandedness` confirmed at `0` (with strand-test comment).
   - `featurecounts.use_filtered_bam` added (default `false`, tested `true`, no effect).
   - Removed duplicate `extra_args` under `multimapper` option set.

2. **Code `scripts/07_featurecounts.py`:**
   - `build_mapping_units()` accepts `use_filtered_bam` parameter (reads `filtered_bam_path`
     from mapping_summary when true).
   - Column-name cleaning regex updated to handle `.filtered.bam` suffix.

3. **Added `scripts/run_strand_test.py`:** empirical strandedness test (steps 7-9 with
   strand 0/1/2, comparison table).

4. **Added strand-test configs:** `config_GSE104853_strand1.yaml`,
   `config_GSE104853_strand2.yaml`.

5. **Config template:** `use_filtered_bam: false` added to `config_template.yaml`.

## Recommendations for reporting

- Report 0 DEGs at FDR < 0.05 as the primary result.
- Note that IL10RA (log2FC = -0.80) and CD163 (log2FC = -1.19) trend in the direction
  reported in the paper, but are not individually significant.
- If pathway analysis is needed, use a rank-based method (e.g. GSEA on the full ranked
  gene list by stat or log2FC) rather than relying on a DEG cutoff.
- The pipeline settings (strandedness, BAM type, counting method) have been validated
  and are analytically correct for this dataset.
