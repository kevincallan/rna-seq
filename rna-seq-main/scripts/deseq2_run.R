#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 analysis script for the RNA-seq pipeline
# =============================================================================
# Usage:
#   Rscript deseq2_run.R \
#     --counts count_matrix.tsv \
#     --samples sample_description.txt \
#     --contrast_name tet1_vs_wt \
#     --numerator tet1 \
#     --denominator wt \
#     --fdr 0.05 \
#     --lfc 0.0 \
#     --outdir results/method/deseq2/tet1_vs_wt \
#     --reference_level wt
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  return(args[idx + 1])
}

counts_file     <- parse_arg("--counts")
samples_file    <- parse_arg("--samples")
contrast_name   <- parse_arg("--contrast_name", "contrast")
numerator       <- parse_arg("--numerator")
denominator     <- parse_arg("--denominator")
fdr_threshold   <- as.numeric(parse_arg("--fdr", "0.05"))
lfc_threshold   <- as.numeric(parse_arg("--lfc", "0.0"))
outdir          <- parse_arg("--outdir", "deseq2_out")
reference_level <- parse_arg("--reference_level", denominator)

cat(sprintf("\n=== DESeq2 Analysis: %s ===\n", contrast_name))
cat(sprintf("Counts:    %s\n", counts_file))
cat(sprintf("Samples:   %s\n", samples_file))
cat(sprintf("Contrast:  %s vs %s\n", numerator, denominator))
cat(sprintf("FDR:       %s\n", fdr_threshold))
cat(sprintf("Output:    %s\n\n", outdir))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
count_data <- read.delim(counts_file, row.names = 1, check.names = FALSE)
sample_info <- read.delim(samples_file, row.names = 1, check.names = FALSE)

# Ensure column names in count matrix match row names in sample info
common <- intersect(colnames(count_data), rownames(sample_info))
if (length(common) == 0) {
  stop("No matching sample names between count matrix and sample description!")
}

count_data  <- count_data[, common, drop = FALSE]
sample_info <- sample_info[common, , drop = FALSE]

cat(sprintf("Samples matched: %d\n", length(common)))
cat(sprintf("Genes in matrix: %d\n", nrow(count_data)))

# Ensure counts are integers
count_data <- round(count_data)

# Set reference level for condition factor
sample_info$condition <- factor(sample_info$condition)
if (reference_level %in% levels(sample_info$condition)) {
  sample_info$condition <- relevel(sample_info$condition, ref = reference_level)
}

# Filter to only conditions in the contrast
keep_samples <- sample_info$condition %in% c(numerator, denominator)
if (sum(keep_samples) < 2) {
  stop(sprintf("Not enough samples for contrast: %s vs %s", numerator, denominator))
}
count_data  <- count_data[, keep_samples, drop = FALSE]
sample_info <- sample_info[keep_samples, , drop = FALSE]
sample_info$condition <- droplevels(sample_info$condition)

cat(sprintf("Samples in contrast: %d\n", ncol(count_data)))

# ---------------------------------------------------------------------------
# Run DESeq2
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = sample_info,
  design    = ~ condition
)

dds <- DESeq(dds)

cat("\nDESeq2 analysis complete.\n")
cat(sprintf("Contrast: %s\n", paste(resultsNames(dds), collapse = ", ")))

# Results
res <- results(dds, contrast = c("condition", numerator, denominator),
               alpha = fdr_threshold)

cat(sprintf("\nSummary for %s:\n", contrast_name))
summary(res)

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

# All results
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, file = file.path(outdir, "de_all.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Significant results
if (lfc_threshold > 0) {
  sig <- res_df[!is.na(res_df$padj) &
                res_df$padj < fdr_threshold &
                abs(res_df$log2FoldChange) > lfc_threshold, ]
} else {
  sig <- res_df[!is.na(res_df$padj) & res_df$padj < fdr_threshold, ]
}
write.table(sig, file = file.path(outdir, "de_significant.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat(sprintf("\nSignificant DEGs (FDR < %s): %d\n", fdr_threshold, nrow(sig)))

# Size factors
sf <- sizeFactors(dds)
sf_df <- data.frame(sample = names(sf), size_factor = sf)
write.table(sf_df, file = file.path(outdir, "size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also write size factors to the parent deseq2 directory for BigWig step
parent_sf <- file.path(dirname(outdir), "size_factors.tsv")
write.table(sf_df, file = parent_sf,
            sep = "\t", quote = FALSE, row.names = FALSE)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file = file.path(outdir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# PCA plot
tryCatch({
  vsd <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pct_var <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle(paste("PCA -", contrast_name)) +
    theme_minimal()

  ggsave(file.path(outdir, "pca.pdf"), p, width = 8, height = 6)
  cat("PCA plot saved.\n")
}, error = function(e) {
  cat(sprintf("Warning: PCA plot failed: %s\n", e$message))
})

# MA plot
tryCatch({
  pdf(file.path(outdir, "ma_plot.pdf"), width = 8, height = 6)
  plotMA(res, main = paste("MA Plot -", contrast_name), ylim = c(-5, 5))
  dev.off()
  cat("MA plot saved.\n")
}, error = function(e) {
  cat(sprintf("Warning: MA plot failed: %s\n", e$message))
})

# ---------------------------------------------------------------------------
# Session info
# ---------------------------------------------------------------------------
sink(file.path(outdir, "session_info.txt"))
cat(sprintf("DESeq2 analysis: %s\n", contrast_name))
cat(sprintf("Date: %s\n\n", Sys.time()))
sessionInfo()
sink()

cat(sprintf("\nAll outputs written to: %s\n", outdir))
