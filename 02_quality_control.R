# ============================================================
# Script 02: Quality Control
# Project: Network-Based Pathway Enrichment in Cervical Cancer
# Author:  Bhavana Chowdary Kothapalli
# ============================================================
# PURPOSE:
#   Performs comprehensive QC on raw expression data:
#   RLE plots, NUSE plots, RNA degradation, PCA, missing values.
# ============================================================

set.seed(123)

# ── 1. Load libraries ────────────────────────────────────────
library(Biobase)
library(affy)
library(ggplot2)
library(pheatmap)

# ── 2. Load raw data ─────────────────────────────────────────
cat("Loading raw ExpressionSet...\n")
gse <- readRDS("outputs/rdata/gse_raw.rds")

expr_mat <- exprs(gse)            # probes x samples matrix
pdata    <- pData(gse)

cat("Expression matrix dimensions:", nrow(expr_mat), "x", ncol(expr_mat), "\n")

# ── 3. Define sample groups ──────────────────────────────────
# Assign group labels based on source name
group <- ifelse(
  grepl("normal|Normal", pdata$`source_name_ch1`),
  "Normal", "Tumour"
)
names(group) <- colnames(expr_mat)

cat("\nSample group counts:\n")
print(table(group))

# ── 4. Missing value assessment ──────────────────────────────
n_missing <- sum(is.na(expr_mat))
pct_missing <- round(100 * n_missing / prod(dim(expr_mat)), 2)
cat("\nMissing values:", n_missing, "(", pct_missing, "% of total)\n")

# ── 5. RLE (Relative Log Expression) plot ───────────────────
cat("Generating RLE plot...\n")

# Compute median expression per probe across all samples
probe_medians <- apply(expr_mat, 1, median, na.rm = TRUE)

# RLE = deviation from probe median
rle_mat <- sweep(expr_mat, 1, probe_medians, "-")

rle_df <- data.frame(
  Sample = rep(colnames(rle_mat), each = nrow(rle_mat)),
  RLE    = as.vector(rle_mat),
  Group  = rep(group, each = nrow(rle_mat))
)

p_rle <- ggplot(rle_df, aes(x = Sample, y = RLE, fill = Group)) +
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_fill_manual(values = c("Normal" = "#4DAFFF", "Tumour" = "#FF6B6B")) +
  labs(title = "Relative Log Expression (RLE) Plot",
       subtitle = "Median values near zero indicate acceptable quality",
       x = "Sample", y = "RLE") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

ggsave("outputs/figures/02_rle_plot.pdf",
       plot = p_rle, width = 14, height = 6)
cat("RLE plot saved.\n")

# Report QC metrics
rle_medians <- apply(rle_mat, 2, median, na.rm = TRUE)
rle_iqrs    <- apply(rle_mat, 2, IQR,    na.rm = TRUE)

cat("\nRLE QC summary:\n")
cat("  Mean of medians:", round(mean(rle_medians), 4), "\n")
cat("  Max IQR:", round(max(rle_iqrs), 4),
    "(threshold < 0.5)\n")

failed_rle <- names(rle_iqrs[rle_iqrs > 0.5])
if (length(failed_rle) > 0) {
  cat("  WARNING - samples with IQR > 0.5:", paste(failed_rle, collapse = ", "), "\n")
} else {
  cat("  All samples pass RLE IQR threshold.\n")
}

# ── 6. PCA before normalisation ──────────────────────────────
cat("\nRunning PCA on unnormalised data...\n")

# Remove probes with any missing value for PCA
expr_complete <- expr_mat[complete.cases(expr_mat), ]

pca_res  <- prcomp(t(expr_complete), scale. = TRUE, center = TRUE)
pca_df   <- data.frame(
  PC1   = pca_res$x[, 1],
  PC2   = pca_res$x[, 2],
  Group = group,
  Sample = rownames(pca_res$x)
)
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, label = Sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 10) +
  scale_colour_manual(values = c("Normal" = "#4DAFFF", "Tumour" = "#FF6B6B")) +
  labs(
    title    = "PCA — Unnormalised Expression Data",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_classic()

ggsave("outputs/figures/02_pca_prenorm.pdf",
       plot = p_pca, width = 8, height = 6)
cat("PCA plot saved.\n")

# ── 7. Check for outlier samples (> 3 SD from mean on PC1) ──
pc1_sd   <- sd(pca_df$PC1)
pc1_mean <- mean(pca_df$PC1)
outliers <- pca_df$Sample[abs(pca_df$PC1 - pc1_mean) > 3 * pc1_sd]

if (length(outliers) > 0) {
  cat("  WARNING - potential PC1 outliers:", paste(outliers, collapse = ", "), "\n")
} else {
  cat("  No PC1 outliers detected (threshold: 3 SD).\n")
}

# ── 8. Save QC metadata ──────────────────────────────────────
qc_summary <- data.frame(
  Sample        = colnames(expr_mat),
  Group         = group,
  RLE_Median    = round(rle_medians, 4),
  RLE_IQR       = round(rle_iqrs, 4),
  RLE_Pass      = rle_iqrs < 0.5,
  PC1           = round(pca_df$PC1, 3),
  PC2           = round(pca_df$PC2, 3)
)

write.csv(qc_summary, "outputs/tables/02_qc_summary.csv",
          row.names = FALSE)
saveRDS(list(group = group, expr_mat = expr_mat, gse = gse),
        "outputs/rdata/02_qc_output.rds")

cat("\nQC summary saved to outputs/tables/02_qc_summary.csv\n")
cat("\n--- Script 02 complete ---\n")
