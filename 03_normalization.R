# ============================================================
# Script 03: Normalisation
# Project: Network-Based Pathway Enrichment in Cervical Cancer
# Author:  Bhavana Chowdary Kothapalli
# ============================================================
# PURPOSE:
#   Imputes missing values, applies quantile normalisation,
#   maps probes to gene symbols, and saves the normalised matrix.
# ============================================================

set.seed(123)

# ── 1. Load libraries ────────────────────────────────────────
library(preprocessCore)   # normalize.quantiles()
library(hgu133plus2.db)   # Affymetrix annotation
library(AnnotationDbi)
library(ggplot2)

# ── 2. Load QC output ────────────────────────────────────────
cat("Loading QC output...\n")
qc_out   <- readRDS("outputs/rdata/02_qc_output.rds")
expr_mat <- qc_out$expr_mat
group    <- qc_out$group
gse      <- qc_out$gse

cat("Dimensions before processing:", nrow(expr_mat), "x", ncol(expr_mat), "\n")

# ── 3. Missing value imputation ──────────────────────────────
cat("\nImputing missing values using row (probe) means...\n")
n_missing <- sum(is.na(expr_mat))
cat("Missing values before imputation:", n_missing, "\n")

# Replace NA with the mean of non-NA values in each row (probe)
row_means <- rowMeans(expr_mat, na.rm = TRUE)
for (i in seq_len(nrow(expr_mat))) {
  na_idx <- is.na(expr_mat[i, ])
  if (any(na_idx)) {
    expr_mat[i, na_idx] <- row_means[i]
  }
}

cat("Missing values after imputation:", sum(is.na(expr_mat)), "\n")

# ── 4. Quantile normalisation ────────────────────────────────
cat("\nApplying quantile normalisation...\n")

expr_norm <- normalize.quantiles(expr_mat)
dimnames(expr_norm) <- dimnames(expr_mat)   # restore row/col names

cat("Normalisation complete.\n")

# ── 5. Verify normalisation with boxplots ────────────────────
# Sample 6 samples for a readable boxplot
sample_idx <- seq(1, ncol(expr_norm), length.out = min(6, ncol(expr_norm)))
sample_idx <- round(sample_idx)

plot_df_before <- data.frame(
  Value  = as.vector(expr_mat[, sample_idx]),
  Sample = rep(colnames(expr_mat)[sample_idx], each = nrow(expr_mat)),
  Stage  = "Before normalisation"
)
plot_df_after <- data.frame(
  Value  = as.vector(expr_norm[, sample_idx]),
  Sample = rep(colnames(expr_norm)[sample_idx], each = nrow(expr_norm)),
  Stage  = "After normalisation"
)
plot_df <- rbind(plot_df_before, plot_df_after)

p_box <- ggplot(plot_df, aes(x = Sample, y = Value, fill = Stage)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Stage, scales = "free_y") +
  labs(title = "Expression Distribution Before and After Quantile Normalisation",
       x = "Sample", y = "Expression Value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("outputs/figures/03_normalisation_boxplot.pdf",
       plot = p_box, width = 10, height = 6)
cat("Normalisation boxplot saved.\n")

# ── 6. Probe-to-gene annotation ──────────────────────────────
cat("\nMapping probe IDs to gene symbols and Entrez IDs...\n")

probe_ids <- rownames(expr_norm)

# Map to gene symbols
symbols <- mapIds(
  hgu133plus2.db,
  keys      = probe_ids,
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

# Map to Entrez IDs
entrez <- mapIds(
  hgu133plus2.db,
  keys      = probe_ids,
  column    = "ENTREZID",
  keytype   = "PROBEID",
  multiVals = "first"
)

# Summarise annotation coverage
n_mapped <- sum(!is.na(symbols))
cat("Probes mapped to gene symbols:", n_mapped, "/", length(probe_ids),
    "(", round(100 * n_mapped / length(probe_ids), 1), "%)\n")

# ── 7. Collapse to gene-level (mean per gene symbol) ─────────
cat("\nCollapsing probes to gene level (mean per gene symbol)...\n")

annotation_df <- data.frame(
  ProbeID  = probe_ids,
  Symbol   = symbols,
  EntrezID = entrez,
  stringsAsFactors = FALSE
)

# Keep only probes that map to a gene symbol
keep     <- !is.na(annotation_df$Symbol)
ann_keep <- annotation_df[keep, ]
expr_keep <- expr_norm[keep, ]

# Aggregate by mean expression per gene symbol
gene_expr <- rowsum(expr_keep, group = ann_keep$Symbol) /
             as.vector(table(ann_keep$Symbol)[unique(ann_keep$Symbol)])

cat("Gene-level matrix dimensions:", nrow(gene_expr), "x", ncol(gene_expr), "\n")

# Build gene annotation lookup
gene_ann <- ann_keep[!duplicated(ann_keep$Symbol), c("Symbol", "EntrezID")]
rownames(gene_ann) <- gene_ann$Symbol

# ── 8. PCA on normalised gene-level data ─────────────────────
cat("\nRunning post-normalisation PCA...\n")

pca_norm <- prcomp(t(gene_expr), scale. = TRUE, center = TRUE)
pct_var  <- round(100 * pca_norm$sdev^2 / sum(pca_norm$sdev^2), 1)

pca_df <- data.frame(
  PC1    = pca_norm$x[, 1],
  PC2    = pca_norm$x[, 2],
  PC3    = pca_norm$x[, 3],
  Group  = group,
  Sample = rownames(pca_norm$x)
)

p_pca_norm <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, label = Sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 10) +
  scale_colour_manual(values = c("Normal" = "#4DAFFF", "Tumour" = "#FF6B6B")) +
  labs(
    title = "PCA — Normalised Expression Data",
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_classic()

ggsave("outputs/figures/03_pca_postnorm.pdf",
       plot = p_pca_norm, width = 8, height = 6)

cat("Post-normalisation PCA: PC1 explains", pct_var[1], "% of variance.\n")

# ── 9. Save outputs ──────────────────────────────────────────
saveRDS(
  list(
    expr_norm = expr_norm,
    gene_expr = gene_expr,
    gene_ann  = gene_ann,
    group     = group,
    pct_var   = pct_var
  ),
  "outputs/rdata/03_normalised.rds"
)

cat("\nNormalised data saved to outputs/rdata/03_normalised.rds\n")
cat("\n--- Script 03 complete ---\n")
