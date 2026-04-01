# ============================================================
# Script 04: Differential Expression Analysis
# Project: Network-Based Pathway Enrichment in Cervical Cancer
# Author:  Bhavana Chowdary Kothapalli
# ============================================================
# PURPOSE:
#   Uses limma with empirical Bayes moderation to identify DEGs
#   between tumour and normal samples. Applies BH FDR correction.
#   Thresholds: |log2FC| > 1, adj.p < 0.05.
# ============================================================

set.seed(123)

# ── 1. Load libraries ────────────────────────────────────────
library(limma)
library(ggplot2)
library(pheatmap)

# ── 2. Load normalised data ──────────────────────────────────
cat("Loading normalised data...\n")
norm_out  <- readRDS("outputs/rdata/03_normalised.rds")
gene_expr <- norm_out$gene_expr
gene_ann  <- norm_out$gene_ann
group     <- norm_out$group

cat("Gene expression matrix:", nrow(gene_expr), "genes x",
    ncol(gene_expr), "samples\n")

# ── 3. Build design matrix ───────────────────────────────────
group_factor <- factor(group, levels = c("Normal", "Tumour"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Normal", "Tumour")

cat("\nDesign matrix summary:\n")
print(colSums(design))

# ── 4. Fit linear model ──────────────────────────────────────
cat("\nFitting linear model with limma...\n")
fit <- lmFit(gene_expr, design)

# Define contrast: Tumour vs Normal
contrast_mat <- makeContrasts(
  Tumour_vs_Normal = Tumour - Normal,
  levels           = design
)

fit2 <- contrasts.fit(fit, contrast_mat)

# Empirical Bayes moderation of variance estimates
fit2 <- eBayes(fit2)

cat("Model fitting complete.\n")

# ── 5. Extract results ───────────────────────────────────────
all_results <- topTable(
  fit2,
  coef      = "Tumour_vs_Normal",
  number    = Inf,        # return all genes
  adjust    = "BH",       # Benjamini-Hochberg FDR
  sort.by   = "adj.P.Val"
)

# Add gene symbol column
all_results$Symbol <- rownames(all_results)

# Report total significant
sig_results <- all_results[
  all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1, ]

cat("\n--- Differential Expression Summary ---\n")
cat("Total genes tested:           ", nrow(all_results), "\n")
cat("Significant DEGs (adj.p<0.05, |log2FC|>1):", nrow(sig_results), "\n")
cat("  Upregulated (log2FC > 1):  ",
    sum(sig_results$logFC > 1), "\n")
cat("  Downregulated (log2FC < -1):",
    sum(sig_results$logFC < -1), "\n")

# Most extreme genes
top_up   <- sig_results[which.max(sig_results$logFC), ]
top_down <- sig_results[which.min(sig_results$logFC), ]
cat("\nMost upregulated gene:  ", rownames(top_up),
    "  log2FC =", round(top_up$logFC, 3), "\n")
cat("Most downregulated gene:", rownames(top_down),
    " log2FC =", round(top_down$logFC, 3), "\n")

# Add significance labels
all_results$DEG_Status <- "Not significant"
all_results$DEG_Status[
  all_results$adj.P.Val < 0.05 & all_results$logFC >  1] <- "Upregulated"
all_results$DEG_Status[
  all_results$adj.P.Val < 0.05 & all_results$logFC < -1] <- "Downregulated"

# ── 6. Volcano plot ──────────────────────────────────────────
cat("\nGenerating volcano plot...\n")

# Label top 10 up and top 10 down genes
top_genes <- c(
  rownames(sig_results[order(-sig_results$logFC), ])[1:10],
  rownames(sig_results[order( sig_results$logFC), ])[1:10]
)
all_results$Label <- ifelse(
  rownames(all_results) %in% top_genes,
  rownames(all_results), ""
)

p_volcano <- ggplot(all_results,
                    aes(x = logFC, y = -log10(adj.P.Val),
                        colour = DEG_Status, label = Label)) +
  geom_point(alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20,
                            colour = "black", fontface = "italic") +
  geom_vline(xintercept = c(-1, 1),  linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Upregulated"     = "#D73027",
    "Downregulated"   = "#4575B4",
    "Not significant" = "grey70"
  )) +
  labs(
    title    = "Volcano Plot: Cervical Cancer vs Normal",
    subtitle = "Thresholds: |log\u2082FC| > 1 and adj.p < 0.05 (BH)",
    x        = "log\u2082 Fold Change (Tumour vs Normal)",
    y        = "-log\u2081\u2080 (Adjusted p-value)",
    colour   = "Status"
  ) +
  theme_classic(base_size = 12)

ggsave("outputs/figures/04_volcano_plot.pdf",
       plot = p_volcano, width = 10, height = 8)
cat("Volcano plot saved.\n")

# ── 7. Heatmap of top 50 DEGs ────────────────────────────────
cat("Generating DEG heatmap...\n")

top50_genes <- head(rownames(sig_results), 50)
hm_mat      <- gene_expr[top50_genes, ]

# Scale per gene (row)
hm_scaled <- t(scale(t(hm_mat)))

annotation_col <- data.frame(
  Group = factor(group, levels = c("Normal", "Tumour")),
  row.names = colnames(hm_mat)
)
ann_colours <- list(
  Group = c("Normal" = "#4DAFFF", "Tumour" = "#FF6B6B")
)

pdf("outputs/figures/04_heatmap_top50_DEGs.pdf", width = 12, height = 14)
pheatmap(
  hm_scaled,
  annotation_col  = annotation_col,
  annotation_colors = ann_colours,
  show_colnames   = TRUE,
  show_rownames   = TRUE,
  fontsize_row    = 7,
  fontsize_col    = 8,
  cluster_cols    = TRUE,
  cluster_rows    = TRUE,
  color           = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  main            = "Top 50 Differentially Expressed Genes\n(z-scored expression)"
)
dev.off()
cat("Heatmap saved.\n")

# ── 8. Save results ──────────────────────────────────────────
write.csv(all_results,  "outputs/tables/04_all_DEGs.csv",         row.names = TRUE)
write.csv(sig_results,  "outputs/tables/04_significant_DEGs.csv", row.names = TRUE)

# Save Entrez IDs for enrichment
sig_entrez <- gene_ann[rownames(sig_results), "EntrezID"]
sig_entrez <- sig_entrez[!is.na(sig_entrez)]

bg_entrez  <- gene_ann$EntrezID
bg_entrez  <- bg_entrez[!is.na(bg_entrez)]

saveRDS(
  list(
    all_results = all_results,
    sig_results = sig_results,
    sig_entrez  = sig_entrez,
    bg_entrez   = bg_entrez,
    gene_expr   = gene_expr,
    gene_ann    = gene_ann,
    group       = group
  ),
  "outputs/rdata/04_deg_results.rds"
)

cat("\nDEG results saved.\n")
cat("\n--- Script 04 complete ---\n")
