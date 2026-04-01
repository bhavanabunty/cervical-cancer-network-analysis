# Script 05: GO Enrichment Analysis
# Project: Network-Based Pathway Enrichment in Cervical Cancer
# Author:  Bhavana Chowdary Kothapalli
# ============================================================

set.seed(123)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(ggplot2)

cat("Loading DEG results...\n")
deg_out    <- readRDS("outputs/rdata/04_deg_results.rds")
sig_entrez <- deg_out$sig_entrez
bg_entrez  <- deg_out$bg_entrez

# ── GO Biological Process ────────────────────────────────────
cat("Running GO Biological Process enrichment...\n")
go_bp <- enrichGO(
  gene          = sig_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ── GO Molecular Function ─────────────────────────────────────
cat("Running GO Molecular Function enrichment...\n")
go_mf <- enrichGO(
  gene          = sig_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ── GO Cellular Component ─────────────────────────────────────
cat("Running GO Cellular Component enrichment...\n")
go_cc <- enrichGO(
  gene          = sig_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("\nGO-BP enriched terms:", nrow(as.data.frame(go_bp)), "\n")
cat("GO-MF enriched terms:", nrow(as.data.frame(go_mf)), "\n")
cat("GO-CC enriched terms:", nrow(as.data.frame(go_cc)), "\n")

# ── Semantic similarity & redundancy reduction ─────────────────
cat("\nComputing semantic similarity for GO-BP (may take a few minutes)...\n")
go_bp_simplified <- simplify(
  go_bp,
  cutoff        = 0.7,
  by            = "p.adjust",
  select_fun    = min
)
cat("GO-BP terms after simplification:", nrow(as.data.frame(go_bp_simplified)), "\n")

# ── Dot plot ──────────────────────────────────────────────────
p_go_bp <- dotplot(go_bp, showCategory = 20, font.size = 10) +
  labs(title = "Top 20 Enriched GO Biological Processes",
       subtitle = "Cervical Cancer DEGs (adj.p < 0.05)") +
  theme_classic(base_size = 11)

ggsave("outputs/figures/05_go_bp_dotplot.pdf",
       plot = p_go_bp, width = 10, height = 9)

# ── Save ──────────────────────────────────────────────────────
write.csv(as.data.frame(go_bp), "outputs/tables/05_go_bp_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "outputs/tables/05_go_mf_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "outputs/tables/05_go_cc_results.csv", row.names = FALSE)
saveRDS(list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc,
             go_bp_simplified = go_bp_simplified),
        "outputs/rdata/05_go_results.rds")

cat("\n--- Script 05 complete ---\n")


# ============================================================
