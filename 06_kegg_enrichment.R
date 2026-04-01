# Script 06: KEGG Pathway Enrichment
# ============================================================

cat("\n========== Script 06: KEGG Enrichment ==========\n")

kegg_res <- enrichKEGG(
  gene          = sig_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 5,
  maxGSSize     = 500
)

cat("KEGG pathways enriched:", nrow(as.data.frame(kegg_res)), "\n")

# Report top pathways
kegg_df <- as.data.frame(kegg_res)
if (nrow(kegg_df) > 0) {
  cat("\nTop 10 KEGG pathways:\n")
  top_kegg <- head(kegg_df[order(kegg_df$p.adjust), ], 10)
  print(top_kegg[, c("Description", "GeneRatio", "p.adjust")])
}

p_kegg <- dotplot(kegg_res, showCategory = 20, font.size = 10) +
  labs(title = "Top 20 Enriched KEGG Pathways",
       subtitle = "Cervical Cancer DEGs (adj.p < 0.05)") +
  theme_classic(base_size = 11)

ggsave("outputs/figures/06_kegg_dotplot.pdf",
       plot = p_kegg, width = 10, height = 9)

write.csv(kegg_df, "outputs/tables/06_kegg_results.csv", row.names = FALSE)
saveRDS(kegg_res, "outputs/rdata/06_kegg_results.rds")

cat("\n--- Script 06 complete ---\n")


# ============================================================
