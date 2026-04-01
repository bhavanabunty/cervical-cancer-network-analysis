# Script 09: Hub Gene Identification
# ============================================================

cat("\n========== Script 09: Hub Gene Identification ==========\n")

topo_out <- readRDS("outputs/rdata/08_topology.rds")
g        <- topo_out$g

# ── Centrality measures ──────────────────────────────────────
deg_cent    <- degree(g)
between_cent <- betweenness(g, normalized = TRUE)
close_cent  <- closeness(g, normalized = TRUE)
eigen_cent  <- eigen_centrality(g)$vector

centrality_df <- data.frame(
  StringID    = V(g)$name,
  GeneName    = V(g)$gene_name,
  Degree      = deg_cent,
  Betweenness = between_cent * vcount(g) * (vcount(g) - 1) / 2,  # raw betweenness
  Closeness   = close_cent,
  Eigenvector = eigen_cent
)

# ── Rank and define hub genes (top 15 by degree) ─────────────
centrality_df <- centrality_df[order(-centrality_df$Degree), ]
hub_threshold <- 15
hub_genes <- centrality_df[1:hub_threshold, ]

cat("\nTop 15 Hub Genes by Degree Centrality:\n")
print(hub_genes[, c("GeneName", "Degree", "Betweenness")])

# ── Visualise hub gene rankings ──────────────────────────────
p_hubs <- ggplot(hub_genes,
                 aes(x = reorder(GeneName, Degree), y = Degree, fill = Betweenness)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#AED6F1", high = "#1A5276",
                      name = "Betweenness\nCentrality") +
  labs(title    = "Top 15 Hub Genes by Degree Centrality",
       subtitle = "Fill colour shows betweenness centrality",
       x = "Gene", y = "Degree (number of interactions)") +
  theme_classic(base_size = 12)

ggsave("outputs/figures/09_hub_genes.pdf",
       plot = p_hubs, width = 8, height = 6)

write.csv(centrality_df, "outputs/tables/09_centrality_all_genes.csv", row.names = FALSE)
write.csv(hub_genes,     "outputs/tables/09_hub_genes.csv",             row.names = FALSE)
saveRDS(list(g = g, centrality_df = centrality_df, hub_genes = hub_genes),
        "outputs/rdata/09_hub_genes.rds")

cat("\n--- Script 09 complete ---\n")


# ============================================================
