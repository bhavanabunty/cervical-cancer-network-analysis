# Script 08: Network Topological Analysis
# ============================================================

cat("\n========== Script 08: Network Topology ==========\n")
library(igraph)

net_out      <- readRDS("outputs/rdata/07_string_network.rds")
interactions <- net_out$interactions
mapped       <- net_out$mapped

# Build igraph object
edges <- data.frame(
  from  = interactions$from,
  to    = interactions$to,
  score = interactions$combined_score
)

g <- graph_from_data_frame(edges, directed = FALSE,
                            vertices = mapped$STRING_id)

# Add gene name as vertex attribute
gene_lookup <- setNames(mapped$gene, mapped$STRING_id)
V(g)$gene_name <- gene_lookup[V(g)$name]

cat("Network: ", vcount(g), "nodes,", ecount(g), "edges\n")

# ── Global topology metrics ──────────────────────────────────
metrics <- list(
  n_nodes            = vcount(g),
  n_edges            = ecount(g),
  mean_degree        = mean(degree(g)),
  median_degree      = median(degree(g)),
  diameter           = diameter(g, unconnected = TRUE),
  avg_path_length    = mean_distance(g, unconnected = TRUE),
  density            = graph.density(g),
  clustering_coeff   = transitivity(g, type = "global")
)

cat("\n--- Network Global Metrics ---\n")
for (nm in names(metrics)) {
  cat(sprintf("  %-22s %s\n", nm, round(metrics[[nm]], 4)))
}

# ── Degree distribution & power law fit ─────────────────────
deg     <- degree(g)
deg_df  <- as.data.frame(table(deg), stringsAsFactors = FALSE)
deg_df$deg <- as.numeric(deg_df$deg)
deg_df  <- deg_df[deg_df$deg > 0, ]

# Log-log linear fit for power law assessment
fit_lm <- lm(log(Freq) ~ log(deg), data = deg_df)
gamma  <- -coef(fit_lm)[2]   # exponent
r2     <- summary(fit_lm)$r.squared

cat(sprintf("\nPower law: gamma = %.3f, R² = %.3f\n", gamma, r2))

p_degree <- ggplot(deg_df, aes(x = deg, y = Freq)) +
  geom_point(colour = "#2E5090", size = 2) +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(title = "Degree Distribution (log-log scale)",
       subtitle = sprintf("Power law fit: \u03B3 = %.2f, R\u00B2 = %.2f", gamma, r2),
       x = "Degree (k)", y = "Frequency P(k)") +
  theme_classic()

ggsave("outputs/figures/08_degree_distribution.pdf",
       plot = p_degree, width = 7, height = 5)

# Save
write.csv(data.frame(Metric = names(metrics), Value = unlist(metrics)),
          "outputs/tables/08_network_metrics.csv", row.names = FALSE)
saveRDS(list(g = g, metrics = metrics, gamma = gamma, r2 = r2),
        "outputs/rdata/08_topology.rds")

cat("\n--- Script 08 complete ---\n")


# ============================================================
