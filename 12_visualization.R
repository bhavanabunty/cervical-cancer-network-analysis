# Script 12: Final Visualisation
# ============================================================

cat("\n========== Script 12: Final Visualisation ==========\n")
library(igraph)
library(ggplot2)

hub_out  <- readRDS("outputs/rdata/09_hub_genes.rds")
mod_out  <- readRDS("outputs/rdata/10_modules.rds")
comp_out <- readRDS("outputs/rdata/11_comparative.rds")
g        <- mod_out$g

# ── Network graph coloured by module ─────────────────────────
set.seed(42)

module_colours <- c(
  "1" = "#E74C3C", "2" = "#3498DB", "3" = "#2ECC71",
  "4" = "#F39C12", "5" = "#9B59B6", "6" = "#1ABC9C"
)
node_colours <- module_colours[as.character(V(g)$module)]
node_colours[is.na(node_colours)] <- "grey70"

hub_names    <- hub_out$hub_genes$GeneName[1:10]
node_sizes   <- ifelse(V(g)$gene_name %in% hub_names, 12, 5)
node_labels  <- ifelse(V(g)$gene_name %in% hub_names, V(g)$gene_name, NA)

pdf("outputs/figures/12_ppi_network_modules.pdf", width = 12, height = 10)
plot(
  g,
  vertex.color  = node_colours,
  vertex.size   = node_sizes,
  vertex.label  = node_labels,
  vertex.label.cex  = 0.8,
  vertex.label.color = "black",
  vertex.frame.color = "white",
  edge.color    = "grey80",
  edge.width    = 0.5,
  layout        = layout_with_fr(g),
  main          = "PPI Network: Functional Modules in Cervical Cancer",
  sub           = "Node size proportional to hub status; colours denote modules"
)
legend("bottomleft",
       legend = paste("Module", 1:3),
       fill   = module_colours[1:3],
       title  = "Module", bty = "n", cex = 0.9)
dev.off()
cat("PPI network plot saved.\n")

# ── Summary bar chart of comparative metrics ──────────────────
comp_df <- data.frame(
  Method  = c("Classical", "Network"),
  Coherence = c(comp_out$bc_classical, comp_out$bc_network)
)

p_comp <- ggplot(comp_df, aes(x = Method, y = Coherence, fill = Method)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = Coherence), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Classical" = "#4DAFFF", "Network" = "#FF6B6B")) +
  ylim(0, 1) +
  labs(title    = "Biological Coherence Score: Classical vs Network Enrichment",
       subtitle = paste0("Jaccard similarity = ",
                         round(comp_out$jaccard, 3),
                         "  |  Novel network pathways: ",
                         round(comp_out$novelty_rate * 100, 1), "%"),
       x = "Enrichment Method",
       y = "Biological Coherence Score (0–1)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave("outputs/figures/12_comparative_coherence.pdf",
       plot = p_comp, width = 7, height = 5)
cat("Comparative chart saved.\n")

# ── Session info ──────────────────────────────────────────────
sink("session_info.txt")
cat("=== R Session Information ===\n")
cat("Generated:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()
cat("Session info saved to session_info.txt\n")

cat("\n========================================\n")
cat("ALL 12 SCRIPTS COMPLETE\n")
cat("Results saved in outputs/figures/ and outputs/tables/\n")
cat("========================================\n")
