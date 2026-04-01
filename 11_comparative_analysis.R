# Script 11: Comparative Analysis (Classical vs Network)
# ============================================================

cat("\n========== Script 11: Comparative Analysis ==========\n")

go_out   <- readRDS("outputs/rdata/05_go_results.rds")
kegg_out <- readRDS("outputs/rdata/06_kegg_results.rds")
mod_out  <- readRDS("outputs/rdata/10_modules.rds")

# Pathway sets
go_terms   <- as.data.frame(go_out$go_bp)$Description
kegg_paths <- as.data.frame(kegg_out)$Description
all_mod_df <- do.call(rbind, mod_out$module_enrichment)
net_terms  <- unique(all_mod_df$Description)

classical_set <- unique(c(go_terms, kegg_paths))
network_set   <- net_terms

# ── Jaccard Similarity ────────────────────────────────────────
intersection <- length(intersect(classical_set, network_set))
union_set    <- length(union(classical_set, network_set))
jaccard      <- intersection / union_set

cat("\n--- Comparative Metrics ---\n")
cat("Classical pathways:          ", length(classical_set), "\n")
cat("Network pathways:            ", length(network_set), "\n")
cat("Intersection:                ", intersection, "\n")
cat("Jaccard Similarity:          ", round(jaccard, 4), "\n")

# ── Novelty Rate ──────────────────────────────────────────────
novel_net      <- setdiff(network_set, classical_set)
novelty_rate   <- length(novel_net) / length(network_set)
cat("Novel network-only pathways: ", length(novel_net), "\n")
cat("Novelty Rate:                ", round(novelty_rate * 100, 1), "%\n")

# ── Biological Coherence Score ────────────────────────────────
# Known CC-relevant process keywords
cc_keywords <- c(
  "cell cycle", "dna replication", "dna repair", "apoptosis",
  "chromosome segregation", "mitotic", "proliferation",
  "p53", "tp53", "hpv", "carcinogenesis", "tumour",
  "immune", "signalling", "kinase", "keratinocyte"
)

score_coherence <- function(terms) {
  terms_lower <- tolower(terms)
  scores <- sapply(terms_lower, function(t) {
    if (any(sapply(cc_keywords[1:5], function(k) grepl(k, t)))) return(3)
    if (any(sapply(cc_keywords[6:10], function(k) grepl(k, t)))) return(2)
    if (any(sapply(cc_keywords[11:16], function(k) grepl(k, t)))) return(1)
    return(0)
  })
  round(sum(scores) / (3 * length(scores)), 3)
}

bc_classical <- score_coherence(classical_set)
bc_network   <- score_coherence(network_set)

cat("\nBiological Coherence Score:\n")
cat("  Classical: ", bc_classical, "\n")
cat("  Network:   ", bc_network, "\n")

# ── Summary table ─────────────────────────────────────────────
comparison_summary <- data.frame(
  Metric    = c("Total pathways", "Novel pathways",
                "Novelty Rate (%)", "Jaccard Similarity",
                "Biological Coherence Score"),
  Classical = c(length(classical_set), NA, NA, round(jaccard, 4), bc_classical),
  Network   = c(length(network_set), length(novel_net),
                round(novelty_rate * 100, 1), round(jaccard, 4), bc_network)
)

write.csv(comparison_summary, "outputs/tables/11_comparative_summary.csv",
          row.names = FALSE)

saveRDS(list(jaccard = jaccard, novelty_rate = novelty_rate,
             bc_classical = bc_classical, bc_network = bc_network,
             novel_net = novel_net),
        "outputs/rdata/11_comparative.rds")

cat("\n--- Script 11 complete ---\n")


# ============================================================
