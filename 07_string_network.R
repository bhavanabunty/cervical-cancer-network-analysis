# Script 07: STRING PPI Network Construction
# ============================================================

cat("\n========== Script 07: STRING Network ==========\n")
library(STRINGdb)

deg_out    <- readRDS("outputs/rdata/04_deg_results.rds")
sig_results <- deg_out$sig_results

# Top 100 DEGs by adjusted p-value
top100 <- head(rownames(sig_results[order(sig_results$adj.P.Val), ]), 100)
top100_df <- data.frame(
  gene   = top100,
  logFC  = sig_results[top100, "logFC"],
  adj.p  = sig_results[top100, "adj.P.Val"]
)

cat("Top 100 DEGs selected for network analysis.\n")

# Initialise STRINGdb (human, version 11.5, medium confidence >= 400)
string_db <- STRINGdb$new(
  version     = "11.5",
  species     = 9606,
  score_threshold = 400,
  input_directory = "outputs/rdata/"
)

# Map gene symbols to STRING IDs
mapped <- string_db$map(top100_df, "gene", removeUnmappedRows = TRUE)
cat("Genes mapped to STRING:", nrow(mapped), "/", nrow(top100_df), "\n")

# Retrieve interactions
interactions <- string_db$get_interactions(mapped$STRING_id)
cat("Interactions retrieved:", nrow(interactions), "\n")

# Save
write.csv(mapped,       "outputs/tables/07_string_mapped_genes.csv",  row.names = FALSE)
write.csv(interactions, "outputs/tables/07_string_interactions.csv",   row.names = FALSE)
saveRDS(list(mapped = mapped, interactions = interactions,
             string_db = string_db, top100_df = top100_df),
        "outputs/rdata/07_string_network.rds")

cat("\n--- Script 07 complete ---\n")


# ============================================================
