# Script 10: Functional Module Detection
# ============================================================

cat("\n========== Script 10: Module Detection ==========\n")
library(clusterProfiler)
library(org.Hs.eg.db)

hub_out <- readRDS("outputs/rdata/09_hub_genes.rds")
g       <- hub_out$g

# ── Louvain community detection ──────────────────────────────
set.seed(123)
communities <- cluster_louvain(g)
V(g)$module <- membership(communities)

n_modules <- length(unique(membership(communities)))
cat("Modules detected:", n_modules, "\n")
cat("Module sizes:\n")
print(table(membership(communities)))
cat("Modularity:", round(modularity(communities), 4), "\n")

# ── GO enrichment per module ─────────────────────────────────
module_enrichment <- list()
deg_out <- readRDS("outputs/rdata/04_deg_results.rds")
gene_ann <- deg_out$gene_ann

for (mod in sort(unique(V(g)$module))) {
  mod_genes <- V(g)$gene_name[V(g)$module == mod]
  mod_genes <- mod_genes[!is.na(mod_genes)]

  mod_entrez <- gene_ann$EntrezID[gene_ann$Symbol %in% mod_genes]
  mod_entrez <- mod_entrez[!is.na(mod_entrez)]

  if (length(mod_entrez) < 3) next

  go_mod <- tryCatch(
    enrichGO(gene = mod_entrez, OrgDb = org.Hs.eg.db,
             ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL
  )

  if (!is.null(go_mod) && nrow(as.data.frame(go_mod)) > 0) {
    df <- as.data.frame(go_mod)
    df$Module <- mod
    module_enrichment[[as.character(mod)]] <- df
    cat("Module", mod, ":", length(mod_genes), "genes,",
        nrow(df), "enriched GO terms\n")
  }
}

all_module_enrich <- do.call(rbind, module_enrichment)

write.csv(all_module_enrich, "outputs/tables/10_module_enrichment.csv",
          row.names = FALSE)

module_membership <- data.frame(
  StringID = V(g)$name,
  GeneName = V(g)$gene_name,
  Module   = V(g)$module
)
write.csv(module_membership, "outputs/tables/10_module_membership.csv",
          row.names = FALSE)

saveRDS(list(g = g, communities = communities,
             module_enrichment = module_enrichment),
        "outputs/rdata/10_modules.rds")

cat("\n--- Script 10 complete ---\n")


# ============================================================
