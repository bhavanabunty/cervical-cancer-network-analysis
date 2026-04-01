# ============================================================
# Script 01: Data Acquisition
# Project: Network-Based Pathway Enrichment in Cervical Cancer
# Author:  Bhavana Chowdary Kothapalli
# Date:    2025
# ============================================================
# PURPOSE:
#   Downloads the GSE63514 dataset from NCBI GEO and saves
#   the raw expression set for downstream processing.
# ============================================================

set.seed(123)

# ── 1. Load libraries ────────────────────────────────────────
library(GEOquery)   # For downloading from GEO
library(Biobase)    # ExpressionSet class

# ── 2. Create output directories ────────────────────────────
dir.create("outputs",         showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/rdata",   showWarnings = FALSE, recursive = TRUE)

# ── 3. Download GSE63514 from GEO ───────────────────────────
cat("Downloading GSE63514 from NCBI GEO...\n")

gse <- getGEO(
  GEO        = "GSE63514",
  GSEMatrix  = TRUE,
  AnnotGPL   = TRUE,
  destdir    = "outputs/rdata/"
)

# GEOquery may return a list; extract the first element
if (is.list(gse)) {
  gse <- gse[[1]]
}

cat("Dataset downloaded successfully.\n")
cat("Class:", class(gse), "\n")
cat("Dimensions:", nrow(gse), "probes x", ncol(gse), "samples\n")

# ── 4. Inspect phenotype data ────────────────────────────────
pdata <- pData(gse)
cat("\nSample phenotype columns:\n")
print(colnames(pdata))

cat("\nSample types:\n")
print(table(pdata$`source_name_ch1`))

# ── 5. Save raw ExpressionSet ────────────────────────────────
saveRDS(gse, file = "outputs/rdata/gse_raw.rds")
cat("\nRaw ExpressionSet saved to outputs/rdata/gse_raw.rds\n")

cat("\n--- Script 01 complete ---\n")
