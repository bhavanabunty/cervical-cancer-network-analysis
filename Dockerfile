FROM rocker/bioconductor:3.17

LABEL maintainer="Bhavana Chowdary Kothapalli <a24bhako@student.his.se>"
LABEL description="Cervical Cancer Network-Based Pathway Enrichment Analysis"
LABEL version="1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN Rscript -e "install.packages(c( \
    'ggplot2', \
    'ggrepel', \
    'pheatmap', \
    'igraph', \
    'rmarkdown', \
    'knitr' \
), repos='https://cloud.r-project.org')"

# Install Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'GEOquery', \
    'affy', \
    'limma', \
    'preprocessCore', \
    'hgu133plus2.db', \
    'org.Hs.eg.db', \
    'AnnotationDbi', \
    'clusterProfiler', \
    'GOSemSim', \
    'STRINGdb', \
    'EnhancedVolcano', \
    'Biobase' \
), update=FALSE)"

# Set working directory
WORKDIR /workspace

# Default command
CMD ["Rscript", "--version"]
