# ============================================================================
# BASE IMAGE DEFINITION
# Uses rocker/r2u tailored for ultra-fast Ubuntu binary installation of CRAN packages
# ============================================================================
FROM rocker/r2u:26.04@sha256:5dce84d9d79eaa937c1478455915c09bf27c159a10b90f2a5968cfd604f7b8d3

# ============================================================================
# SYSTEM DEPENDENCIES
# Installs core build utilities and documentation parsers with explicit pinning
# ============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake=4.2* \
    make=4.4* \
    pandoc=3.7* \
    && rm -rf /var/lib/apt/lists/*

# ============================================================================
# ENVIRONMENT & REPOSITORY CONFIGURATION
# Pins the CRAN repository snapshot mirror to maintain reproducibility
# ============================================================================
ENV CRAN_REPO="https://packagemanager.posit.co/cran/__linux__/resolute/2026-08-10"

RUN echo 'options(repos = c(CRAN = Sys.getenv("CRAN_REPO")))' >> /usr/lib/R/etc/Rprofile.site

# ============================================================================
# BIOINFORMATICS & R PACKAGES
# Stages Bioconductor environment and required data science libraries
# ============================================================================
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(version = '3.23', ask = FALSE); \
    options(warn = 2); \
    BiocManager::install(c('GDCRNATools', 'DESeq2', 'apeglm'), ask = FALSE); \
    install.packages(c('RobustRankAggreg', 'ggplot2', 'ggrepel', 'data.table'))"

# ============================================================================
# RUNTIME ENTRYPOINT
# Defaults container lifecycle to an interactive or scripted non-saving R session
# ============================================================================
CMD ["R", "--no-save", "--no-restore", "--no-environ"]
