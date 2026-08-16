# TCGA_pipeline End-to-End Pipeline Engine

This repository hosts an automated bioinformatics pipeline engine built with **Nextflow** and containerised using **Docker**. The entire workflow is managed through an automated continuous integration pipeline using **GitHub Actions**.

---

## 🚀 CI/CD Pipeline Workflow

The automation workflow (`.github/workflows/ci.yml`) triggers on any **push** or **pull request** targeting the `main` or `master` branches. It specifically watches for changes in the `docker/` directory, core Nextflow scripts (`main.nf`, `nextflow.config`), and the workflow file itself.

### Key Steps in the CI/CD Pipeline
* **Code Checkout**: Fetches the code repo using `actions/checkout@v5`.
* **Docker Setup**: Configures standard and cached environment layers using `docker/setup-buildx-action@v3`.
* **Container Registry Login**: Authenticates into the GitHub Container Registry (GHCR) using a GitHub Personal Access Token (`secrets.GH_PAT`).
* **Image Management**: Lowercases repository variables, builds the container from `docker/Dockerfile`, tags it with the current commit SHA and `latest`, and pushes it to GHCR.
* **AWS Authentication**: Requests a secure OIDC Token and assumes the `EC2-S3-Bioinformatics-Role` in `ap-south-1` (Mumbai) region.
* **Nextflow Engine Run**: Downloads the latest Nextflow runtime engine and executes a validation pass using the integrated dry-run simulation mode (`-profile test -stub-run`).
* **Automated Commit**: Captures runtime outputs and updates the `results/` path directly inside your repository.

---

## 📦 Container Environment (Reproducibility Stack)

The runtime context is fully isolated inside a Docker container configured in `docker/Dockerfile`. This layout is designed for maximum speed, consistency, and analytical reproducibility:

* **Base Layer**: Built on top of `rocker/r2u` (Ubuntu base) utilizing optimized binary execution pools to install CRAN components.
* **System Utilities**: Packages core development stacks explicitly pinned to system variants (`cmake=4.2*`, `make=4.4*`, `pandoc=3.7*`).
* **Pinned Mirror Snapshots**: Binds runtime dependency pulls to a specific Posit Package Manager mirror timestamp (`2026-08-10`) to guarantee deterministic cross-environment builds.
* **Pre-bundled Toolsets**: Pre-installs full data science frameworks and specialized packages:
  * **Bioconductor (v3.23)**: Orchestrates core packages `GDCRNATools`, `DESeq2`, and `apeglm`.
  * **CRAN Packages**: Integrates algorithmic dependencies `RobustRankAggreg`, `ggplot2`, `ggrepel`, and `data.table`.

---

## 🧬 Data Processing & Meta-Analysis Details

The pipeline integrates complex Nextflow automation processes, R data-wrangling, statistical meta-analysis, and plotting routines to parse, aggregate, and report metrics across biological profiles:

### 1. Automated Cohort Fetching & Processing (`GDCRNATOOLS` Module)
* **Automated Batch Isolation**: Downloads specified TCGA cancer and normal data files using the embedded GDC API client payload layout via `gdcRNADownload`.
* **Parallel Matrix Generation**: Consolidates separate sample profile structures across distinct cohorts concurrently, routing raw expression arrays out as structural CSV records (`${cancer}_RNA_counts.csv`).
* **GDC Matrix Consolidation**: Runs underlying subroutines to combine raw sample data, parse pre-miRNA counts, and isolate mature microRNA sequences by filtering precursor elements.

### 2. Differential Expression Modeling (`DESEQ2` Module)
* **Design Matrix Assembly**: Automatically builds experimental condition maps (`${cancer}_coldata.csv`) mapping cohorts explicitly into `SolidTissueNormal` vs `PrimaryTumor`.
* **Log2 Fold-Change Shrinkage**: Runs negative binomial distributions via `DESeq2` and uses `apeglm` empirical Bayes shrinkage to reduce noise in low-count genes.
* **Annotation & Target Stratification**: Appends structural genomic markers (`GeneSymbol`, `Class`, `Chromosome`), enforcing an absolute filtration threshold ($\text{adjusted } p\text{-value} < 0.05$, $|\log_2\text{FC}| \ge 1$). 
* **Directional Splits**: Outputs standalone, sorted statistical subsets for total signature hits (`_DEGs_sig.csv`), upregulated markers (`_DEGs_Up.csv`), and downregulated markers (`_DEGs_Dn.csv`).

### 3. Meta-Analysis & Data Integration Orchestration (`RRA` Module)
* **Cloud Storage Sync**: Automatically exports rendering assets directly to a centralized AWS S3 Bucket environment (`s3://tcga-pipeline-data-lake/results`) using Nextflow `publishDir` tracking directives.
* **Sequential Execution Sync**: Pairs your entire asset tracking stack together, processing scripts in order through a single execution runtime channel:
  1. **Multi-Study Aggregation**: Combines signature paths (`${joined_files}`) using Ensembl IDs into a multi-study matrix tracking global directional variations.
  2. **Robust Rank Aggregation**: Computes non-parametric rankings sorted by ascending probability score fields (smallest RRA score indicates highest conservation and priority) and tags hits (`Not_Significant`, `Significant_Low_Freq`, or `Significant_Hub`).
  3. **Visual Target Reporting**: Feeds RRA outputs dynamically into plotting engines to map high-resolution vector figures.

---

## 🗺️ Pipeline Orchestration & Workflow Topology

The main workflow script (`main.nf`) orchestrates data movement sequentially across parallel asynchronous execution channels:

```text
[input_manifest] ──┐
                   ├──> [Channel Join] ──> GDCRNATOOLS ──> DESEQ2 ──> .collect() ──> RRA
[input_metadata] ──┘                        (Stage 1)     (Stage 2)     (Stage 3)     (Stage 4)
```

* **Dynamic Joins**: Ingests input manifest and clinical metadata CSV tracks independently, validating and pairing records on their matching `cancer_type` key before runtime initialization.
* **Array Aggregation**: Monitors step states across independent asynchronously evaluated tasks, aggregating all variable length `DESEQ2.out.degs` files into a static flat array (`.collect()`) before passing them to the final `RRA` meta-analysis module.
* **Structural Result Routing**: Maps execution outcomes systematically into targeted downstream physical folders (e.g., `STAR_count_matrices/`, `Significant_DEGs/`, `Upregulated_DEGs/`, and `Downregulated_DEGs/`).

---

## 🛠️ Prerequisites & Setup

To replicate or run this pipeline locally, you will need to prepare the following components:

### Local Development Requirements
* **Nextflow**: Required to parse and orchestrate data tasks.
* **Docker**: Needed to replicate the local runtime image.
* **Java 11 or later**: Required to support Nextflow executions.
* **R (v4.0+)**: Needed if executing the parsing scripts outside of Docker.

```bash
# Quick install Nextflow locally
wget -qO- https://nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### GitHub Environment Secrets
Before running the actions pipeline successfully, ensure you have provisioned this repository variable secret:
* `GH_PAT`: A GitHub Personal Access Token with permissions to write packages (`packages:write`) to the GitHub Container Registry.

---

## 💻 How to Run the Pipeline

### Running with Docker Container Profiles
The environment is pre-configured to run with isolated Docker images hosted on the GitHub Container Registry. To trigger the pipeline using the built-in test profile (which pulls input manifests directly from your AWS S3 data lake), execute:

```bash
nextflow run main.nf -profile test
```

### Local Dry-Run Simulation (Stub Run)
To validate execution routes and verify workflow plumbing logic without downloading full cohort source arrays:
```bash
nextflow run main.nf -profile test -stub-run
```

### Manual Command-Line Script Execution
The individual modules can be run manually inside an R-configured environment:
```bash
# 1. Multi-Cohort Matrix Aggregation
Rscript R_scripts/merge_log2fc.R "cohort1_DEGs_sig.csv cohort2_DEGs_sig.csv" gencode_v36_gene_annotation_table.csv

# 2. Robust Rank Aggregation Meta-Analysis
Rscript R_scripts/run_rra.R "cohort1_DEGs_sig.csv cohort2_DEGs_sig.csv" gencode_v36_gene_annotation_table.csv Differential_Combined_Log2FC.csv

# 3. Visual Reporting Plot Generation
Rscript R_scripts/rra_plot.R RRA_Consensus_Upregulated.csv RRA_Consensus_Downregulated.csv
```

---

## 📂 Project Structure

* **`.github/workflows/ci.yml`**: GitHub Actions automated pipeline setup.
* **`main.nf`**: Main Nextflow pipeline coordination and output routing workflow file.
* **`modules/GDCRNATools.nf`**: Nextflow process encapsulating raw cohort pulling and parsing matrices.
* **`modules/DESeq2.nf`**: Nextflow process managing conditional differential models and adaptive fold shrinkage.
* **`modules/RobustRankAggreg.nf`**: Nextflow meta-analysis connector, multi-study matrix compiler, and S3 publisher node.
* **`nextflow.config`**: Container orchestration profile definitions and environment file paths.
* **`docker/Dockerfile`**: Optimized reproducibility recipe utilizing speed-enhanced r2u binary layers.
* **`gencode_v36_gene_annotation_table.csv`**: Central transcript-to-gene translation mapping dataset.
* **`R_scripts/GDC_Merge_Refactored.R`**: R data wrangling script for GDC matrix assembly and mature miRNA processing.
* **`R_scripts/merge_log2fc.R`**: Matrix compiler calculating consensus multi-study Log2FC profiles.
* **`R_scripts/run_rra.R`**: Robust Rank Aggregation computation, sorting engine, and hub-gene classifier.
* **`R_scripts/rra_plot.R`**: R graphics routine for RRA Hub marker scatterplot prioritization maps.
* **`results/`**: Output directory for generated analytical metrics and pipeline artifacts.
