#!/usr/bin/env Rscript

# ============================================================================
# MATRIX AGGREGATION ENGINE
# Description: Compiles independent cohort log2 fold-change vectors into a 
#              unified multi-study matrix for downstream prioritization.
# Parameters:
#   - args[1]: Space-separated string of path targets to signature files
#   - args[2]: Centralized CSV gene translation annotation matrix
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
deg_files_str <- args[1]
annotate_file <- args[2]

anno <- read.csv(annotate_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
master_matrix <- data.frame(EnsemblID = anno$Geneid, stringsAsFactors = FALSE)
master_matrix$GeneSymbol <- anno$GeneSymbol

file_vector <- unlist(strsplit(deg_files_str, " "))
raw_cancer_names <- c()

# ============================================================================
# SECTION 1: COHORT EXTRPOLATION LOOP
# Iterates through available files to parse and merge local expression trends
# ============================================================================
for (i in seq_along(file_vector)) {
    current_file <- file_vector[i]
    current_cancer <- gsub("_DEGs_sig.csv", "", basename(current_file))
    raw_cancer_names <- c(raw_cancer_names, current_cancer)

    data <- read.csv(current_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    match_idx <- match(master_matrix$EnsemblID, data$EnsemblID)
    master_matrix[[current_cancer]] <- data$log2FoldChange[match_idx]
}

# ============================================================================
# SECTION 2: CONSENSUS METRICS CALCULATIONS
# Tabulates cumulative observation frequencies to highlight global markers
# ============================================================================
alphabetical_order <- sort(raw_cancer_names)
numeric_grid <- master_matrix[, alphabetical_order, drop = FALSE]

master_matrix$Total_DEG <- rowSums(!is.na(numeric_grid))    
master_matrix$Total_Up <- rowSums(!is.na(numeric_grid) & numeric_grid >= 1)    
master_matrix$Total_Dn <- rowSums(!is.na(numeric_grid) & numeric_grid <= -1)

meta_cols <- c("EnsemblID", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn")
ordered_log2fc_df <- master_matrix[, c(meta_cols, alphabetical_order), drop = FALSE]
write.csv(ordered_log2fc_df, file = "Differential_Combined_Log2FC.csv", row.names = FALSE)