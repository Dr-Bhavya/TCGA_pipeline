#!/usr/bin/env Rscript

# ============================================================================
# SCRIPT: rra_meta_analysis.R
# Description: Segregates differentially expressed gene lists by direction,
#              consolidates identifiers, executes Robust Rank Aggregation,
#              and strictly ranks outcomes by ascending RRA score values.
# Arguments:
#   - args[1]: Space-separated string of file paths to DEG results
#   - args[2]: Path to the master genomic annotation table
#   - args[3]: Path to the combined matrix file
# ============================================================================

library(RobustRankAggreg)

args <- commandArgs(trailingOnly = TRUE)
deg_files_str     <- args[1]
annotate_file     <- args[2]
combined_matrix_f <- args[3]

anno <- read.csv(annotate_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
file_vector <- unlist(strsplit(deg_files_str, " "))
alphabetical_order <- sort(gsub("_DEGs_sig.csv", "", basename(file_vector)))

ordered_log2fc_df <- read.csv(combined_matrix_f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

up_ids_list <- list()
dn_ids_list <- list()

# ============================================================================
# SECTION 1: DIRECTIONAL CANDIDATE SORTING
# ============================================================================
for (i in seq_along(file_vector)) {
    current_file <- file_vector[i]
    current_cancer <- gsub("_DEGs_sig.csv", "", basename(current_file))
    data <- read.csv(current_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    
    up_genes_only <- data[data$log2FoldChange >= 1, ]
    up_ids_list[[current_cancer]] <- up_genes_only[order(up_genes_only$log2FoldChange, decreasing = TRUE), ]$EnsemblID
    
    dn_genes_only <- data[data$log2FoldChange <= -1, ]
    dn_ids_list[[current_cancer]] <- dn_genes_only[order(dn_genes_only$log2FoldChange, decreasing = FALSE), ]$EnsemblID
}

# ============================================================================
# SECTION 2: IDENTIFIER CONSOLIDATION MATRIX EXPORT
# ============================================================================
max_len_up <- max(sapply(up_ids_list, length))
write.csv(as.data.frame(lapply(up_ids_list, function(x) c(x, rep("", max_len_up - length(x)))), check.names = FALSE)[, alphabetical_order], file = "All_Cancers_Upregulated_EnsemblIDs.csv", row.names = FALSE)

max_len_dn <- max(sapply(dn_ids_list, length))
write.csv(as.data.frame(lapply(dn_ids_list, function(x) c(x, rep("", max_len_dn - length(x)))), check.names = FALSE)[, alphabetical_order], file = "All_Cancers_Downregulated_EnsemblIDs.csv", row.names = FALSE)

# ============================================================================
# SECTION 3: ROBUST RANK AGGREGATION WRAPPER ENGINE
# ============================================================================
execute_rra <- function(ids_list, output_name, mode) {
    clean_lists <- lapply(ids_list, function(x) x[x != "" & !is.na(x)])
    rra_res <- aggregateRanks(glist = clean_lists, method = "RRA")
    
    rra_res$padj <- p.adjust(rra_res$Score, method = "BH")
    rra_res$EnsemblID <- rra_res$Name
    rra_res$GeneSymbol <- anno$GeneSymbol[match(rra_res$EnsemblID, anno$Geneid)]
    
    match_master <- match(rra_res$EnsemblID, ordered_log2fc_df$EnsemblID)
    rra_res$Total_DEG <- ordered_log2fc_df$Total_DEG[match_master]
    rra_res$Total_Up <- ordered_log2fc_df$Total_Up[match_master]
    rra_res$Total_Dn <- ordered_log2fc_df$Total_Dn[match_master]
    
    # Clean NA counts to prevent logical filter evaluation dropouts
    if (is.null(rra_res$Total_Up)) rra_res$Total_Up <- 0
    if (is.null(rra_res$Total_Dn)) rra_res$Total_Dn <- 0
    rra_res$Total_Up[is.na(rra_res$Total_Up)] <- 0
    rra_res$Total_Dn[is.na(rra_res$Total_Dn)] <- 0
    
    # ------------------------------------------------------------------------
    # CRITICAL FIX: EXPLICIT ORDER BY ASCENDING SCORE VALUE
    # Smallest RRA probability score = highest significance priority candidate
    # ------------------------------------------------------------------------
    rra_res <- rra_res[order(rra_res$Score, decreasing = FALSE), ]
    
    # Assign sequential Rank tracking metrics based on the sorted scores
    rra_res$Rank <- 1:nrow(rra_res)
    
    # Apply categorization rules (Green if Freq <= 3, Red/Blue if Freq > 3)
    rra_res$Category <- "Not_Significant"
    sig_idx <- !is.na(rra_res$padj) & rra_res$padj < 0.05
    
    if (mode == "UP") {
        low_freq_idx <- sig_idx & (rra_res$Total_Up <= 3)
        rra_res$Category[low_freq_idx] <- "Significant_Low_Freq"
        
        hub_idx <- sig_idx & (rra_res$Total_Up > 3)
        rra_res$Category[hub_idx] <- "Significant_Hub"
    } else {
        low_freq_idx <- sig_idx & (rra_res$Total_Dn <= 3)
        rra_res$Category[low_freq_idx] <- "Significant_Low_Freq"
        
        hub_idx <- sig_idx & (rra_res$Total_Dn > 3)
        rra_res$Category[hub_idx] <- "Significant_Hub"
    }
    
    # Generate log p-value column safely
    rra_res$log_p <- -log10(rra_res$padj)
    rra_res$log_p[is.infinite(rra_res$log_p)] <- max(rra_res$log_p[!is.infinite(rra_res$log_p)], na.rm = TRUE) + 2
    
    # Structure columns and export final file
    rra_res <- rra_res[, c("EnsemblID", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn", 
                           "Score", "padj", "Rank", "Category", "log_p")]
    colnames(rra_res) <- c("EnsemblID", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn", 
                           "RRA_Score_pvalue", "padj", "Rank", "Category", "log_p")
    
    write.csv(rra_res, file = output_name, row.names = FALSE)
}

execute_rra(up_ids_list, "RRA_Consensus_Upregulated.csv", "UP")
execute_rra(dn_ids_list, "RRA_Consensus_Downregulated.csv", "DN")