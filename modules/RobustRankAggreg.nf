process RRA {
    container 'docker_tcga_env:latest'

    input:
    path(deg_files)
    path(annotate)

    output:
    path("Differential_Combined_Log2FC.csv")           , emit: log2fc_matrix
    path("All_Cancers_Upregulated_EnsemblIDs.csv")     , emit: up_ids_matrix
    path("All_Cancers_Downregulated_EnsemblIDs.csv")   , emit: dn_ids_matrix
    path("RRA_Consensus_Upregulated.csv")              , emit: rra_up
    path("RRA_Consensus_Downregulated.csv")            , emit: rra_dn
    path("RRA_Upregulated_Significance.pdf")           , emit: pdf_up
    path("RRA_Upregulated_Significance.tiff")          , emit: tiff_up
    path("RRA_Downregulated_Significance.pdf")          , emit: pdf_dn
    path("RRA_Downregulated_Significance.tiff")         , emit: tiff_dn


    script:
    """
    #!/usr/bin/env Rscript
    library(RobustRankAggreg)
    
    anno <- read.csv("${annotate}", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    
    master_matrix <- data.frame(EnsemblID = anno\$Geneid, stringsAsFactors = FALSE)
    master_matrix\$GeneSymbol <- anno\$GeneSymbol

    file_vector <- unlist(strsplit("${deg_files}", " "))
    
    up_ids_list <- list()
    dn_ids_list <- list()
    raw_cancer_names <- c()

    for (i in 1:length(file_vector)) {
        current_file   <- file_vector[i]
        current_cancer <- gsub("_DEGs_sig.csv", "", basename(current_file))
        raw_cancer_names <- c(raw_cancer_names, current_cancer)

        
        data <- read.csv(current_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
        
        # --- 1. ISOLATE AND SORT UPREGULATED IDs (Descending) ---
        up_genes_only   <- data[data\$log2FoldChange >= 1, ]
        up_genes_sorted <- up_genes_only[order(up_genes_only\$log2FoldChange, decreasing = TRUE), ]
        up_ids_list[[current_cancer]] <- up_genes_sorted\$EnsemblID
        
        # --- 2. ISOLATE AND SORT DOWNREGULATED IDs (Ascending) ---
        # Most negative values (highest fold change silencing) sit at the top
        dn_genes_only   <- data[data\$log2FoldChange <= -1, ]
        dn_genes_sorted <- dn_genes_only[order(dn_genes_only\$log2FoldChange, decreasing = FALSE), ]
        dn_ids_list[[current_cancer]] <- dn_genes_sorted\$EnsemblID
        
        # --- 3. MASTER LOG2FC MATRIX AGGREGATION ---
        sub_data <- data[, c("EnsemblID", "log2FoldChange")]
        colnames(sub_data) <- c("EnsemblID", current_cancer)
        
        match_idx <- match(master_matrix\$EnsemblID, data\$EnsemblID)
        master_matrix[[current_cancer]] <- data\$log2FoldChange[match_idx]
    }
        alphabetical_order <- sort(raw_cancer_names)

    numeric_grid <- master_matrix[, alphabetical_order, drop = FALSE]
    master_matrix\$Total_DEG <- rowSums(!is.na(numeric_grid))    
    master_matrix\$Total_Up  <- rowSums(!is.na(numeric_grid) & numeric_grid >= 1)    
    master_matrix\$Total_Dn  <- rowSums(!is.na(numeric_grid) & numeric_grid <= -1)
    
    meta_cols <- c("EnsemblID", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn")
    ordered_log2fc_df <- master_matrix[, c(meta_cols, alphabetical_order), drop = FALSE]
    write.csv(ordered_log2fc_df, file = "Differential_Combined_Log2FC.csv", row.names = FALSE)

    # =========================================================================
    # PHASE 2: CONSOLIDATING SORTED ID SHEETS (Unequal Row Length Padding)
    # =========================================================================
    print(">>> Phase 2: Building consolidated sorted Ensembl ID sheets...")
    
    # Consolidate Upregulated Sheets
    max_len_up <- max(sapply(up_ids_list, length))
    padded_list_up <- lapply(up_ids_list, function(x) c(x, rep("", max_len_up - length(x))))
    id_matrix_up <- as.data.frame(padded_list_up, check.names = FALSE)
    id_matrix_up <- id_matrix_up[, alphabetical_order, drop = FALSE]
    write.csv(id_matrix_up, file = "All_Cancers_Upregulated_EnsemblIDs.csv", row.names = FALSE)
    
    # Consolidate Downregulated Sheets
    max_len_dn <- max(sapply(dn_ids_list, length))
    padded_list_dn <- lapply(dn_ids_list, function(x) c(x, rep("", max_len_dn - length(x))))
    id_matrix_dn <- as.data.frame(padded_list_dn, check.names = FALSE)
    id_matrix_dn <- id_matrix_dn[, alphabetical_order, drop = FALSE]
    write.csv(id_matrix_dn, file = "All_Cancers_Downregulated_EnsemblIDs.csv", row.names = FALSE)
    
    print(">>> Matrix consolidation pipeline complete!")

library(ggplot2)
library(ggrepel)
    run_list_rra <- function(compiled_ids_list, output_name, mode) {
    # RRA list engine requires unpadded lists (strip out empty string padding elements)
        clean_lists <- list()
        for (cancer in names(compiled_ids_list)) {
            vector_vals <- compiled_ids_list[[cancer]]
            clean_vector <- vector_vals[vector_vals != "" & !is.na(vector_vals)]
            clean_lists[[cancer]] <- clean_vector
        }

        rra_res <- aggregateRanks(glist = clean_lists, method = "RRA")

        # Add a Benjamini-Hochberg False Discovery Rate correction column
        rra_res\$padj <- p.adjust(rra_res\$Score, method = "BH")
                
        # Match Ensembl IDs back to our master annotation map to preserve Gene Symbols
        rra_res\$GeneSymbol <- anno\$GeneSymbol[match(rra_res\$Name, anno\$Geneid)]

        match_master <- match(rra_res\$Name, ordered_log2fc_df\$EnsemblID)
        rra_res\$Total_DEG <- ordered_log2fc_df\$Total_DEG[match_master]
        rra_res\$Total_Up  <- ordered_log2fc_df\$Total_Up[match_master]
        rra_res\$Total_Dn  <- ordered_log2fc_df\$Total_Dn[match_master]
        
        # Restructure columns and export final priority candidates sheet
        rra_res <- rra_res[, c("Name", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn", "Score", "padj"), drop = FALSE]
        colnames(rra_res) <- c("EnsemblID", "GeneSymbol", "Total_DEG", "Total_Up", "Total_Dn", "RRA_Score_pvalue", "padj")
        
        rra_res\$Rank <- 1:nrow(rra_res)

        # --- APPLIED LOOKUP CATEGORIZATION BASED EXACTLY ON YOUR SPECIFIED DATA ---
        rra_res\$Category <- "Not_Significant" # Grey default base
        
        sig_idx <- !is.na(rra_res\$padj) & rra_res\$padj < 0.05
        rra_res\$Category[sig_idx] <- "Significant_Low_Freq" # Green4 default
        
        if (mode == "UP") {
            hub_idx <- sig_idx & !is.na(rra_res\$Total_Up) & rra_res\$Total_Up >= 3
            rra_res\$Category[hub_idx] <- "Significant_Hub" # Red
            color_palette <- c("Not_Significant" = "grey", "Significant_High_Freq" = "green4", "Significant_Hub" = "red")
            legend_labels <- c("adj p-value >= 0.05", "Freq < 3, adj p-value < 0.05", "Freq <= 3, adj p-value < 0.05")
            plot_title <- "Target Prioritization: Upregulated RRA Hub Genes"
            pdf_file <- "RRA_Upregulated_Significance.pdf"
            tiff_file <- "RRA_Upregulated_Significance.tiff"
        } else {
            hub_idx <- sig_idx & !is.na(rra_res\$Total_Dn) & rra_res\$Total_Dn >= 3
            rra_res\$Category[hub_idx] <- "Significant_Hub" # Blue
            color_palette <- c("Not_Significant" = "grey", "Significant_Low_Freq" = "green4", "Significant_Hub" = "blue")
            legend_labels <- c("adj p-value >= 0.05", "Freq < 3, adj p-value < 0.05", "Freq >= 3, adj p-value < 0.05")
            plot_title <- "Target Prioritization: Downregulated RRA Hub Genes"
            pdf_file <- "RRA_Downregulated_Significance.pdf"
            tiff_file <- "RRA_Downregulated_Significance.tiff"
        }

        # Handle extreme significance values (removes infinity limits)
        rra_res\$log_p <- -log10(rra_res\$padj)
        rra_res\$log_p[is.infinite(rra_res\$log_p)] <- max(rra_res\$log_p[!is.infinite(rra_res\$log_p)], na.rm=TRUE) + 2

        # Construct Plot Core Layout Configurations
        g_plot <- ggplot(rra_res, aes(x = Rank, y = log_p)) +
            geom_point(shape = 16, size = 2, aes(color = Category)) +
            scale_color_manual(values = color_palette, labels = legend_labels) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                legend.position = c(0.70, 0.75),
                legend.title = element_blank(),
                legend.box.background = element_rect(colour = "black")
            ) +
            labs(
                title = plot_title,
                y = "-log10 (RRA adjusted p-value)",
                x = "RRA Rank"
            ) +
            geom_text_repel(
                data = rra_res[1:20, ],
                aes(label = GeneSymbol),
                size = 3.5,
                segment.color = "grey10",
                max.overlaps = Inf
            )
        write.csv(rra_res, file = output_name, row.names = FALSE)
        ggsave(pdf_file, g_plot, width = 7, height = 7, device = "pdf")
        ggsave(tiff_file, g_plot, width = 7, height = 7, device = "tiff", dpi = 300)
    }

    run_list_rra(up_ids_list, "RRA_Consensus_Upregulated.csv", "UP")
    run_list_rra(dn_ids_list, "RRA_Consensus_Downregulated.csv", "DN")
    """
}