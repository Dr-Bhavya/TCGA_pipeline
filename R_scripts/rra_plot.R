#!/usr/bin/env Rscript

# ============================================================================
# SCRIPT: generate_plots.R
# Description: Reads enriched RRA data structures and renders prioritization plots.
#              Forces all categories to display in the legend keys.
# ============================================================================

library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
rra_up_csv   <- args[1]
rra_down_csv <- args[2]

generate_scatterplot <- function(csv_file, mode, pdf_name, tiff_name) {
    rra_res <- read.csv(csv_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    
    # Convert category column into a strict factor matching the preset palette keys
    rra_res$Category <- factor(rra_res$Category, levels = c("Not_Significant", "Significant_Low_Freq", "Significant_Hub"))
    
    if (mode == "UP") {
        color_palette <- c("Not_Significant" = "grey", "Significant_Low_Freq" = "green4", "Significant_Hub" = "red")
        legend_labels <- c("adj p-value >= 0.05", "Freq <= 3, adj p-value < 0.05", "Freq > 3, adj p-value < 0.05")
        plot_title <- "Target Prioritization: Upregulated RRA Hub Genes"
    } else {
        color_palette <- c("Not_Significant" = "grey", "Significant_Low_Freq" = "green4", "Significant_Hub" = "blue")
        legend_labels <- c("adj p-value >= 0.05", "Freq <= 3, adj p-value < 0.05", "Freq > 3, adj p-value < 0.05")
        plot_title <- "Target Prioritization: Downregulated RRA Hub Genes"
    }

    g_plot <- ggplot(rra_res, aes(x = Rank, y = log_p)) +
        # CRITICAL FIX: show.legend = TRUE forces ggplot to build legend items for all levels 
        geom_point(shape = 16, size = 2, aes(color = Category), show.legend = TRUE) + 
        scale_color_manual(values = color_palette, labels = legend_labels, drop = FALSE) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        theme_bw() +
        theme(
            panel.grid = element_blank(), panel.border = element_blank(),
            axis.line = element_line(colour = "black"), axis.text = element_text(size = 12),
            axis.title = element_text(size = 14), legend.position = c(0.70, 0.75),
            legend.title = element_blank(), legend.box.background = element_rect(colour = "black")
        ) +
        labs(title = plot_title, y = "-log10 (RRA adjusted p-value)", x = "RRA Rank") +
        geom_text_repel(data = rra_res[1:20, ], aes(label = GeneSymbol), size = 3.5, segment.color = "grey10", max.overlaps = Inf)

    ggsave(pdf_name, g_plot, width = 6, height = 6, device = "pdf")
    ggsave(tiff_name, g_plot, width = 6, height = 6, device = "tiff", dpi = 300)
}

generate_scatterplot(rra_up_csv, "UP", "RRA_Upregulated_Significance.pdf", "RRA_Upregulated_Significance.tiff")
generate_scatterplot(rra_down_csv, "DN", "RRA_Downregulated_Significance.pdf", "RRA_Downregulated_Significance.tiff")