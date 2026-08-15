#!/usr/bin/env nextflow

/* 
 * ============================================================================
 * MODULE: DESEQ2
 * Description: Conducts differential gene expression analysis across sample
 *              cohorts and exports direction-specific significance tables.
 * ============================================================================
 */
process DESEQ2 {
    container 'docker_tcga_env:latest'

    input:
        tuple val(cancer), path(RNA_count), path(metadata_path)
        path annotate
    
    output:
        tuple val(cancer), path("${cancer}_coldata.csv"), emit: coldata
        tuple val(cancer), path("${cancer}_DEGs_sig.csv"), emit: degs
        tuple val(cancer), path("${cancer}_DEGs_Up.csv"), emit: degs_up
        tuple val(cancer), path("${cancer}_DEGs_Dn.csv"), emit: degs_dn

    script:
        // Execution Context: Embedded R Script 
    
        """
        #!/usr/bin/env Rscript
        library(DESeq2)

        counts <- read.csv("${RNA_count}", row.names = 1, header = TRUE, check.names = FALSE)
        metadata <- read.csv("${metadata_path}", row.names = 1, header = TRUE, check.names = FALSE)

        sample_names <- colnames(counts)
        matched_metadata <- metadata[match(sample_names, metadata\$sample), ]

        colData <- data.frame(
            row.names = sample_names,
            condition = matched_metadata\$sample_type,
            stringsAsFactors = FALSE
        )
        
        write.csv(colData, file = "${cancer}_coldata.csv")

        colData <- colData["condition"]
        colData\$condition <- factor(colData\$condition)

        all(rownames(colData) %in% colnames(counts))
        all(rownames(colData) == colnames(counts))
 
        # ============================================================================
        # SECTION 1: DESEQ2 MODELING & LOG2FC SHRINKAGE
        # ============================================================================
         
        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData = colData,
            design = ~ condition
        )
        dds\$condition <- factor(dds\$condition, levels = c("SolidTissueNormal", "PrimaryTumor"))
        dds <- DESeq(dds)
        res <- results(dds)
        resLFC <- lfcShrink(dds, coef = "condition_PrimaryTumor_vs_SolidTissueNormal", type = "apeglm")
        resLFC_df <- as.data.frame(resLFC)
        
        
        # ============================================================================
        # SECTION 2: TRANSCRIPT TRANSLATION & FILTERING
        # ============================================================================
        
        anno <- read.csv("${annotate}", header = TRUE, check.names = FALSE)
        resLFC_df\$EnsemblID <- rownames(resLFC_df)

        matched_idx <- match(resLFC_df\$EnsemblID, anno\$Geneid)

        resLFC_df\$GeneSymbol <- anno\$GeneSymbol[matched_idx]
        resLFC_df\$Class <- anno\$Class[matched_idx]
        resLFC_df\$Chromosome <- anno\$Chromosome[matched_idx]

        base_cols <- c("EnsemblID", "GeneSymbol", "Class", "Chromosome")
        stat_cols <- c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
        resLFC_df <- resLFC_df[, c(base_cols, stat_cols)]

        resLFC_filtered <- resLFC_df[!is.na(resLFC_df\$padj), ]
        resLFC_filtered <- resLFC_filtered[resLFC_filtered\$padj < 0.05 & abs(resLFC_filtered\$log2FoldChange) >= 1, ]
        write.csv(resLFC_filtered, file = "${cancer}_DEGs_sig.csv", row.names = FALSE)

        # ============================================================================
        # SECTION 3: DIRECTIONAL SIGNIFICANCE SORTING
        # ============================================================================
        
        resLFC_up <- resLFC_filtered[resLFC_filtered\$log2FoldChange >= 1, ]
        resLFC_up_sorted <- resLFC_up[order(resLFC_up\$log2FoldChange, decreasing = TRUE), ]
        write.csv(resLFC_up_sorted, file = "${cancer}_DEGs_Up.csv", row.names = FALSE)

        resLFC_dn <- resLFC_filtered[resLFC_filtered\$log2FoldChange < -1, ]
        resLFC_dn_sorted <- resLFC_dn[order(resLFC_dn\$log2FoldChange, decreasing = FALSE), ]
        write.csv(resLFC_dn_sorted, file = "${cancer}_DEGs_Dn.csv", row.names = FALSE)
        """
}