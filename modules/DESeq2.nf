#!/usr/bin/env nextflow

/* 
 * ============================================================================
 * MODULE: DESEQ2
 * Description: Conducts differential expression analysis using DESeq2 and 
 *              performs downstream log2 fold-change shrinkage profiling.
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
        /* 
         * Execution Context: Embedded R Script 
         */
        """
        #!/usr/bin/env Rscript
        library(DESeq2)

        # Parse count matrix and corresponding sample tracking data
        counts <- read.csv("${RNA_count}", row.names = 1, header = TRUE, check.names = FALSE)
        metadata <- read.csv("${metadata_path}", row.names = 1, header = TRUE, check.names = FALSE)

        # Align sample metadata observations with matrix target names
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

        # Verify integrity of matrix layout match configurations
        all(rownames(colData) %in% colnames(counts))
        all(rownames(colData) == colnames(counts))

        # Instantiate core DESeq2 structure and calculate properties
        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData = colData,
            design = ~ condition
        )
        dds\$condition <- factor(dds\$condition, levels = c("SolidTissueNormal", "PrimaryTumor"))
        dds <- DESeq(dds)
        res <- results(dds)
        
        # Apply empirical shrinkage for accurate low-count filtering
        resLFC <- lfcShrink(dds, coef = "condition_PrimaryTumor_vs_SolidTissueNormal", type = "apeglm")
        resLFC_df <- as.data.frame(resLFC)
        
        # Map genomic annotation data with Ensembl IDs
        anno <- read.csv("${annotate}", header = TRUE, check.names = FALSE)
        resLFC_df\$EnsemblID <- rownames(resLFC_df)

        matched_idx <- match(resLFC_df\$EnsemblID, anno\$Geneid)

        resLFC_df\$GeneSymbol <- anno\$GeneSymbol[matched_idx]
        resLFC_df\$Class <- anno\$Class[matched_idx]
        resLFC_df\$Chromosome <- anno\$Chromosome[matched_idx]

        # Restructure and clean columns
        base_cols <- c("EnsemblID", "GeneSymbol", "Class", "Chromosome")
        stat_cols <- c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
        resLFC_df <- resLFC_df[, c(base_cols, stat_cols)]

        # Isolate and filter significant Differentially Expressed Genes (DEGs)
        resLFC_filtered <- resLFC_df[!is.na(resLFC_df\$padj), ]
        resLFC_filtered <- resLFC_filtered[resLFC_filtered\$padj < 0.05 & abs(resLFC_filtered\$log2FoldChange) >= 1, ]
        write.csv(resLFC_filtered, file = "${cancer}_DEGs_sig.csv", row.names = FALSE)

        # Segment and sort upregulated datasets
        resLFC_up <- resLFC_filtered[resLFC_filtered\$log2FoldChange >= 1, ]
        resLFC_up_sorted <- resLFC_up[order(resLFC_up\$log2FoldChange, decreasing = TRUE), ]
        write.csv(resLFC_up_sorted, file = "${cancer}_DEGs_Up.csv", row.names = FALSE)

        # Segment and sort downregulated datasets
        resLFC_dn <- resLFC_filtered[resLFC_filtered\$log2FoldChange < -1, ]
        resLFC_dn_sorted <- resLFC_dn[order(resLFC_dn\$log2FoldChange, decreasing = FALSE), ]
        write.csv(resLFC_dn_sorted, file = "${cancer}_DEGs_Dn.csv", row.names = FALSE)
        """
}
