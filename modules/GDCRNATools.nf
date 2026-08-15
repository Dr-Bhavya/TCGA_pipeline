#!/usr/bin/env nextflow

/* 
 * ============================================================================
 * MODULE: GDCRNATOOLS
 * Description: Downloads TCGA manifest components and assembles count matrices
 * ============================================================================
 */
process GDCRNATOOLS {
    container 'docker_tcga_env:latest'

    input:
        tuple val(cancer), path(manifest_cancer), path(manifest_normal), path(metadata_path)
        path merge_script
    
    output:
        tuple val(cancer), path("${cancer}_Data"), emit: tcga_folder
        tuple val(cancer), path("${cancer}_RNA_counts.csv"), emit: count

    script:
        /* 
         * Execution Context: Embedded R Script 
         */
        """
        #!/usr/bin/env Rscript

        source("${merge_script}")

        # Import demographic and clinical annotations
        metadata_file <- read.csv("${metadata_path}", row.names = 1, header = TRUE, check.names = FALSE)

        # Consolidate target paths into a single mapping vector
        manifest_vector <- c("${manifest_cancer}", "${manifest_normal}")

        # Applies gdcRNADownload to each manifest in the vector
        lapply(manifest_vector, function(x) {              
            gdcRNADownload(
                manifest = x,
                directory = "${cancer}_Data",
                method = "gdc-client"
            )
        })

        # Assemble individual sample files into a single unified matrix
        rnaCounts <- gdcRNAMerge(
            metadata = metadata_file, 
            path = "${cancer}_Data",
            organized = FALSE,
            data.type = "RNAseq"
        )
                        
        # Export structured matrix to persistent storage
        write.csv(rnaCounts, file = "${cancer}_RNA_counts.csv")
        """
}
