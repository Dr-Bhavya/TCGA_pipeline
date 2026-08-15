#!/usr/bin/env nextflow

/* 
 * ============================================================================
 * MODULE INCLUDE STATEMENTS
 * ============================================================================
 */
include { GDCRNATOOLS } from './modules/GDCRNATools.nf'
include { DESEQ2 } from './modules/DESeq2.nf'
include { RRA } from './modules/RobustRankAggreg.nf'

/* 
 * ============================================================================
 * PIPELINE INPUT PARAMETERS
 * ============================================================================
 */
params {
    // Primary input (file of input manifest files, one per line)
    input_manifest: Path
    input_metadata: Path
    custom_merge_script: Path
    annotation: Path
}

/* 
 * ============================================================================
 * MAIN WORKFLOW BLOCK
 * ============================================================================
 */
workflow {
    main:
        // Initialize manifest channel and parse CSV contents
        read_manifest_ch = channel.fromPath(params.input_manifest)
        read_manifest_ch = channel.fromPath(params.input_manifest)
            .splitCsv(header: true)
            .map { row -> 
                tuple(
                    row.cancer_type, 
                    file(row.cancer_manifest_file_path), 
                    file(row.normal_manifest_file_path)
                ) 
            }

        // Initialize metadata channel and parse CSV contents
        read_metadata_ch = channel.fromPath(params.input_metadata)
            .splitCsv(header: true)
            .map { row -> 
                tuple(row.cancer_type, file(row.metadata_file_path)) 
            }

        // Join manifest data and clinical metadata on the common cancer type key
        combined_input_ch = read_manifest_ch.join(read_metadata_ch)

        // Process data downloading and parsing using GDCRNATools
        GDCRNATOOLS(combined_input_ch, file(params.custom_merge_script))
            
        // Construct DESeq2 input channel by linking counts with metadata
        deseq2_input_ch = GDCRNATOOLS.out.count.join(read_metadata_ch)
        DESEQ2(deseq2_input_ch, params.annotation)

        // Extract and aggregate all differentially expressed gene files from all cohorts
        all_deg_up_files_ch = DESEQ2.out.degs.map { cancer, file_path -> file_path }.collect()
        
        // Compute Robust Rank Aggregation cross-study comparison
        RRA(all_deg_up_files_ch, params.annotation)

    publish:
        // Expose internal process output channels to the workflow output context
        downloaded_files = GDCRNATOOLS.out.tcga_folder
        count_matrix = GDCRNATOOLS.out.count
        condition_files = DESEQ2.out.coldata
        degs_results = DESEQ2.out.degs
        degs_up = DESEQ2.out.degs_up
        degs_dn = DESEQ2.out.degs_dn
        degs_matrix = RRA.out.log2fc_matrix
        ranked_up_mat = RRA.out.up_ids_matrix
        ranked_dn_mat = RRA.out.dn_ids_matrix
        rra_cons_up = RRA.out.rra_up
        rra_cons_dn = RRA.out.rra_dn
        pdf_rra_up = RRA.out.pdf_up
        tiff_rra_up = RRA.out.tiff_up
        pdf_rra_dn = RRA.out.pdf_dn
        tiff_rra_dn = RRA.out.tiff_dn
}

/* 
 * ============================================================================
 * WORKFLOW OUTPUTS ROUTING BLOCK
 * ============================================================================
 */
output {
    downloaded_files {
        path 'STAR_files'
    }
    count_matrix {
        path 'STAR_count_matrices'
    }
    condition_files {
        path 'STAR_files'
    }
    degs_results {
        path 'Significant_DEGs'
    }
    degs_up {
        path 'Upregulated_DEGs'
    }
    degs_dn {
        path 'Downregulated_DEGs'
    }  
    degs_matrix {
        path 'Significant_DEGs'
    }  
    ranked_up_mat {
        path 'Upregulated_DEGs'
    }  
    ranked_dn_mat {
        path 'Downregulated_DEGs'
    }
    rra_cons_up {
        path 'Upregulated_DEGs'
    }  
    rra_cons_dn {
        path 'Downregulated_DEGs'
    }  
    pdf_rra_up { 
        path 'Upregulated_DEGs' 
    }
    tiff_rra_up { 
        path 'Upregulated_DEGs' 
    }   
    pdf_rra_dn { 
        path 'Downregulated_DEGs' 
    }   
    tiff_rra_dn { 
        path 'Downregulated_DEGs' 
    }   
}
