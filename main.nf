#!/usr/bin/env nextflow

/* 
 * ============================================================================
 * PIPELINE ORCHESTRATION WORKFLOW
 * Description: Main pipeline script integrating TCGA data fetching, DESeq2 
 *              differential analysis, and a structured RRA meta-analysis.
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
    input_manifest: Path
    input_metadata: Path
    annotation: Path
    custom_merge_script: Path
    rra_script_merge: Path
    rra_script_run: Path
    rra_script_plot: Path
}

/* 
 * ============================================================================
 * MAIN WORKFLOW BLOCK
 * ============================================================================
 */
workflow {

    main:
        // Parse the sample input manifest CSV file into structural tuple data
        read_manifest_ch = channel.fromPath(params.input_manifest)
            .splitCsv(header: true)
            .map { row -> 
                tuple(row.cancer_type, file(row.cancer_manifest_file_path), file(row.normal_manifest_file_path)) 
            }

        // Parse matching sample metadata details
        read_metadata_ch = channel.fromPath(params.input_metadata)
            .splitCsv(header: true)
            .map { row -> 
                tuple(row.cancer_type, file(row.metadata_file_path)) 
            }

        // Combine manifest data and clinical metadata on the common cancer type key
        combined_input_ch = read_manifest_ch.join(read_metadata_ch)

        // Stage 1: Execute downloading and counts construction
        GDCRNATOOLS(combined_input_ch, file(params.custom_merge_script))
            
        // Stage 2: Bind generated expression matrices with sample metadata for DESeq2
        deseq2_input_ch = GDCRNATOOLS.out.count.join(read_metadata_ch)
        DESEQ2(deseq2_input_ch, params.annotation)

        // Stage 3: Isolate and collect all downstream targets into a global array
        all_deg_files_ch = DESEQ2.out.degs.map { cancer, file_path -> file_path }.collect()
        
        // Stage 4: Run Robust Rank Aggregation tracking cross-cohort consensus signals
        RRA(
            all_deg_files_ch, 
            params.annotation, 
            file(params.rra_script_merge), 
            file(params.rra_script_run), 
            file(params.rra_script_plot)
        )

    publish:
        // Expose nested pipeline output channels to the target deployment block
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
        //pdf_rra_up = RRA.out.pdf_up
        tiff_rra_up = RRA.out.tiff_up
        //pdf_rra_dn = RRA.out.pdf_dn
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
    tiff_rra_up { 
        path 'Upregulated_DEGs' 
    }    
    tiff_rra_dn { 
        path 'Downregulated_DEGs' 
    }
    /*pdf_rra_up {
         path "s3://tcga-pipeline-data-lake/results"
        }   
    pdf_rra_dn { 
        path "s3://tcga-pipeline-data-lake/results" 
    }   
    */
}
