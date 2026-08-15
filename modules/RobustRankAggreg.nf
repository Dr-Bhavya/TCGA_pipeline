process RRA {
    container 'docker_tcga_env:latest'

    input:
        path deg_files
        path annotate
        path script_merge
        path script_rra
        path script_plot

    output:
        path "Differential_Combined_Log2FC.csv"        , emit: log2fc_matrix
        path "All_Cancers_Upregulated_EnsemblIDs.csv"  , emit: up_ids_matrix
        path "All_Cancers_Downregulated_EnsemblIDs.csv", emit: dn_ids_matrix
        path "RRA_Consensus_Upregulated.csv"           , emit: rra_up
        path "RRA_Consensus_Downregulated.csv"         , emit: rra_dn
        path "RRA_Upregulated_Significance.pdf"        , emit: pdf_up
        path "RRA_Upregulated_Significance.tiff"       , emit: tiff_up
        path "RRA_Downregulated_Significance.pdf"      , emit: pdf_dn
        path "RRA_Downregulated_Significance.tiff"     , emit: tiff_dn

    script: 
        def joined_files = deg_files.join(' ')
        
        """
        chmod +x ${script_merge} ${script_rra} ${script_plot}
        
        # Step 1: Outputs 'Differential_Combined_Log2FC.csv'
        ./${script_merge} "${joined_files}" "${annotate}"
        
        # Step 2: Passes the dynamic filename output of Step 1 into args[3]
        ./${script_rra} "${joined_files}" "${annotate}" "Differential_Combined_Log2FC.csv"
        
        # Step 3: Passes the generated consensus outputs dynamically to the plotting script
        ./${script_plot} "RRA_Consensus_Upregulated.csv" "RRA_Consensus_Downregulated.csv"
        """
}