load_data_for_HTML_report_from_result_files <- function(output_paths_file_path, has_run_tf_binding_site_analysis) {
    cat("Importing data...\n")

    result_paths <- data.table::fread(output_paths_file_path)
    n_samples <- length(result_paths$sample_id)

    # initialize empty lists ------------------------------------
    # by sample
    number_of_markers_vector <- vector(length = n_samples)
    expression_data_list <- vector("list", n_samples)
    somatic_SNV_data_list <- vector("list", n_samples)
    SV_data_list <- vector("list", n_samples)
    CNA_data_list <- vector("list", n_samples)
    cis_activation_summary_list <- vector("list", n_samples)
    SV_genes_all_list <- vector("list", n_samples)
    SV_genehancer_all_list <- vector("list", n_samples)
    SV_chipseq_all_list <- vector("list", n_samples)
    CNA_genes_all_list <- vector("list", n_samples)
    CNA_genehancer_all_list <- vector("list", n_samples)
    CNA_chipseq_all_list <- vector("list", n_samples)
    somatic_SNV_genes_all_list <- vector("list", n_samples)
    somatic_SNV_genehancer_all_list <- vector("list", n_samples)
    somatic_SNV_chipseq_all_list <- vector("list", n_samples)
    if (has_run_tf_binding_site_analysis) {
        somatic_SNV_tf_binding_data_genes_cis_activated_only_list <- vector("list", n_samples)
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_list <- vector("list", n_samples)
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff_list <- vector("list", n_samples)
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_list <- vector("list", n_samples)
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_list <- vector("list", n_samples)
        cis_activation_summary_with_n_relevant_somatic_SNVs_list <- vector("list", n_samples)
    }
    genehancer_with_cis_activation_summary_list <- vector("list", n_samples)


    # by subgroup
    chipseq_by_genes_list <- vector("list")


    # import data ----------------------------------------------
    for (i in seq_len(n_samples)) {
        sample_id <- result_paths$sample_id[i]


        number_of_markers_vector[i] <- nrow(readRDS(result_paths$markers_Rds[i]))
        expression_data_list[[sample_id]] <- readRDS(result_paths$expression_Rds[i])
        somatic_SNV_data_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_Rds[i])
        SV_data_list[[sample_id]] <- readRDS(result_paths$SV_Rds[i])
        CNA_data_list[[sample_id]] <- readRDS(result_paths$CNA_Rds[i])
        cis_activation_summary_list[[sample_id]] <- readRDS(result_paths$cis_activation_summary_Rds[i])
        SV_genes_all_list[[sample_id]] <- readRDS(result_paths$SV_genes_all_Rds[i])
        SV_genehancer_all_list[[sample_id]] <- readRDS(result_paths$SV_genehancer_all_Rds[i])
        SV_chipseq_all_list[[sample_id]] <- readRDS(result_paths$SV_chipseq_all_Rds[i])
        CNA_genes_all_list[[sample_id]] <- readRDS(result_paths$CNA_genes_all_Rds[i])
        CNA_genehancer_all_list[[sample_id]] <- readRDS(result_paths$CNA_genehancer_all_Rds[i])
        CNA_chipseq_all_list[[sample_id]] <- readRDS(result_paths$CNA_chipseq_all_Rds[i])
        somatic_SNV_genes_all_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_genes_all_Rds[i])
        somatic_SNV_genehancer_all_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_genehancer_all_Rds[i])
        somatic_SNV_chipseq_all_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_chipseq_all_Rds[i])
        if (has_run_tf_binding_site_analysis) {
            somatic_SNV_tf_binding_data_genes_cis_activated_only_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_tf_binding_data_genes_cis_activated_only_Rds[i])
            somatic_SNV_tf_binding_data_genehancer_cis_activated_only_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_tf_binding_data_genehancer_cis_activated_only_Rds[i])
            somatic_SNV_tf_binding_data_gene_summary_max_score_diff_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_tf_binding_data_gene_summary_max_score_diff_Rds[i])
            somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_Rds[i])
            somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_list[[sample_id]] <- readRDS(result_paths$somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_Rds[i])
            cis_activation_summary_with_n_relevant_somatic_SNVs_list[[sample_id]] <- readRDS(result_paths$cis_activation_summary_with_n_relevant_somatic_SNVs_Rds[i])
        }
        genehancer_with_cis_activation_summary_list[[sample_id]] <- readRDS(result_paths$genehancer_with_cis_activation_summary_Rds[i])

        if (!(result_paths$subgroup[i] %in% names(chipseq_by_genes_list))) {
            chipseq_by_genes_list[[result_paths$subgroup[i]]] <- readRDS(result_paths$chipseq_by_gene_Rds[i])
        }
    }

    # convert data to "long" tables -----------------------------
    marker_data_table <- data.table::setDT(data.frame(sample_ID = result_paths$sample_id, n_markers = number_of_markers_vector))
    expression_data_table <- data.table::rbindlist(expression_data_list, idcol = "sample_ID", fill = T)
    somatic_SNV_data_table <- data.table::rbindlist(somatic_SNV_data_list, idcol = "sample_ID", fill = T)
    SV_data_table <- data.table::rbindlist(SV_data_list, idcol = "sample_ID", fill = T)
    CNA_data_table <- data.table::rbindlist(CNA_data_list, idcol = "sample_ID", fill = T)
    cis_activation_summary_table <- data.table::rbindlist(cis_activation_summary_list, idcol = "sample_ID", fill = T) # sample_id already exists

    SV_genes_all_table <- data.table::rbindlist(SV_genes_all_list, idcol = "sample_ID", fill = T)
    SV_genehancer_all_table <- data.table::rbindlist(SV_genehancer_all_list, idcol = "sample_ID", fill = T)
    SV_chipseq_all_table <- data.table::rbindlist(SV_chipseq_all_list, idcol = "sample_ID", fill = T)
    CNA_genes_all_table <- data.table::rbindlist(CNA_genes_all_list, idcol = "sample_ID", fill = T)
    CNA_genehancer_all_table <- data.table::rbindlist(CNA_genehancer_all_list, idcol = "sample_ID", fill = T)
    CNA_chipseq_all_table <- data.table::rbindlist(CNA_chipseq_all_list, idcol = "sample_ID", fill = T)
    somatic_SNV_genes_all_table <- data.table::rbindlist(somatic_SNV_genes_all_list, idcol = "sample_ID", fill = T)
    somatic_SNV_genehancer_all_table <- data.table::rbindlist(somatic_SNV_genehancer_all_list, idcol = "sample_ID", fill = T)
    somatic_SNV_chipseq_all_table <- data.table::rbindlist(somatic_SNV_chipseq_all_list, idcol = "sample_ID", fill = T)

    if (has_run_tf_binding_site_analysis) {
        somatic_SNV_tf_binding_data_genes_cis_activated_only_table <- data.table::rbindlist(somatic_SNV_tf_binding_data_genes_cis_activated_only_list, idcol = "sample_ID", fill = T)
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table <- data.table::rbindlist(somatic_SNV_tf_binding_data_genehancer_cis_activated_only_list, idcol = "sample_ID", fill = T)
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table <- data.table::rbindlist(somatic_SNV_tf_binding_data_gene_summary_max_score_diff_list, idcol = "sample_ID", fill = T)
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table <- data.table::rbindlist(somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_list, idcol = "sample_ID", fill = T)
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table <- data.table::rbindlist(somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_list, idcol = "sample_ID", fill = T)
        cis_activation_summary_with_n_relevant_somatic_SNVs_table <- data.table::rbindlist(cis_activation_summary_with_n_relevant_somatic_SNVs_list, idcol = "sample_ID", fill = T)
    }

    genehancer_with_cis_activation_summary_table <- data.table::rbindlist(genehancer_with_cis_activation_summary_list, idcol = "sample_ID", fill = T)

    chipseq_by_genes_table <- data.table::rbindlist(chipseq_by_genes_list, idcol = "subgroup", fill = T)


    # create new data list (e.g. for cache) -----------------------------
    if (has_run_tf_binding_site_analysis) {
        return(
            list(
                result_paths = result_paths,
                marker_data_table = marker_data_table,
                expression_data_table = expression_data_table,
                somatic_SNV_data_table = somatic_SNV_data_table,
                SV_data_table = SV_data_table,
                CNA_data_table = CNA_data_table,
                cis_activation_summary_table = cis_activation_summary_table,
                SV_genes_all_table = SV_genes_all_table,
                SV_genehancer_all_table = SV_genehancer_all_table,
                SV_chipseq_all_table = SV_chipseq_all_table,
                CNA_genes_all_table = CNA_genes_all_table,
                CNA_genehancer_all_table = CNA_genehancer_all_table,
                CNA_chipseq_all_table = CNA_chipseq_all_table,
                somatic_SNV_genes_all_table = somatic_SNV_genes_all_table,
                somatic_SNV_genehancer_all_table = somatic_SNV_genehancer_all_table,
                somatic_SNV_chipseq_all_table = somatic_SNV_chipseq_all_table,
                somatic_SNV_tf_binding_data_genes_cis_activated_only_table = somatic_SNV_tf_binding_data_genes_cis_activated_only_table,
                somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table = somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table,
                somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table = somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table,
                somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table = somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table,
                somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table = somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table,
                cis_activation_summary_with_n_relevant_somatic_SNVs_table = cis_activation_summary_with_n_relevant_somatic_SNVs_table,
                genehancer_with_cis_activation_summary_table = genehancer_with_cis_activation_summary_table,
                chipseq_by_genes_table = chipseq_by_genes_table
            )
        )
    } else {
        return(
            list(
                result_paths = result_paths,
                marker_data_table = marker_data_table,
                expression_data_table = expression_data_table,
                somatic_SNV_data_table = somatic_SNV_data_table,
                SV_data_table = SV_data_table,
                CNA_data_table = CNA_data_table,
                cis_activation_summary_table = cis_activation_summary_table,
                SV_genes_all_table = SV_genes_all_table,
                SV_genehancer_all_table = SV_genehancer_all_table,
                SV_chipseq_all_table = SV_chipseq_all_table,
                CNA_genes_all_table = CNA_genes_all_table,
                CNA_genehancer_all_table = CNA_genehancer_all_table,
                CNA_chipseq_all_table = CNA_chipseq_all_table,
                somatic_SNV_genes_all_table = somatic_SNV_genes_all_table,
                somatic_SNV_genehancer_all_table = somatic_SNV_genehancer_all_table,
                somatic_SNV_chipseq_all_table = somatic_SNV_chipseq_all_table,

                # commented out to show difference in if and else statement
                # somatic_SNV_tf_binding_data_genes_cis_activated_only_table = somatic_SNV_tf_binding_data_genes_cis_activated_only_table,
                # somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table = somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table,
                # somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table = somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table,
                # somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table = somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table,
                # somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table = somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table,
                # cis_activation_summary_with_n_relevant_somatic_SNVs_table = cis_activation_summary_with_n_relevant_somatic_SNVs_table,

                genehancer_with_cis_activation_summary_table = genehancer_with_cis_activation_summary_table,
                chipseq_by_genes_table = chipseq_by_genes_table
            )
        )
    }
}