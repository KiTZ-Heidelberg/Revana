#' Run Step 5 of the Revana workflow on a cohort/subgroup of tumor sample data
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param subgroup_name Name of the subgroup. The name of the subgroup is used as reference in later created HTML reports
#' @param skip_merging_Rds_chunks Logical value. Should data chunks stored in .Rds-files be merged together. Default is TRUE - this step is required for regular report generation.
#' @param run_tf_binding_site_analysis Was TF binding site analysis conducted?
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#' @details
#' Step 5 of the revana workflow ...
#' * merges .Rds file chunks (if skip_merging_Rds_chunks == TRUE)
#' * creates results paths file
#'
#' @export


run_step5_per_cohort <- function(paths_file_path,
                                 output_dir,
                                 subgroup_name = "DEFAULT_GROUP",
                                 skip_merging_Rds_chunks = FALSE,
                                 run_tf_binding_site_analysis = TRUE,
                                 verbose = FALSE) {
    # import paths file  -------------------------------------
    paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

    # merge RDS chunks  --------------------------------------
    sample_ids <- paths$sample_id

    if (!skip_merging_Rds_chunks) {
        sample_ids %>% purrr::map(
            function(s_id) {
                sample_dir_path <- file.path(output_dir, s_id)

                # TF binding data
                if (run_tf_binding_site_analysis) {

                    # somatic_SNV_tf_binding_data.genehancer.all.Rds
                    merge_Rds_chunks(
                        pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genehancer\\.all\\.", "__chunk_[[:digit:]]+\\.Rds"),
                        folder = sample_dir_path,
                        output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genehancer.all.Rds")),
                        verbose = verbose
                    )
                    # somatic_SNV_tf_binding_data.genehancer.cis_activated_only.Rds
                    merge_Rds_chunks(
                        pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genehancer\\.cis_activated_only\\.", "__chunk_[[:digit:]]+\\.Rds"),
                        folder = sample_dir_path,
                        output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.Rds")),
                        verbose = verbose
                    )
                    # somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.Rds
                    merge_Rds_chunks(
                        pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genehancer\\.cis_activated_only\\.relevant_only\\.", "__chunk_[[:digit:]]+\\.Rds"),
                        folder = sample_dir_path,
                        output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.Rds")),
                        verbose = verbose
                    )

                    # # somatic_SNV_tf_binding_data.genes.all.Rds
                    # merge_Rds_chunks(
                    #     pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genes\\.all\\.", "__chunk_[[:digit:]]+\\.Rds"),
                    #     folder = sample_dir_path,
                    #     output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genes.all.Rds"))
                    #     )
                    # somatic_SNV_tf_binding_data.genes.cis_activated_only.Rds
                    merge_Rds_chunks(
                        pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genes\\.cis_activated_only\\.", "__chunk_[[:digit:]]+\\.Rds"),
                        folder = sample_dir_path,
                        output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.Rds")),
                        verbose = verbose
                    )
                    # somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.Rds
                    merge_Rds_chunks(
                        pattern = paste0(s_id, "\\.somatic_SNV_tf_binding_data\\.genes\\.cis_activated_only\\.relevant_only\\.", "__chunk_[[:digit:]]+\\.Rds"),
                        folder = sample_dir_path,
                        output_file_path = file.path(sample_dir_path, paste0(s_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.Rds")),
                        verbose = verbose
                    )
                }
            }
        )
    }



    # create results path file  -------------------------------------
    results_path_file_path <- file.path(output_dir, "results.paths.txt")
    sample_dir_paths <- file.path(output_dir, sample_ids)
    create_file_path <- function(suffix = "") {
        return(file.path(sample_dir_paths, paste0(sample_ids, suffix)))
    }
    results_paths <- data.frame(
        sample_id = sample_ids,
        subgroup = rep_len(subgroup_name, length(sample_ids)),
        has_tf_binding_site_analysis = run_tf_binding_site_analysis,
        CNA_Rds = create_file_path(suffix = ".CNA.Rds"),
        CNA_chipseq_all_Rds = create_file_path(suffix = ".CNA.chipseq.all.Rds"),
        CNA_chipseq_all_txt = create_file_path(suffix = ".CNA.chipseq.all.txt"),
        CNA_chipseq_cis_activated_only_Rds = create_file_path(suffix = ".CNA.chipseq.cis_activated_only.Rds"),
        CNA_chipseq_cis_activated_only_txt = create_file_path(suffix = ".CNA.chipseq.cis_activated_only.txt"),
        CNA_genehancer_all_Rds = create_file_path(suffix = ".CNA.genehancer.all.Rds"),
        CNA_genehancer_all_txt = create_file_path(suffix = ".CNA.genehancer.all.txt"),
        CNA_genehancer_cis_activated_only_Rds = create_file_path(suffix = ".CNA.genehancer.cis_activated_only.Rds"),
        CNA_genehancer_cis_activated_only_txt = create_file_path(suffix = ".CNA.genehancer.cis_activated_only.txt"),
        CNA_genes_all_Rds = create_file_path(suffix = ".CNA.genes.all.Rds"),
        CNA_genes_all_txt = create_file_path(suffix = ".CNA.genes.all.txt"),
        CNA_genes_cis_activated_only_Rds = create_file_path(suffix = ".CNA.genes.cis_activated_only.Rds"),
        CNA_genes_cis_activated_only_txt = create_file_path(suffix = ".CNA.genes.cis_activated_only.txt"),
        CNA_txt = create_file_path(suffix = ".CNA.txt"),
        copy_number_Rds = create_file_path(suffix = ".copy_number.Rds"),
        copy_number_txt = create_file_path(suffix = ".copy_number.txt"),
        copy_number_by_gene_Rds = create_file_path(suffix = ".copy_number_by_gene.Rds"),
        copy_number_by_gene_txt = create_file_path(suffix = ".copy_number_by_gene.txt"),
        OHE_Rds = create_file_path(suffix = ".OHE.Rds"),
        OHE_txt = create_file_path(suffix = ".OHE.txt"),
        SV_Rds = create_file_path(suffix = ".SV.Rds"),
        SV_chipseq_all_Rds = create_file_path(suffix = ".SV.chipseq.all.Rds"),
        SV_chipseq_all_txt = create_file_path(suffix = ".SV.chipseq.all.txt"),
        SV_chipseq_cis_activated_only_Rds = create_file_path(suffix = ".SV.chipseq.cis_activated_only.Rds"),
        SV_chipseq_cis_activated_only_txt = create_file_path(suffix = ".SV.chipseq.cis_activated_only.txt"),
        SV_genehancer_all_Rds = create_file_path(suffix = ".SV.genehancer.all.Rds"),
        SV_genehancer_all_txt = create_file_path(suffix = ".SV.genehancer.all.txt"),
        SV_genehancer_cis_activated_only_Rds = create_file_path(suffix = ".SV.genehancer.cis_activated_only.Rds"),
        SV_genehancer_cis_activated_only_txt = create_file_path(suffix = ".SV.genehancer.cis_activated_only.txt"),
        SV_genes_all_Rds = create_file_path(suffix = ".SV.genes.all.Rds"),
        SV_genes_all_txt = create_file_path(suffix = ".SV.genes.all.txt"),
        SV_genes_cis_activated_only_Rds = create_file_path(suffix = ".SV.genes.cis_activated_only.Rds"),
        SV_genes_cis_activated_only_txt = create_file_path(suffix = ".SV.genes.cis_activated_only.txt"),
        SV_txt = create_file_path(suffix = ".SV.txt"),
        chipseq_with_cis_activation_summary_Rds = create_file_path(suffix = ".chipseq_with_cis_activation_summary.Rds"),
        chipseq_with_cis_activation_summary_txt = create_file_path(suffix = ".chipseq_with_cis_activation_summary.txt"),
        cis_activated_genes_Rds = create_file_path(suffix = ".cis.activated.genes.Rds"),
        cis_activated_genes_txt = create_file_path(suffix = ".cis.activated.genes.txt"),
        cis_activation_summary_Rds = create_file_path(suffix = ".cis.activation.summary.Rds"),
        cis_activation_summary_txt = create_file_path(suffix = ".cis.activation.summary.txt"),
        expression_Rds = create_file_path(suffix = ".expression.Rds"),
        expression_txt = create_file_path(suffix = ".expression.txt"),
        geneanno_marker_summary_Rds = create_file_path(suffix = ".geneanno_marker_summary.Rds"),
        geneanno_marker_summary_txt = create_file_path(suffix = ".geneanno_marker_summary.txt"),
        genehancer_with_cis_activation_summary_Rds = create_file_path(suffix = ".genehancer_with_cis_activation_summary.Rds"),
        genehancer_with_cis_activation_summary_txt = create_file_path(suffix = ".genehancer_with_cis_activation_summary.txt"),
        markers_Rds = create_file_path(suffix = ".markers.Rds"),
        markers_txt = create_file_path(suffix = ".markers.txt"),
        runs_Rds = create_file_path(suffix = ".runs.Rds"),
        runs_txt = create_file_path(suffix = ".runs.txt"),
        somatic_SNV_Rds = create_file_path(suffix = ".somatic_SNV.Rds"),
        somatic_SNV_chipseq_all_Rds = create_file_path(suffix = ".somatic_SNV.chipseq.all.Rds"),
        somatic_SNV_chipseq_all_txt = create_file_path(suffix = ".somatic_SNV.chipseq.all.txt"),
        somatic_SNV_chipseq_cis_activated_only_Rds = create_file_path(suffix = ".somatic_SNV.chipseq.cis_activated_only.Rds"),
        somatic_SNV_chipseq_cis_activated_only_txt = create_file_path(suffix = ".somatic_SNV.chipseq.cis_activated_only.txt"),
        somatic_SNV_genehancer_all_Rds = create_file_path(suffix = ".somatic_SNV.genehancer.all.Rds"),
        somatic_SNV_genehancer_all_txt = create_file_path(suffix = ".somatic_SNV.genehancer.all.txt"),
        somatic_SNV_genehancer_cis_activated_only_Rds = create_file_path(suffix = ".somatic_SNV.genehancer.cis_activated_only.Rds"),
        somatic_SNV_genehancer_cis_activated_only_txt = create_file_path(suffix = ".somatic_SNV.genehancer.cis_activated_only.txt"),
        somatic_SNV_genes_all_Rds = create_file_path(suffix = ".somatic_SNV.genes.all.Rds"),
        somatic_SNV_genes_all_txt = create_file_path(suffix = ".somatic_SNV.genes.all.txt"),
        somatic_SNV_genes_cis_activated_only_Rds = create_file_path(suffix = ".somatic_SNV.genes.cis_activated_only.Rds"),
        somatic_SNV_genes_cis_activated_only_txt = create_file_path(suffix = ".somatic_SNV.genes.cis_activated_only.txt"),
        somatic_SNV_txt = create_file_path(suffix = ".somatic_SNV.txt"),
        somatic_SNV_tf_binding_data_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.Rds"),
        somatic_SNV_tf_binding_data_genehancer_all_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.all.Rds"),
        somatic_SNV_tf_binding_data_genehancer_all_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.all.txt"),
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.Rds"),
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.txt"),
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.Rds"),
        somatic_SNV_tf_binding_data_genehancer_cis_activated_only_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.txt"),
        somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.Rds"),
        somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.txt"),
        somatic_SNV_tf_binding_data_genes_cis_activated_only_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genes.cis_activated_only.Rds"),
        somatic_SNV_tf_binding_data_genes_cis_activated_only_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.genes.cis_activated_only.txt"),
        somatic_SNV_tf_binding_data_txt = create_file_path(suffix = ".somatic_SNV_tf_binding_data.txt"),
        cis_activation_summary_with_n_relevant_somatic_SNVs_Rds = create_file_path(suffix = ".cis_activation_summary_with_n_relevant_somatic_SNVs.Rds"),
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data_gene_summary_max_score_diff.Rds"),
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff.Rds"),
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_Rds = create_file_path(suffix = ".somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff.Rds"),
        chipseq_by_gene_Rds = file.path(output_dir, "chipseq_by_gene.Rds")
    )

    data.table::fwrite(
        results_paths,
        file = results_path_file_path,
        sep = "\t"
    )
}