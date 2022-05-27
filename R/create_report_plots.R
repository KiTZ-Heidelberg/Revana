plot_subgroups_bar_chart <- function(processed_data, color_palette) {
    ggplot2::ggplot(processed_data$result_paths) +
        ggplot2::geom_bar(ggplot2::aes(x = subgroup, fill = subgroup)) +
        ggplot2::labs(x = "Subgroup", y = "Number of samples", title = "Number of samples per subgroup") +
        color_palette$scale_fill_subgroup +
        ggplot2::theme(legend.position = "none")
}

plot_markers_bar_chart_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$marker_data_table %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_infreq(sample_ID)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_bar(ggplot2::aes(y = sample_ID, x = n_markers, fill = subgroup), stat = "identity") +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("Markers per sample - ", subgroup))
}

plot_somatic_SNV_bar_chart_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$somatic_SNV_data_table %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_infreq(sample_ID)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_bar(ggplot2::aes(y = sample_ID, fill = subgroup)) +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("Somatic SNVs per sample - ", subgroup))
}

plot_SV_bar_chart_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$SV_data_table %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_infreq(sample_ID)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_bar(ggplot2::aes(y = sample_ID, fill = subgroup)) +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("SVs per sample - ", subgroup))
}

plot_CNA_bar_chart_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$CNA_data_table %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_infreq(sample_ID)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_bar(ggplot2::aes(y = sample_ID, fill = subgroup)) +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("CNAs per sample - ", subgroup))
}

plot_expression_box_plot_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$expression_data_table %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_reorder(sample_ID, -FPKM)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_boxplot(ggplot2::aes(x = sample_ID, y = FPKM, fill = subgroup)) +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("Gene Expression per sample - ", subgroup)) +
        ggplot2::scale_y_continuous(trans = "pseudo_log") +
        ggplot2::coord_flip()
}

plot_expression_PCA_plot_by_subgroup <- function(processed_data, subgroup, color_palette) {
    expression_PCA_data_bygroup <-
        processed_data$expression_data_table %>%
        dplyr::distinct(sample_ID, gene_name, .keep_all = TRUE) %>% # replace after process_expression bug is fixed
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        tidyr::replace_na(list(FPKM = 0)) %>%
        # dplyr::mutate(FPKM_log_trans = log(FPKM+0.001)) %>%
        tidyr::pivot_wider(id_cols = "sample_ID", names_from = "gene_name", values_from = "FPKM") %>%
        tibble::column_to_rownames("sample_ID")

    expression_PCA_data_bygroup_NO_0_VARIANCE <- expression_PCA_data_bygroup[, which(apply(expression_PCA_data_bygroup, 2, var) != 0)]


    expression_PCA_by_group <- prcomp(expression_PCA_data_bygroup_NO_0_VARIANCE, scale = TRUE)
    var_explained <- expression_PCA_by_group$sdev^2 / sum(expression_PCA_by_group$sdev^2)

    expression_PCA_plot_by_group <- expression_PCA_by_group$x %>%
        as.data.frame() %>%
        tibble::rownames_to_column("sample_ID") %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, color = subgroup)) +
        color_palette$scale_color_subgroup +
        ggplot2::geom_text(ggplot2::aes(x = PC1, y = PC2, label = sample_ID), hjust = 0, vjust = 0) +
        ggplot2::labs(
            x = paste0("PC1: ", round(var_explained[1] * 100, 1), "%"),
            y = paste0("PC2: ", round(var_explained[2] * 100, 1), "%")
        )

    return(expression_PCA_plot_by_group)
}

plot_cis_activated_genes_bar_chart_by_subgroup <- function(processed_data, subgroup, color_palette) {
    processed_data$cis_activated_genes_by_sample %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        dplyr::mutate(sample_ID = forcats::fct_reorder(sample_ID, n_cis_activated_genes)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_bar(ggplot2::aes(y = sample_ID, x = n_cis_activated_genes, fill = subgroup), stat = "identity") +
        color_palette$scale_fill_subgroup +
        ggplot2::labs(title = paste0("Cis activated genes per sample - ", subgroup))
}

plot_ase_detection_upset_plot_by_subgroup <- function(processed_data, subgroup) {
    processed_data$cis_activation_summary_table_for_upsetr %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
        dplyr::mutate(ones1 = 1, ones2 = 1) %>%
        dplyr::select("set_id", "sample_ID", "ones1", "cis_activated_gene", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap", "ones2") %>%
        UpSetR::upset(
            sets = c("cis_activated_gene", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap")
        )
}

plot_ase_detection_upset_plot_whole_cohort <- function(processed_data) {
    processed_data$cis_activation_summary_table_for_upsetr %>%
        # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
        dplyr::mutate(ones1 = 1, ones2 = 1) %>%
        dplyr::select("set_id", "sample_ID", "ones1", "cis_activated_gene", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap", "ones2") %>%
        UpSetR::upset(
            sets = c("cis_activated_gene", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap")
        )
}

plot_OHE_ASE_q_value_dot_plot <- function(processed_data, color_palette) {
    processed_data$cis_activation_summary_table_mut_info %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = -log10(OHE_used_pvalue), y = -log10(ABH), color = subgroup)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(0.05)), color = "red", linetype = "dashed") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed") +
        ggplot2::labs(title = "OHE and ASE test q-values") +
        color_palette$scale_color_subgroup
}

plot_applied_filters_upset_plot_by_subgroup <- function(processed_data, subgroup) {
    processed_data$cis_activation_summary_table_for_upsetr %>%
        dplyr::filter(sample_ID %in% processed_data$samples_by_subgroup[[subgroup]]) %>%
        # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
        dplyr::mutate(ones1 = 1, ones2 = 1) %>%
        dplyr::select("set_id", "sample_ID", "ones1", "passed_FPKM_filter", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene", "ones2") %>%
        UpSetR::upset(
            sets = c("passed_FPKM_filter", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene"),
            nintersects = NA
        )
}

plot_applied_filters_upset_plot_whole_cohort <- function(processed_data) {
    applied_filters_upset_plot_whole_cohort_plot <- processed_data$cis_activation_summary_table_for_upsetr %>%
        # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
        dplyr::mutate(ones1 = 1, ones2 = 1) %>%
        dplyr::select("set_id", "sample_ID", "ones1", "passed_FPKM_filter", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene", "ones2") %>%
        UpSetR::upset(
            sets = c("passed_FPKM_filter", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene"),
            nintersects = NA
        )
}

save_ROC_curve_plot <- function(rocit_object, path) {
    svglite::svglite(path)
    plot(rocit_object)
    dev.off()
}

plot_sensitivity_and_specificity_curve <- function(rocit_object) {
    tibble::tibble(threshold = rocit_object$Cutoff, TPR = rocit_object$TPR, `1-FPR` = 1 - rocit_object$FPR) %>%
        dplyr::mutate(youdenIndex = TPR + `1-FPR` - 1) %>%
        tidyr::pivot_longer(cols = c(TPR, `1-FPR`, youdenIndex), names_to = "measure", values_to = "value") %>%
        ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(y = value, x = threshold, color = measure)) +
        ggplot2::theme_linedraw()
}