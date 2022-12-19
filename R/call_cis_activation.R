merge_OHE_with_marker_summary <- function(OHE_table, marker_summary_table) {
    # convert to data table
    OHE_dt <- data.table::setDT(OHE_table)
    marker_summary_dt <- data.table::setDT(marker_summary_table)

    # add available data column
    OHE_dt[, OHE_data_available := TRUE]
    marker_summary_dt[, marker_summary_data_available := TRUE]

    merged_table <- merge(OHE_dt, marker_summary_dt, by = "gene_name", all = TRUE)

    merged_table[is.na(OHE_data_available), OHE_data_available := FALSE]
    merged_table[is.na(marker_summary_data_available), marker_summary_data_available := FALSE]

    return(merged_table)
}

add_cis_expression_filter_columns <- function(merged_table,
                                              threshold_min_FPKM_normal_gene = 5,
                                              threshold_min_FPKM_oncogene = 1,
                                              threshold_max_ASE_ABH = 0.05,
                                              rescue_oncogenes_ASE_filter = TRUE,
                                              threshold_max_ASE_pvalue_for_rescued_oncogenes = 0.05,
                                              threshold_min_mean_delta_abs_diploid = 0.3,
                                              threshold_min_mean_delta_abs_cnv = 0.2,
                                              cnv_threshold_percentage_markers_cnv = 0.3,
                                              threshold_max_pvalue_OHE = 0.05 # -> check if pvalue is alright or if q value has to be generated
) {

    # FPKM filter ------------------------------
    merged_table[, passed_FPKM_filter := ifelse(
        is_cancer_gene,
        FPKM >= threshold_min_FPKM_oncogene,
        FPKM >= threshold_min_FPKM_normal_gene
    )]

    # ASE filter ------------------------------

        # multtest ABH is invalid sometimes (-> see calculate_ASE.R)
        # under these circumstances use the more conservative Bonferroni

    merged_table[, passed_ASE_filter := ifelse(is.na(ABH), Bonferroni <= threshold_max_ASE_ABH, ABH <= threshold_max_ASE_ABH) ]
    # data.table::setnafill(merged_table, type="const", fill=FALSE, nan=NA, cols=c("passed_ASE_filter"))

    merged_table[, ASE_filter_rescued_oncogene := FALSE]
    if (rescue_oncogenes_ASE_filter) {
        merged_table[
            is_cancer_gene == TRUE,
            ASE_filter_rescued_oncogene := combined_p_value <= threshold_max_ASE_pvalue_for_rescued_oncogenes
        ]
    }

    # mean_delta_abs filter ------------------------------
    merged_table[, passed_mean_delta_abs_filter := mean_delta_abs > ifelse(
        # gene considered as cnv with at least `cnv_threshold_percentage_markers_cnv` (default = 30%) cnv markers
        n_cnv_markers / n_markers >= cnv_threshold_percentage_markers_cnv,
        threshold_min_mean_delta_abs_cnv,
        threshold_min_mean_delta_abs_diploid
    )]

    # OHE filter ------------------------------
    merged_table[, OHE_used_reference := ifelse(is.na(OHE_pvalue_ref_filtered), "ref_all", "ref_filtered")]
    merged_table[, OHE_used_pvalue := ifelse(is.na(OHE_pvalue_ref_filtered), OHE_pvalue_ref_all, OHE_pvalue_ref_filtered)]
    merged_table[, passed_OHE_filter := OHE_used_pvalue <= threshold_max_pvalue_OHE]

    ### COPY NUMBER NORMALIZED
    merged_table[, OHE_used_reference_copy_number_normalized := ifelse(is.na(OHE_pvalue_ref_filtered_copy_number_normalized), "ref_all", "ref_filtered")]
    merged_table[, OHE_used_pvalue_copy_number_normalized := ifelse(is.na(OHE_pvalue_ref_filtered_copy_number_normalized), OHE_pvalue_ref_all_copy_number_normalized, OHE_pvalue_ref_filtered_copy_number_normalized)]
    merged_table[, passed_OHE_filter_copy_number_normalized := OHE_used_pvalue_copy_number_normalized <= threshold_max_pvalue_OHE]

    # combine filters ------------------------------
    merged_table[, cis_activated_gene := (
        (imprinting_status == "no_imprinting") & # this also excludes all doubtfully imprinted genes
            passed_FPKM_filter &
            (passed_ASE_filter | ASE_filter_rescued_oncogene | ASE_marker_run_overlap) &
            passed_mean_delta_abs_filter &
            passed_OHE_filter
    )][
        , cis_activated_gene := (cis_activated_gene & (!is.na(cis_activated_gene)))
    ][
        ### COPY NUMBER NORMALIZED
        , cis_activated_gene_copy_number_normalized := (
            (imprinting_status == "no_imprinting") & # this also excludes all doubtfully imprinted genes
                passed_FPKM_filter &
                (passed_ASE_filter | ASE_filter_rescued_oncogene | ASE_marker_run_overlap) &
                passed_mean_delta_abs_filter &
                passed_OHE_filter_copy_number_normalized
        )
    ][
        , cis_activated_gene_copy_number_normalized := (cis_activated_gene_copy_number_normalized & (!is.na(cis_activated_gene_copy_number_normalized)))
    ]

    return(merged_table)
}


call_cis_activation <- function(OHE_table,
                                marker_summary_table,
                                threshold_min_FPKM_normal_gene = 5,
                                threshold_min_FPKM_oncogene = 1,
                                threshold_max_ASE_ABH = 0.05,
                                rescue_oncogenes_ASE_filter = TRUE,
                                threshold_max_ASE_pvalue_for_rescued_oncogenes = 0.05,
                                threshold_min_mean_delta_abs_diploid = 0.3,
                                threshold_min_mean_delta_abs_cnv = 0.2,
                                cnv_threshold_percentage_markers_cnv = 0.3,
                                threshold_max_pvalue_OHE = 0.05) {
    merged_table <- merge_OHE_with_marker_summary(OHE_table, marker_summary_table)
    add_cis_expression_filter_columns(
        merged_table,
        threshold_min_FPKM_normal_gene,
        threshold_min_FPKM_oncogene,
        threshold_max_ASE_ABH,
        rescue_oncogenes_ASE_filter,
        threshold_max_ASE_pvalue_for_rescued_oncogenes,
        threshold_min_mean_delta_abs_diploid,
        threshold_min_mean_delta_abs_cnv,
        cnv_threshold_percentage_markers_cnv,
        threshold_max_pvalue_OHE
    )

    return(merged_table)
}