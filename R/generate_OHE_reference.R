convert_cohort_expression_list_to_table <- function(cohort_expression_list) {
    return(data.table::rbindlist(cohort_expression_list, idcol = "sample_id"))
}

convert_marker_summary_list_to_table <- function(cohort_marker_summary_list) {
    return(data.table::rbindlist(cohort_marker_summary_list, idcol = "sample_id"))
}

add_is_biallelic_column_to_expression_table <- function(cohort_expression_table, cohort_marker_summary_table) {
    cohort_expression_table[cohort_marker_summary_table, is_biallelic := i.is_biallelic, on = c("gene_name", "sample_id")]
}

apply_LOO_testing_to_expression_table <- function(expression_table, threshold_min_biallelic_samples = 10) {
    if (threshold_min_biallelic_samples < 3) {
        warning("threshold_min_biallelic_samples is set to less than 3 => using threshold_min_biallelic_samples = 3")
        threshold_min_biallelic_samples_safe <- 3
    } else {
        threshold_min_biallelic_samples_safe <- threshold_min_biallelic_samples
    }

    compute_LOO_tvalue <- function(x_j, y_j) {
        # corrected n is not length of y.j but length of y.j + 1
        return((x_j - mean(y_j)) / ((1 + (length(y_j) + 1 - 2)^-1) * (sd(y_j)^2))^0.5)
    }

    compute_LOO_pvalue <- function(t, n) {
        return(stats::pt(t, n, lower.tail = F))
    }

    compute_LOO_tvalue_others <- function(dbl_vector) {
        n <- length(dbl_vector)
        if (n == 1) {
            return(NA_real_)
        }
        returned_vector <- vector(mode = "double", length = n)
        for (i in seq_len(n)) {
            returned_vector[i] <- compute_LOO_tvalue(dbl_vector[i], dbl_vector[-i])
        }
        return(returned_vector)
    }

    expression_table[
        # count bialleic expressed samples
        , n_biallelic_samples := sum(is_biallelic == TRUE),
        by = gene_name
    ]
    expression_table[
        # apply LOO testing only to biallelic samples
        is_biallelic == TRUE &
            # require at least `threshold_min_biallelic_samples` samples
            n_biallelic_samples >= threshold_min_biallelic_samples_safe,
        # compute tvalue
        LOO_tvalue := compute_LOO_tvalue_others(FPKM),
        by = gene_name
    ][
        # compute pvalue
        , LOO_pvalue := compute_LOO_pvalue(LOO_tvalue, n_biallelic_samples - 2),
        by = gene_name
    ][
        ### COPY NUMBER NORMALIZED
        # apply LOO testing only to biallelic samples
        is_biallelic == TRUE &
            # require at least `threshold_min_biallelic_samples` samples
            n_biallelic_samples >= threshold_min_biallelic_samples_safe,
        # compute tvalue
        LOO_tvalue_copy_number_normalized := compute_LOO_tvalue_others(FPKM_copy_number_normalized),
        by = gene_name
    ][
        ### COPY NUMBER NORMALIZED
        # compute pvalue
        , LOO_pvalue_copy_number_normalized := compute_LOO_pvalue(LOO_tvalue_copy_number_normalized, n_biallelic_samples - 2),
        by = gene_name
    ]

    return(expression_table)
}

apply_ref_filter <- function(cohort_expression_table, min_LOO_pvalue_threshold = 0.05) {
    cohort_expression_table[
        , passed_ref_filter := ((LOO_pvalue >= min_LOO_pvalue_threshold) & is_biallelic)
    ][
        , passed_ref_filter_copy_number_normalized := ((LOO_pvalue_copy_number_normalized >= min_LOO_pvalue_threshold) & is_biallelic)
    ]
    return(cohort_expression_table)
}

generate_OHE_reference <- function(cohort_expression_table, cohort_marker_summary_table, threshold_min_biallelic_samples = 10) {
    # this function updates the cohort_expression_table by reference
    # no assignment of the returned value is needed
    add_is_biallelic_column_to_expression_table(cohort_expression_table, cohort_marker_summary_table)
    apply_LOO_testing_to_expression_table(cohort_expression_table, threshold_min_biallelic_samples = threshold_min_biallelic_samples)
    apply_ref_filter(cohort_expression_table)
    return(cohort_expression_table)
}