calculate_OHE <- function(cohort_expression_table, SAMPLE_ID) {
    compute_LOO_tvalue <- function(x_j, y_j) {
        # corrected n is not length of y.j but length of y.j + 1
        return((x_j - mean(y_j)) / ((1 + (length(y_j) + 1 - 2)^-1) * (sd(y_j)^2))^0.5)
    }

    compute_LOO_tvalue <- function(x_j, y_j) {
        # corrected n is not length of y.j but length of y.j + 1
        return((x_j - mean(y_j, na.rm = TRUE)) / ((1 + (sum(!is.na(y_j)) + 1 - 2)^-1) * (sd(y_j, na.rm = TRUE)^2))^0.5)
    }

    compute_LOO_pvalue <- function(t, n) {
        return(stats::pt(t, n, lower.tail = F))
    }

    compute_LOO_tvalue_of_sample <- function(dbl_vector, index_of_sample_in_dbl_vector) {
        n <- length(dbl_vector)
        if (index_of_sample_in_dbl_vector > n) {
            stop("index of the sample not within dbl_vector")
        }
        returned_vector <- rep(NA_real_, n)
        n_without_na <- sum(!is.na(dbl_vector))
        if (n_without_na <= 3) {
            # at least 3 samples are required to calculate LOO
            return(returned_vector)
        }
        returned_vector[index_of_sample_in_dbl_vector] <-
            compute_LOO_tvalue(
                dbl_vector[index_of_sample_in_dbl_vector],
                dbl_vector[-index_of_sample_in_dbl_vector]
            )

        return(returned_vector)
    }

    compute_LOO_tvalue_of_sample_with_filtered_reference <- function(dbl_vector, index_of_sample_in_dbl_vector, filter_vector) {
        n <- length(dbl_vector)
        if (index_of_sample_in_dbl_vector > n) {
            stop("index of the sample not within dbl_vector")
        }
        returned_vector <- rep(NA_real_, n)
        filtered_dbl_vector_without_index_sample <- dbl_vector[-index_of_sample_in_dbl_vector][filter_vector[-index_of_sample_in_dbl_vector]]
        n_ref_filtered_without_na <- sum(!is.na(filtered_dbl_vector_without_index_sample))
        if (n_ref_filtered_without_na <= 2) {
            # at least 2 filtered samples are required to calculate LOO
            return(returned_vector)
        }
        returned_vector[index_of_sample_in_dbl_vector] <-
            compute_LOO_tvalue(
                dbl_vector[index_of_sample_in_dbl_vector],
                filtered_dbl_vector_without_index_sample
            )

        return(returned_vector)
    }

    cohort_expression_table[
        , c("n_used_in_ref_all", "n_used_in_ref_filtered", "OHE_tvalue_ref_all", "OHE_tvalue_ref_filtered") := list(
            sum(
                !(sample_id == SAMPLE_ID)
            ),
            sum(
                ((!(sample_id == SAMPLE_ID)) &
                    passed_ref_filter &
                    (!is.na(passed_ref_filter)))
            ),
            compute_LOO_tvalue_of_sample(FPKM, which(sample_id == SAMPLE_ID)[1]),
            compute_LOO_tvalue_of_sample_with_filtered_reference(FPKM, which(sample_id == SAMPLE_ID)[1], passed_ref_filter & !is.na(passed_ref_filter))
        ),
        by = gene_name
    ][
        , c("OHE_pvalue_ref_all", "OHE_pvalue_ref_filtered") := list(
            # n-2 degrees of freedom -> n here is n_used_in_ref + 1 (for the compared sample)
            compute_LOO_pvalue(OHE_tvalue_ref_all, n_used_in_ref_all - 1),
            compute_LOO_pvalue(OHE_tvalue_ref_filtered, n_used_in_ref_filtered - 1)
        )
    ][
        ### COPY NUMBER NORMALIZED
        , c("n_used_in_ref_all_copy_number_normalized", "n_used_in_ref_filtered_copy_number_normalized", "OHE_tvalue_ref_all_copy_number_normalized", "OHE_tvalue_ref_filtered_copy_number_normalized") := list(
            sum(
                !(sample_id == SAMPLE_ID)
            ),
            sum(
                ((!(sample_id == SAMPLE_ID)) &
                    passed_ref_filter_copy_number_normalized &
                    (!is.na(passed_ref_filter_copy_number_normalized)))
            ),
            compute_LOO_tvalue_of_sample(FPKM_copy_number_normalized, which(sample_id == SAMPLE_ID)[1]),
            compute_LOO_tvalue_of_sample_with_filtered_reference(FPKM_copy_number_normalized, which(sample_id == SAMPLE_ID)[1], passed_ref_filter_copy_number_normalized & !is.na(passed_ref_filter_copy_number_normalized))
        ),
        by = gene_name
    ][
        ### COPY NUMBER NORMALIZED
        , c("OHE_pvalue_ref_all_copy_number_normalized", "OHE_pvalue_ref_filtered_copy_number_normalized") := list(
            # n-2 degrees of freedom -> n here is n_used_in_ref + 1 (for the compared sample)
            compute_LOO_pvalue(OHE_tvalue_ref_all_copy_number_normalized, n_used_in_ref_all_copy_number_normalized - 1),
            compute_LOO_pvalue(OHE_tvalue_ref_filtered_copy_number_normalized, n_used_in_ref_filtered_copy_number_normalized - 1)
        )
    ]

    return(
        cohort_expression_table[
            sample_id == SAMPLE_ID
        ]
    )
}