create_mut_analysis_HTML <- function(processed_data, has_run_tf_binding_site_analysis, HTML_report_figure_directory) {
    # CREATE HTML FOR PART III ####################################

    # any_muts
    contingency_table_any_muts <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_muts
    )

    cis_activation_frequency_WITHOUT_any_muts <- contingency_table_any_muts[2, 1] / (contingency_table_any_muts[1, 1] + contingency_table_any_muts[2, 1])
    cis_activation_frequency_WITH_any_muts <- contingency_table_any_muts[2, 2] / (contingency_table_any_muts[1, 2] + contingency_table_any_muts[2, 2])

    chi_square_test_any_muts <- chisq.test(contingency_table_any_muts)
    odds_ratio_any_muts <- (contingency_table_any_muts[2, 2] / contingency_table_any_muts[1, 2]) / (contingency_table_any_muts[2, 1] / contingency_table_any_muts[1, 1])

    # any_muts_genehancer
    contingency_table_any_muts_genehancer <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_muts_from_genehancer_TAD
    )

    cis_activation_frequency_WITHOUT_any_muts_genehancer <- contingency_table_any_muts_genehancer[2, 1] / (contingency_table_any_muts_genehancer[1, 1] + contingency_table_any_muts_genehancer[2, 1])
    cis_activation_frequency_WITH_any_muts_genehancer <- contingency_table_any_muts_genehancer[2, 2] / (contingency_table_any_muts_genehancer[1, 2] + contingency_table_any_muts_genehancer[2, 2])

    chi_square_test_any_muts_genehancer <- chisq.test(contingency_table_any_muts_genehancer)
    odds_ratio_any_muts_genehancer <- (contingency_table_any_muts_genehancer[2, 2] / contingency_table_any_muts_genehancer[1, 2]) / (contingency_table_any_muts_genehancer[2, 1] / contingency_table_any_muts_genehancer[1, 1])

    # any_muts_gene
    contingency_table_any_muts_gene <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_muts_from_gene_TAD
    )

    cis_activation_frequency_WITHOUT_any_muts_gene <- contingency_table_any_muts_gene[2, 1] / (contingency_table_any_muts_gene[1, 1] + contingency_table_any_muts_gene[2, 1])
    cis_activation_frequency_WITH_any_muts_gene <- contingency_table_any_muts_gene[2, 2] / (contingency_table_any_muts_gene[1, 2] + contingency_table_any_muts_gene[2, 2])

    chi_square_test_any_muts_gene <- chisq.test(contingency_table_any_muts_gene)
    odds_ratio_any_muts_gene <- (contingency_table_any_muts_gene[2, 2] / contingency_table_any_muts_gene[1, 2]) / (contingency_table_any_muts_gene[2, 1] / contingency_table_any_muts_gene[1, 1])


    # any_muts_chipseq
    contingency_table_any_muts_chipseq <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_muts_with_chipseq
    )

    cis_activation_frequency_WITHOUT_any_muts_chipseq <- contingency_table_any_muts_chipseq[2, 1] / (contingency_table_any_muts_chipseq[1, 1] + contingency_table_any_muts_chipseq[2, 1])
    cis_activation_frequency_WITH_any_muts_chipseq <- contingency_table_any_muts_chipseq[2, 2] / (contingency_table_any_muts_chipseq[1, 2] + contingency_table_any_muts_chipseq[2, 2])

    chi_square_test_any_muts_chipseq <- chisq.test(contingency_table_any_muts_chipseq)
    odds_ratio_any_muts_chipseq <- (contingency_table_any_muts_chipseq[2, 2] / contingency_table_any_muts_chipseq[1, 2]) / (contingency_table_any_muts_chipseq[2, 1] / contingency_table_any_muts_chipseq[1, 1])

    # any_CNAs
    contingency_table_any_CNAs <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_CNAs
    )

    cis_activation_frequency_WITHOUT_any_CNAs <- contingency_table_any_CNAs[2, 1] / (contingency_table_any_CNAs[1, 1] + contingency_table_any_CNAs[2, 1])
    cis_activation_frequency_WITH_any_CNAs <- contingency_table_any_CNAs[2, 2] / (contingency_table_any_CNAs[1, 2] + contingency_table_any_CNAs[2, 2])

    chi_square_test_any_CNAs <- chisq.test(contingency_table_any_CNAs)
    odds_ratio_any_CNAs <- (contingency_table_any_CNAs[2, 2] / contingency_table_any_CNAs[1, 2]) / (contingency_table_any_CNAs[2, 1] / contingency_table_any_CNAs[1, 1])

    # any_CNAs_genehancer
    contingency_table_any_CNAs_genehancer <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_CNAs_from_genehancer_TAD
    )

    cis_activation_frequency_WITHOUT_any_CNAs_genehancer <- contingency_table_any_CNAs_genehancer[2, 1] / (contingency_table_any_CNAs_genehancer[1, 1] + contingency_table_any_CNAs_genehancer[2, 1])
    cis_activation_frequency_WITH_any_CNAs_genehancer <- contingency_table_any_CNAs_genehancer[2, 2] / (contingency_table_any_CNAs_genehancer[1, 2] + contingency_table_any_CNAs_genehancer[2, 2])

    chi_square_test_any_CNAs_genehancer <- chisq.test(contingency_table_any_CNAs_genehancer)
    odds_ratio_any_CNAs_genehancer <- (contingency_table_any_CNAs_genehancer[2, 2] / contingency_table_any_CNAs_genehancer[1, 2]) / (contingency_table_any_CNAs_genehancer[2, 1] / contingency_table_any_CNAs_genehancer[1, 1])

    # any_CNAs_gene
    contingency_table_any_CNAs_gene <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_CNAs_from_gene_TAD
    )

    cis_activation_frequency_WITHOUT_any_CNAs_gene <- contingency_table_any_CNAs_gene[2, 1] / (contingency_table_any_CNAs_gene[1, 1] + contingency_table_any_CNAs_gene[2, 1])
    cis_activation_frequency_WITH_any_CNAs_gene <- contingency_table_any_CNAs_gene[2, 2] / (contingency_table_any_CNAs_gene[1, 2] + contingency_table_any_CNAs_gene[2, 2])

    chi_square_test_any_CNAs_gene <- chisq.test(contingency_table_any_CNAs_gene)
    odds_ratio_any_CNAs_gene <- (contingency_table_any_CNAs_gene[2, 2] / contingency_table_any_CNAs_gene[1, 2]) / (contingency_table_any_CNAs_gene[2, 1] / contingency_table_any_CNAs_gene[1, 1])

    # any_CNAs_chipseq
    contingency_table_any_CNAs_chipseq <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_CNAs_with_chipseq
    )

    cis_activation_frequency_WITHOUT_any_CNAs_chipseq <- contingency_table_any_CNAs_chipseq[2, 1] / (contingency_table_any_CNAs_chipseq[1, 1] + contingency_table_any_CNAs_chipseq[2, 1])
    cis_activation_frequency_WITH_any_CNAs_chipseq <- contingency_table_any_CNAs_chipseq[2, 2] / (contingency_table_any_CNAs_chipseq[1, 2] + contingency_table_any_CNAs_chipseq[2, 2])

    chi_square_test_any_CNAs_chipseq <- chisq.test(contingency_table_any_CNAs_chipseq)
    odds_ratio_any_CNAs_chipseq <- (contingency_table_any_CNAs_chipseq[2, 2] / contingency_table_any_CNAs_chipseq[1, 2]) / (contingency_table_any_CNAs_chipseq[2, 1] / contingency_table_any_CNAs_chipseq[1, 1])

    # any_SVs
    contingency_table_any_SVs <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_SVs
    )

    cis_activation_frequency_WITHOUT_any_SVs <- contingency_table_any_SVs[2, 1] / (contingency_table_any_SVs[1, 1] + contingency_table_any_SVs[2, 1])
    cis_activation_frequency_WITH_any_SVs <- contingency_table_any_SVs[2, 2] / (contingency_table_any_SVs[1, 2] + contingency_table_any_SVs[2, 2])

    chi_square_test_any_SVs <- chisq.test(contingency_table_any_SVs)
    odds_ratio_any_SVs <- (contingency_table_any_SVs[2, 2] / contingency_table_any_SVs[1, 2]) / (contingency_table_any_SVs[2, 1] / contingency_table_any_SVs[1, 1])

    # any_SVs_genehancer
    contingency_table_any_SVs_genehancer <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_SVs_from_genehancer_TAD
    )

    cis_activation_frequency_WITHOUT_any_SVs_genehancer <- contingency_table_any_SVs_genehancer[2, 1] / (contingency_table_any_SVs_genehancer[1, 1] + contingency_table_any_SVs_genehancer[2, 1])
    cis_activation_frequency_WITH_any_SVs_genehancer <- contingency_table_any_SVs_genehancer[2, 2] / (contingency_table_any_SVs_genehancer[1, 2] + contingency_table_any_SVs_genehancer[2, 2])

    chi_square_test_any_SVs_genehancer <- chisq.test(contingency_table_any_SVs_genehancer)
    odds_ratio_any_SVs_genehancer <- (contingency_table_any_SVs_genehancer[2, 2] / contingency_table_any_SVs_genehancer[1, 2]) / (contingency_table_any_SVs_genehancer[2, 1] / contingency_table_any_SVs_genehancer[1, 1])

    # any_SVs_gene
    contingency_table_any_SVs_gene <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_SVs_from_gene_TAD
    )

    cis_activation_frequency_WITHOUT_any_SVs_gene <- contingency_table_any_SVs_gene[2, 1] / (contingency_table_any_SVs_gene[1, 1] + contingency_table_any_SVs_gene[2, 1])
    cis_activation_frequency_WITH_any_SVs_gene <- contingency_table_any_SVs_gene[2, 2] / (contingency_table_any_SVs_gene[1, 2] + contingency_table_any_SVs_gene[2, 2])

    chi_square_test_any_SVs_gene <- chisq.test(contingency_table_any_SVs_gene)
    odds_ratio_any_SVs_gene <- (contingency_table_any_SVs_gene[2, 2] / contingency_table_any_SVs_gene[1, 2]) / (contingency_table_any_SVs_gene[2, 1] / contingency_table_any_SVs_gene[1, 1])

    # any_SVs_chipseq
    contingency_table_any_SVs_chipseq <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_SVs_with_chipseq
    )

    cis_activation_frequency_WITHOUT_any_SVs_chipseq <- contingency_table_any_SVs_chipseq[2, 1] / (contingency_table_any_SVs_chipseq[1, 1] + contingency_table_any_SVs_chipseq[2, 1])
    cis_activation_frequency_WITH_any_SVs_chipseq <- contingency_table_any_SVs_chipseq[2, 2] / (contingency_table_any_SVs_chipseq[1, 2] + contingency_table_any_SVs_chipseq[2, 2])

    chi_square_test_any_SVs_chipseq <- chisq.test(contingency_table_any_SVs_chipseq)
    odds_ratio_any_SVs_chipseq <- (contingency_table_any_SVs_chipseq[2, 2] / contingency_table_any_SVs_chipseq[1, 2]) / (contingency_table_any_SVs_chipseq[2, 1] / contingency_table_any_SVs_chipseq[1, 1])


    # any_somatic_SNVs
    contingency_table_any_somatic_SNVs <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_somatic_SNVs
    )

    cis_activation_frequency_WITHOUT_any_somatic_SNVs <- contingency_table_any_somatic_SNVs[2, 1] / (contingency_table_any_somatic_SNVs[1, 1] + contingency_table_any_somatic_SNVs[2, 1])
    cis_activation_frequency_WITH_any_somatic_SNVs <- contingency_table_any_somatic_SNVs[2, 2] / (contingency_table_any_somatic_SNVs[1, 2] + contingency_table_any_somatic_SNVs[2, 2])

    chi_square_test_any_somatic_SNVs <- chisq.test(contingency_table_any_somatic_SNVs)
    odds_ratio_any_somatic_SNVs <- (contingency_table_any_somatic_SNVs[2, 2] / contingency_table_any_somatic_SNVs[1, 2]) / (contingency_table_any_somatic_SNVs[2, 1] / contingency_table_any_somatic_SNVs[1, 1])

    # any_somatic_SNVs_genehancer
    contingency_table_any_somatic_SNVs_genehancer <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_somatic_SNVs_from_genehancer_TAD
    )

    cis_activation_frequency_WITHOUT_any_somatic_SNVs_genehancer <- contingency_table_any_somatic_SNVs_genehancer[2, 1] / (contingency_table_any_somatic_SNVs_genehancer[1, 1] + contingency_table_any_somatic_SNVs_genehancer[2, 1])
    cis_activation_frequency_WITH_any_somatic_SNVs_genehancer <- contingency_table_any_somatic_SNVs_genehancer[2, 2] / (contingency_table_any_somatic_SNVs_genehancer[1, 2] + contingency_table_any_somatic_SNVs_genehancer[2, 2])

    chi_square_test_any_somatic_SNVs_genehancer <- chisq.test(contingency_table_any_somatic_SNVs_genehancer)
    odds_ratio_any_somatic_SNVs_genehancer <- (contingency_table_any_somatic_SNVs_genehancer[2, 2] / contingency_table_any_somatic_SNVs_genehancer[1, 2]) / (contingency_table_any_somatic_SNVs_genehancer[2, 1] / contingency_table_any_somatic_SNVs_genehancer[1, 1])

    # any_somatic_SNVs_gene
    contingency_table_any_somatic_SNVs_gene <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_somatic_SNVs_from_gene_TAD
    )

    cis_activation_frequency_WITHOUT_any_somatic_SNVs_gene <- contingency_table_any_somatic_SNVs_gene[2, 1] / (contingency_table_any_somatic_SNVs_gene[1, 1] + contingency_table_any_somatic_SNVs_gene[2, 1])
    cis_activation_frequency_WITH_any_somatic_SNVs_gene <- contingency_table_any_somatic_SNVs_gene[2, 2] / (contingency_table_any_somatic_SNVs_gene[1, 2] + contingency_table_any_somatic_SNVs_gene[2, 2])

    chi_square_test_any_somatic_SNVs_gene <- chisq.test(contingency_table_any_somatic_SNVs_gene)
    odds_ratio_any_somatic_SNVs_gene <- (contingency_table_any_somatic_SNVs_gene[2, 2] / contingency_table_any_somatic_SNVs_gene[1, 2]) / (contingency_table_any_somatic_SNVs_gene[2, 1] / contingency_table_any_somatic_SNVs_gene[1, 1])

    # any_somatic_SNVs_chipseq
    contingency_table_any_somatic_SNVs_chipseq <- cont_table(
        processed_data$cis_activation_summary_table_ANY_mut_info$cis_activated_gene,
        processed_data$cis_activation_summary_table_ANY_mut_info$any_somatic_SNVs_with_chipseq
    )

    cis_activation_frequency_WITHOUT_any_somatic_SNVs_chipseq <- contingency_table_any_somatic_SNVs_chipseq[2, 1] / (contingency_table_any_somatic_SNVs_chipseq[1, 1] + contingency_table_any_somatic_SNVs_chipseq[2, 1])
    cis_activation_frequency_WITH_any_somatic_SNVs_chipseq <- contingency_table_any_somatic_SNVs_chipseq[2, 2] / (contingency_table_any_somatic_SNVs_chipseq[1, 2] + contingency_table_any_somatic_SNVs_chipseq[2, 2])

    chi_square_test_any_somatic_SNVs_chipseq <- chisq.test(contingency_table_any_somatic_SNVs_chipseq)
    odds_ratio_any_somatic_SNVs_chipseq <- (contingency_table_any_somatic_SNVs_chipseq[2, 2] / contingency_table_any_somatic_SNVs_chipseq[1, 2]) / (contingency_table_any_somatic_SNVs_chipseq[2, 1] / contingency_table_any_somatic_SNVs_chipseq[1, 1])



    if (has_run_tf_binding_site_analysis) {
        # any_relevant_somatic_SNVs
        contingency_table_any_relevant_somatic_SNVs <- cont_table(
            (!is.na(processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene)) & processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene,
            (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_gene_TAD > 0) |
                (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_chipseq > 0) |
                (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_genehancer > 0)
        )

        cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs <- contingency_table_any_relevant_somatic_SNVs[2, 1] / (contingency_table_any_relevant_somatic_SNVs[1, 1] + contingency_table_any_relevant_somatic_SNVs[2, 1])
        cis_activation_frequency_WITH_any_relevant_somatic_SNVs <- contingency_table_any_relevant_somatic_SNVs[2, 2] / (contingency_table_any_relevant_somatic_SNVs[1, 2] + contingency_table_any_relevant_somatic_SNVs[2, 2])

        chi_square_test_any_relevant_somatic_SNVs <- chisq.test(contingency_table_any_relevant_somatic_SNVs)
        odds_ratio_any_relevant_somatic_SNVs <- (contingency_table_any_relevant_somatic_SNVs[2, 2] / contingency_table_any_relevant_somatic_SNVs[1, 2]) / (contingency_table_any_relevant_somatic_SNVs[2, 1] / contingency_table_any_relevant_somatic_SNVs[1, 1])

        # any_relevant_somatic_SNVs_genehancer
        contingency_table_any_relevant_somatic_SNVs_genehancer <- cont_table(
            (!is.na(processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene)) & processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene,
            (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_genehancer > 0)
        )

        cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_genehancer <- contingency_table_any_relevant_somatic_SNVs_genehancer[2, 1] / (contingency_table_any_relevant_somatic_SNVs_genehancer[1, 1] + contingency_table_any_relevant_somatic_SNVs_genehancer[2, 1])
        cis_activation_frequency_WITH_any_relevant_somatic_SNVs_genehancer <- contingency_table_any_relevant_somatic_SNVs_genehancer[2, 2] / (contingency_table_any_relevant_somatic_SNVs_genehancer[1, 2] + contingency_table_any_relevant_somatic_SNVs_genehancer[2, 2])

        chi_square_test_any_relevant_somatic_SNVs_genehancer <- chisq.test(contingency_table_any_relevant_somatic_SNVs_genehancer)
        odds_ratio_any_relevant_somatic_SNVs_genehancer <- (contingency_table_any_relevant_somatic_SNVs_genehancer[2, 2] / contingency_table_any_relevant_somatic_SNVs_genehancer[1, 2]) / (contingency_table_any_relevant_somatic_SNVs_genehancer[2, 1] / contingency_table_any_relevant_somatic_SNVs_genehancer[1, 1])

        # any_relevant_somatic_SNVs_gene
        contingency_table_any_relevant_somatic_SNVs_gene <- cont_table(
            (!is.na(processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene)) & processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene,
            (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_gene_TAD > 0)
        )

        cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_gene <- contingency_table_any_relevant_somatic_SNVs_gene[2, 1] / (contingency_table_any_relevant_somatic_SNVs_gene[1, 1] + contingency_table_any_relevant_somatic_SNVs_gene[2, 1])
        cis_activation_frequency_WITH_any_relevant_somatic_SNVs_gene <- contingency_table_any_relevant_somatic_SNVs_gene[2, 2] / (contingency_table_any_relevant_somatic_SNVs_gene[1, 2] + contingency_table_any_relevant_somatic_SNVs_gene[2, 2])

        chi_square_test_any_relevant_somatic_SNVs_gene <- chisq.test(contingency_table_any_relevant_somatic_SNVs_gene)
        odds_ratio_any_relevant_somatic_SNVs_gene <- (contingency_table_any_relevant_somatic_SNVs_gene[2, 2] / contingency_table_any_relevant_somatic_SNVs_gene[1, 2]) / (contingency_table_any_relevant_somatic_SNVs_gene[2, 1] / contingency_table_any_relevant_somatic_SNVs_gene[1, 1])

        # any_relevant_somatic_SNVs_chipseq ###CHECK gene TAD?
        contingency_table_any_relevant_somatic_SNVs_chipseq <- cont_table(
            (!is.na(processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene)) & processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$cis_activated_gene,
            (processed_data$cis_activation_summary_with_n_relevant_somatic_SNVs_table$n_relevant_somatic_SNVs_via_chipseq > 0)
        )

        cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_chipseq <- contingency_table_any_relevant_somatic_SNVs_chipseq[2, 1] / (contingency_table_any_relevant_somatic_SNVs_chipseq[1, 1] + contingency_table_any_relevant_somatic_SNVs_chipseq[2, 1])
        cis_activation_frequency_WITH_any_relevant_somatic_SNVs_chipseq <- contingency_table_any_relevant_somatic_SNVs_chipseq[2, 2] / (contingency_table_any_relevant_somatic_SNVs_chipseq[1, 2] + contingency_table_any_relevant_somatic_SNVs_chipseq[2, 2])

        chi_square_test_any_relevant_somatic_SNVs_chipseq <- chisq.test(contingency_table_any_relevant_somatic_SNVs_chipseq)
        odds_ratio_any_relevant_somatic_SNVs_chipseq <- (contingency_table_any_relevant_somatic_SNVs_chipseq[2, 2] / contingency_table_any_relevant_somatic_SNVs_chipseq[1, 2]) / (contingency_table_any_relevant_somatic_SNVs_chipseq[2, 1] / contingency_table_any_relevant_somatic_SNVs_chipseq[1, 1])


        # ROC CURVES - gene
        class_rocit_gene <- (!is.na(processed_data$somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table$cis_activated_gene)) & processed_data$somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table$cis_activated_gene
        score_rocit_gene <- processed_data$somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table$score_diff

        if (length(unique(class_rocit_gene)) == 2) {
            # MAKE ROC PLOTS
            ROCit::rocit(
                score = score_rocit_gene,
                class = class_rocit_gene
            ) -> rocitdata_gene
            save_ROC_curve_plot(rocitdata_gene, file.path(HTML_report_figure_directory, paste0("TF_binding_ROC_via_gene.svg")))

            # MAKE SENS_SPEC PLOTS
            plot_sensitivity_and_specificity_curve(rocitdata_gene) %>%
                save_ggplot(file.path(HTML_report_figure_directory, paste0("TF_binding_SENS_SPEC_via_gene.png")))
        }


        # ROC CURVES - genehancer
        class_rocit_genehancer <- (!is.na(processed_data$somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table$cis_activated_gene)) & processed_data$somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table$cis_activated_gene
        score_rocit_genehancer <- processed_data$somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table$score_diff
        if (length(unique(class_rocit_genehancer)) == 2) {
            # MAKE ROC PLOTS
            ROCit::rocit(
                score = score_rocit_genehancer,
                class = class_rocit_genehancer
            ) -> rocitdata_genehancer
            save_ROC_curve_plot(rocitdata_genehancer, file.path(HTML_report_figure_directory, paste0("TF_binding_ROC_via_genehancer.svg")))
            # MAKE SENS_SPEC PLOTS
            plot_sensitivity_and_specificity_curve(rocitdata_genehancer) %>%
                save_ggplot(file.path(HTML_report_figure_directory, paste0("TF_binding_SENS_SPEC_via_genehancer.png")))
        }

        # ROC CURVES - chipseq
        class_rocit_chipseq <- (!is.na(processed_data$somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table$cis_activated_gene)) & processed_data$somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table$cis_activated_gene
        score_rocit_chipseq <- processed_data$somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table$score_diff
        if (length(unique(class_rocit_chipseq)) == 2) {
            # MAKE ROC PLOTS
            ROCit::rocit(
                score = score_rocit_chipseq,
                class = class_rocit_chipseq
            ) -> rocitdata_chipseq
            save_ROC_curve_plot(rocitdata_chipseq, file.path(HTML_report_figure_directory, paste0("TF_binding_ROC_via_chipseq.svg")))
            # MAKE SENS_SPEC PLOTS
            plot_sensitivity_and_specificity_curve(rocitdata_chipseq) %>%
                save_ggplot(file.path(HTML_report_figure_directory, paste0("TF_binding_SENS_SPEC_via_chipseq.png")))
        }
    }


    # SV subcategories - genehancer
    cat("create SV subcategory - genehancer table...\n")
    SV_genehancer_position_category_summary_cis_activated <- processed_data$SV_genehancer_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::group_by(position_category_SV) %>%
        dplyr::summarize(n_cis_activated = dplyr::n(), n_unique_sample_gene_combinations_cis_activated = dplyr::n_distinct(sample_ID, connected_gene))

    SV_genehancer_position_category_summary_all <- processed_data$SV_genehancer_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::group_by(position_category_SV) %>%
        dplyr::summarize(n = dplyr::n(), n_unique_sample_gene_combinations = dplyr::n_distinct(sample_ID, connected_gene))

    SV_genehancer_position_category_summary_sorted_freq <- dplyr::left_join(SV_genehancer_position_category_summary_all, SV_genehancer_position_category_summary_cis_activated, by = c("position_category_SV")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(freq))

    SV_genehancer_position_category_summary_sorted_n <- dplyr::left_join(SV_genehancer_position_category_summary_all, SV_genehancer_position_category_summary_cis_activated, by = c("position_category_SV")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(n))

    SV_genehancer_position_category_summary_HTML_sorted_freq <- knitr::kable(SV_genehancer_position_category_summary_sorted_freq, format = "html", table.attr = "class=\"table sv-pos-table table-striped table-bordered\"")
    SV_genehancer_position_category_summary_HTML_sorted_n <- knitr::kable(SV_genehancer_position_category_summary_sorted_n, format = "html", table.attr = "class=\"table sv-pos-table table-striped table-bordered\"")
    cat("created SV subcategory - genehancer table!\n")

    # SV subcategories - chipseq
    cat("create SV subcategory - chipseq table...\n")
    SV_chipseq_position_category_summary_cis_activated <- processed_data$SV_chipseq_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::group_by(position_category_SV) %>%
        dplyr::summarize(n_cis_activated = dplyr::n(), n_unique_sample_gene_combinations_cis_activated = dplyr::n_distinct(sample_ID, gene_name))

    SV_chipseq_position_category_summary_all <- processed_data$SV_chipseq_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::group_by(position_category_SV) %>%
        dplyr::summarize(n = dplyr::n(), n_unique_sample_gene_combinations = dplyr::n_distinct(sample_ID, gene_name))

    SV_chipseq_position_category_summary_sorted_freq <- dplyr::left_join(SV_chipseq_position_category_summary_all, SV_chipseq_position_category_summary_cis_activated, by = c("position_category_SV")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(freq))

    SV_chipseq_position_category_summary_sorted_n <- dplyr::left_join(SV_chipseq_position_category_summary_all, SV_chipseq_position_category_summary_cis_activated, by = c("position_category_SV")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(n))

    SV_chipseq_position_category_summary_HTML_sorted_freq <- knitr::kable(SV_chipseq_position_category_summary_sorted_freq, format = "html", table.attr = "class=\"table sv-pos-table table-striped table-bordered\"")
    SV_chipseq_position_category_summary_HTML_sorted_n <- knitr::kable(SV_chipseq_position_category_summary_sorted_n, format = "html", table.attr = "class=\"table sv-pos-table table-striped table-bordered\"")
    cat("created SV subcategory - chipseq table!\n")

    # CNA subcategories - genehancer
    cat("create CNA subcategory - genehancer table...\n")
    CNA_genehancer_position_category_summary_cis_activated <- processed_data$CNA_genehancer_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::group_by(position_category_CNA) %>%
        dplyr::summarize(n_cis_activated = dplyr::n(), n_unique_sample_gene_combinations_cis_activated = dplyr::n_distinct(sample_ID, connected_gene))

    CNA_genehancer_position_category_summary_all <- processed_data$CNA_genehancer_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::group_by(position_category_CNA) %>%
        dplyr::summarize(n = dplyr::n(), n_unique_sample_gene_combinations = dplyr::n_distinct(sample_ID, connected_gene))

    CNA_genehancer_position_category_summary_sorted_freq <- dplyr::left_join(CNA_genehancer_position_category_summary_all, CNA_genehancer_position_category_summary_cis_activated, by = c("position_category_CNA")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(freq))

    CNA_genehancer_position_category_summary_sorted_n <- dplyr::left_join(CNA_genehancer_position_category_summary_all, CNA_genehancer_position_category_summary_cis_activated, by = c("position_category_CNA")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(n))

    CNA_genehancer_position_category_summary_HTML_sorted_freq <- knitr::kable(CNA_genehancer_position_category_summary_sorted_freq, format = "html", table.attr = "class=\"table cna-pos-table table-striped table-bordered\"")
    CNA_genehancer_position_category_summary_HTML_sorted_n <- knitr::kable(CNA_genehancer_position_category_summary_sorted_n, format = "html", table.attr = "class=\"table cna-pos-table table-striped table-bordered\"")
    cat("created CNA subcategory - genehancer table!\n")

    # CNA subcategories - chipseq
    cat("create CNA subcategory - chipseq table...\n")
    CNA_chipseq_position_category_summary_cis_activated <- processed_data$CNA_chipseq_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::group_by(position_category_CNA) %>%
        dplyr::summarize(n_cis_activated = dplyr::n(), n_unique_sample_gene_combinations_cis_activated = dplyr::n_distinct(sample_ID, gene_name))

    CNA_chipseq_position_category_summary_all <- processed_data$CNA_chipseq_all_table %>%
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!is.na(cis_activated_gene)))) %>%
        dplyr::group_by(position_category_CNA) %>%
        dplyr::summarize(n = dplyr::n(), n_unique_sample_gene_combinations = dplyr::n_distinct(sample_ID, gene_name))

    CNA_chipseq_position_category_summary_sorted_freq <- dplyr::left_join(CNA_chipseq_position_category_summary_all, CNA_chipseq_position_category_summary_cis_activated, by = c("position_category_CNA")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(freq))

    CNA_chipseq_position_category_summary_sorted_n <- dplyr::left_join(CNA_chipseq_position_category_summary_all, CNA_chipseq_position_category_summary_cis_activated, by = c("position_category_CNA")) %>%
        tidyr::replace_na(list(n_cis_activated = 0, n_unique_sample_gene_combinations_cis_activated = 0)) %>%
        dplyr::mutate(freq = (n_cis_activated / n), freq_unique_sample_gene_combinations = (n_unique_sample_gene_combinations_cis_activated / n_unique_sample_gene_combinations)) %>%
        dplyr::arrange(desc(n))

    CNA_chipseq_position_category_summary_HTML_sorted_freq <- knitr::kable(CNA_chipseq_position_category_summary_sorted_freq, format = "html", table.attr = "class=\"table cna-pos-table table-striped table-bordered\"")
    CNA_chipseq_position_category_summary_HTML_sorted_n <- knitr::kable(CNA_chipseq_position_category_summary_sorted_n, format = "html", table.attr = "class=\"table cna-pos-table table-striped table-bordered\"")
    cat("created CNA subcategory - chipseq table!\n")


    # util function
    chi_sq_test_as_HTML <- function(chi_sq_test) {
        return(paste0(
            "<p>Pearson's Chi-squared test with Yates' continuity correction</p>
            <p>X-squared = ", chi_sq_test$statistic, ", df = ", chi_sq_test$parameter, ", p-value = ", chi_sq_test$p.value, "</p>"
        ))
    }

    HTML <- paste0('
        <div class = "card my-4">
            <h5 class="card-header">Relation between mutations and cis activation
                <span class="badge bg-primary float-end">Mutation Type</span>
            </h5>
            <div class="p-3">
                <h3>All</h3>
                <h5>Overview</h5>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no mutation </th>
                        <th style="text-align:right;"> any mutation </th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_muts[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_muts[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_muts, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_muts, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_muts, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_muts), '
                <h3>Genehancer</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no mutation (via genehancer)</th>
                        <th style="text-align:right;"> any mutation (via genehancer)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_muts_genehancer[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_genehancer[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_muts_genehancer[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_genehancer[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_muts_genehancer, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_muts_genehancer, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_muts_genehancer, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_muts_genehancer), '
                <h3>Gene</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no mutation (via gene TAD)</th>
                        <th style="text-align:right;"> any mutation (via gene TAD)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_muts_gene[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_gene[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_muts_gene[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_gene[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_muts_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_muts_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_muts_gene, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_muts_gene), '
                <h3>ChIP-Seq</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no mutation (with affected ChIP-Seq region)</th>
                        <th style="text-align:right;"> any mutation (with affected ChIP-Seq region)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_muts_chipseq[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_chipseq[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_muts_chipseq[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_muts_chipseq[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_muts_chipseq, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_muts_chipseq, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_muts_chipseq, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_muts_chipseq), '
            </div>
        </div>

        <div class = "card my-4">
            <h5 class="card-header">Relation between SVs and cis activation
                <span class="badge bg-primary float-end">Mutation Type</span>
            </h5>
            <div class="p-3">
                <h3>All</h3>
                <h5>Overview</h5>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no SV </th>
                        <th style="text-align:right;"> any SV </th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_SVs[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_SVs[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_SVs, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_SVs, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_SVs, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_SVs), '
                <h3>Genehancer</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no SV (via genehancer)</th>
                        <th style="text-align:right;"> any SV (via genehancer)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_SVs_genehancer[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_genehancer[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_genehancer[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_genehancer[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_SVs_genehancer, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_SVs_genehancer, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_SVs_genehancer, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_SVs_genehancer), '
                <h3>Gene</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no SV (via gene TAD)</th>
                        <th style="text-align:right;"> any SV (via gene TAD)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_SVs_gene[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_gene[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_gene[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_gene[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_SVs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_SVs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_SVs_gene, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_SVs_gene), '
                <h3>ChIP-Seq</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no SV (with affected ChIP-Seq region)</th>
                        <th style="text-align:right;"> any SV (with affected ChIP-Seq region)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_SVs_chipseq[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_chipseq[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_chipseq[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_SVs_chipseq[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_SVs_chipseq, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_SVs_chipseq, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_SVs_chipseq, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_SVs_chipseq), '
                <!-- EXCLUDED AS FOR RIGHT NOW
                <h3>SV position categories (relative to gene and genehancer region)</h3>
                <h5>sorted by descending n</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', SV_genehancer_position_category_summary_HTML_sorted_n, '
                </div>
                <h5 class="mt-3">sorted by descending freq</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', SV_genehancer_position_category_summary_HTML_sorted_freq, '
                </div>
                <h3>SV position categories (relative to gene and chipseq region)</h3>
                <h5>sorted by descending n</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', SV_chipseq_position_category_summary_HTML_sorted_n, '
                </div>
                <h5 class="mt-3">sorted by descending freq</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', SV_chipseq_position_category_summary_HTML_sorted_freq, '
                </div>
                -->
            </div>
        </div>

        <div class = "card my-4">
            <h5 class="card-header">Relation between CNAs and cis activation
                <span class="badge bg-primary float-end">Mutation Type</span>
            </h5>
            <div class="p-3">
                <h3>All</h3>
                <h5>Overview</h5>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no CNA </th>
                        <th style="text-align:right;"> any CNA </th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_CNAs[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_CNAs, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_CNAs, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_CNAs, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_CNAs), '
                <h3>Genehancer</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no CNA (via genehancer)</th>
                        <th style="text-align:right;"> any CNA (via genehancer)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_genehancer[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_genehancer[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_genehancer[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_genehancer[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_CNAs_genehancer, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_CNAs_genehancer, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_CNAs_genehancer, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_CNAs_genehancer), '
                <h3>Gene</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no CNA (via gene TAD)</th>
                        <th style="text-align:right;"> any CNA (via gene TAD)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_gene[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_gene[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_gene[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_gene[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_CNAs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_CNAs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_CNAs_gene, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_CNAs_gene), '
                <h3>ChIP-Seq</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no CNA (with affected ChIP-Seq region)</th>
                        <th style="text-align:right;"> any CNA (with affected ChIP-Seq region)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_chipseq[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_chipseq[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_chipseq[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_CNAs_chipseq[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_CNAs_chipseq, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_CNAs_chipseq, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_CNAs_chipseq, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_CNAs_chipseq), '
                <!-- EXCLUDED AS FOR RIGHT NOW

                <h3>CNA position categories (relative to gene and genehancer region)</h3>
                <h5>sorted by descending n</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', CNA_genehancer_position_category_summary_HTML_sorted_n, '
                </div>
                <h5 class="mt-3">sorted by descending freq</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', CNA_genehancer_position_category_summary_HTML_sorted_freq, '
                </div>
                <h3>CNA position categories (relative to gene and chipseq region)</h3>
                <h5>sorted by descending n</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', CNA_chipseq_position_category_summary_HTML_sorted_n, '
                </div>
                <h5 class="mt-3">sorted by descending freq</h5>
                <div style="display: block;
                    overflow-y: scroll;
                    height: 30em;">
                ', CNA_chipseq_position_category_summary_HTML_sorted_freq, '
                </div>
            -->
            </div>
        </div>

        <div class = "card my-4">
            <h5 class="card-header">Relation between somatic SNVs and cis activation
                <span class="badge bg-primary float-end">Mutation Type</span>
            </h5>
            <div class="p-3">
                <h3>All</h3>
                <h5>Overview</h5>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no somatic_SNV </th>
                        <th style="text-align:right;"> any somatic_SNV </th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_somatic_SNVs, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_somatic_SNVs, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_somatic_SNVs, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_somatic_SNVs), '
                <h3>Genehancer</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no somatic_SNV (via genehancer)</th>
                        <th style="text-align:right;"> any somatic_SNV (via genehancer)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_genehancer[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_genehancer[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_genehancer[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_genehancer[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_somatic_SNVs_genehancer, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_somatic_SNVs_genehancer, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_somatic_SNVs_genehancer, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_somatic_SNVs_genehancer), '
                <h3>Gene</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no somatic_SNV (via gene TAD)</th>
                        <th style="text-align:right;"> any somatic_SNV (via gene TAD)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_gene[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_gene[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_gene[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_gene[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_somatic_SNVs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_somatic_SNVs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_somatic_SNVs_gene, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_somatic_SNVs_gene), '
                <h3>ChIP-Seq</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no somatic_SNV (with affected ChIP-Seq region)</th>
                        <th style="text-align:right;"> any somatic_SNV (with affected ChIP-Seq region)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_chipseq[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_chipseq[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_chipseq[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_somatic_SNVs_chipseq[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_somatic_SNVs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_somatic_SNVs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_somatic_SNVs_chipseq, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_somatic_SNVs_chipseq), "
            </div>
        </div>")


    if (has_run_tf_binding_site_analysis) {
        HTML <- paste0(
            HTML,
            '<div class = "card my-4">
            <h5 class="card-header">Relation between \"relevant\" somatic SNVs and cis activation
                <span class="badge bg-primary float-end">Mutation Type</span>
            </h5>
            <div class="p-3 border-bottom">
              \"Relevant\" somatic SNVs/InDels in this context describe SNVs, that introduce new significant transcription factor binding sites to transcription factors with a minimum gene expression of > 5 FPKM in the same tumor sample.
            </div>
            <div class="p-3">
                <h3>All</h3>
                <h5>Overview</h5>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no relevant somatic SNV </th>
                        <th style="text-align:right;"> any relevant somatic SNV </th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_relevant_somatic_SNVs, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_relevant_somatic_SNVs, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_relevant_somatic_SNVs), '
                <h3>Genehancer</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no relevant somatic SNV (via genehancer)</th>
                        <th style="text-align:right;"> any relevant somatic SNV (via genehancer)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_genehancer[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_genehancer[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_genehancer[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_genehancer[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_genehancer, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_relevant_somatic_SNVs_genehancer, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_relevant_somatic_SNVs_genehancer, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_relevant_somatic_SNVs_genehancer), '
                <h5>TF binding site analysis ROC</h5>
                <div class="d-flex">
                <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_ROC_via_genehancer.svg"></div>
                <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_SENS_SPEC_via_genehancer.png"></div>
                </div>
                <h3>Gene</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no relevant somatic SNV (via gene TAD)</th>
                        <th style="text-align:right;"> any relevant somatic SNV (via gene TAD)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_gene[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_gene[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_gene[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_gene[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_relevant_somatic_SNVs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_relevant_somatic_SNVs_gene, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_relevant_somatic_SNVs_gene), '
                <h5>TF binding site analysis ROC</h5>
                <div class="d-flex">
                <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_ROC_via_gene.svg"></div>
                <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_SENS_SPEC_via_gene.png"></div>
                </div>
                <h3>ChIP-Seq</h3>
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                        <th style="text-align:left;">   </th>
                        <th style="text-align:right;"> no relevant somatic SNV (with affected ChIP-Seq region)</th>
                        <th style="text-align:right;"> any relevant somatic SNV (with affected ChIP-Seq region)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align:left;"> gene NOT cis activated </td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_chipseq[1, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_chipseq[1, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> gene cis activated</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_chipseq[2, 1], '</td>
                            <td style="text-align:right;">', contingency_table_any_relevant_somatic_SNVs_chipseq[2, 2], '</td>
                        </tr>
                        <tr>
                            <td style="text-align:left;"> cis activation frequency</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITHOUT_any_relevant_somatic_SNVs_gene, '</td>
                            <td style="text-align:right;">', cis_activation_frequency_WITH_any_relevant_somatic_SNVs_gene, "</td>
                        </tr>
                    </tbody>
                </table>
                <p>ODDS RATIO: ", odds_ratio_any_relevant_somatic_SNVs_chipseq, "</p>
                <h5>Chi-Square Test</h5>
                ", chi_sq_test_as_HTML(chi_square_test_any_relevant_somatic_SNVs_chipseq), '
                <h5>TF binding site analysis ROC</h5>
                <div class="d-flex">
                    <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_ROC_via_chipseq.svg"></div>
                    <div class="p-3 w-50"><img class="w-100" src="figures/TF_binding_SENS_SPEC_via_chipseq.png"></div>
                </div>
            </div>
        </div>'
        )
    }

    return(HTML)
}