generate_alpha_scale_label_for_cis_activated_samples_per_gene_bar_plot <- function(matched_TAD, mutation_type) {
    dplyr::case_when(
        (matched_TAD == "anywhere") & (mutation_type == "all") ~ "Mutations (related)",
        (matched_TAD == "anywhere") & (mutation_type == "SV_only") ~ "SVs (related)",
        (matched_TAD == "anywhere") & (mutation_type == "CNA_only") ~ "CNAs (related)",
        (matched_TAD == "anywhere") & (mutation_type == "somatic_SNV_only") ~ "somatic SNVs (related)",
        (matched_TAD == "anywhere") & (mutation_type == "relevant_somatic_SNV_only") ~ "somatic SNVs with relevant TF binding site (related)",
        (matched_TAD == "gene") & (mutation_type) == "all" ~ "Mutations in gene TAD",
        (matched_TAD == "gene") & (mutation_type) == "SV_only" ~ "SVs in gene TAD",
        (matched_TAD == "gene") & (mutation_type) == "CNA_only" ~ "CNAs in gene TAD",
        (matched_TAD == "gene") & (mutation_type) == "somatic_SNV_only" ~ "somatic SNVs in gene TAD",
        (matched_TAD == "gene") & (mutation_type) == "relevant_somatic_SNV_only" ~ "somatic SNVs with relevant TF binding site in gene TAD",
        (matched_TAD == "genehancer") & (mutation_type) == "all" ~ "Mutations in genehancer TADs",
        (matched_TAD == "genehancer") & (mutation_type) == "SV_only" ~ "SVs in genehancer TADs",
        (matched_TAD == "genehancer") & (mutation_type) == "CNA_only" ~ "CNAs in genehancer TADs",
        (matched_TAD == "genehancer") & (mutation_type) == "somatic_SNV_only" ~ "somatic SNVs in genehancer TADs",
        (matched_TAD == "genehancer") & (mutation_type) == "relevant_somatic_SNV_only" ~ "somatic SNVs with relevant TF binding site in genehancer TADs",
        (matched_TAD == "chipseq") & (mutation_type) == "all" ~ "Mutations with affected ChIP-Seq regions",
        (matched_TAD == "chipseq") & (mutation_type) == "SV_only" ~ "SVs with affected ChIP-Seq regions",
        (matched_TAD == "chipseq") & (mutation_type) == "CNA_only" ~ "CNAs with affected ChIP-Seq regions",
        (matched_TAD == "chipseq") & (mutation_type) == "somatic_SNV_only" ~ "somatic SNVs with affected ChIP-Seq regions",
        (matched_TAD == "chipseq") & (mutation_type) == "relevant_somatic_SNV_only" ~ "somatic SNVs with relevant TF binding site with affected ChIP-Seq regions",
        TRUE ~ " - "
    )
}

generate_modified_cis_activation_summary_table_mut_info_for_cis_activated_samples_per_gene_bar_plot <- function(processed_data, copy_number_normalized = FALSE, matched_TAD, mutation_type, filter_occurrence, filter_mutations, sort_by, pc_only) {
    if (copy_number_normalized == TRUE) {
        cis_act_df <- processed_data$cis_activation_summary_table_mut_info %>%
            dplyr::filter(cis_activated_gene_copy_number_normalized == TRUE)
    } else {
        cis_act_df <- processed_data$cis_activation_summary_table_mut_info %>%
            dplyr::filter(cis_activated_gene == TRUE)
    }

    temp_df <- cis_act_df %>%
        dplyr::mutate(n_muts_of_mut_type = dplyr::case_when(
            (matched_TAD == "anywhere") & (mutation_type == "all") ~ n_muts,
            (matched_TAD == "anywhere") & (mutation_type == "SV_only") ~ n_SVs,
            (matched_TAD == "anywhere") & (mutation_type == "CNA_only") ~ n_CNAs,
            (matched_TAD == "anywhere") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs,
            (matched_TAD == "anywhere") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_valid_for_cis_activated_genes,
            (matched_TAD == "gene") & (mutation_type == "all") ~ n_muts_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "SV_only") ~ n_SVs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "CNA_only") ~ n_CNAs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes,
            (matched_TAD == "genehancer") & (mutation_type == "all") ~ n_muts_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "SV_only") ~ n_SVs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "CNA_only") ~ n_CNAs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes,
            (matched_TAD == "chipseq") & (mutation_type == "all") ~ n_muts_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "SV_only") ~ n_SVs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "CNA_only") ~ n_CNAs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "relevant_somatic_SNV_only") ~ 0,
            TRUE ~ n_muts
        )) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(n_cis_activated_samples_across_cohort = dplyr::n()) %>%
        dplyr::mutate(n_cis_activated_samples_with_mut_across_cohort = sum(n_muts_of_mut_type > 0)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n_cis_activated_samples_across_cohort >= filter_occurrence) %>%
        dplyr::filter(n_cis_activated_samples_with_mut_across_cohort >= filter_mutations) %>%
        dplyr::filter(gene_type == "protein_coding" | (!pc_only))

    if (sort_by == "occurrence") {
        temp_df_sorted_factors <- temp_df %>%
            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_across_cohort))
    }
    if (sort_by == "mutations") {
        temp_df_sorted_factors <- temp_df %>%
            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_with_mut_across_cohort))
    }

    return(temp_df_sorted_factors)
}

plot_cis_activated_samples_per_gene_bar_plot_integrated <- function(modified_cis_activation_summary_table_mut_info, color_palette, alpha_scale_name = "Mutations (related)", axis.text.y.size = 1.7) {
    # a.k. toggle plot
    ggplot2::ggplot(modified_cis_activation_summary_table_mut_info) +
        ggplot2::geom_bar(ggplot2::aes(y = gene_name, fill = subgroup, alpha = (n_muts_of_mut_type > 0))) +
        color_palette$scale_fill_subgroup +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = axis.text.y.size)) +
        ggplot2::scale_alpha_manual(values = c(1, 0.5), breaks = c(TRUE, FALSE), labels = c("mutations detected", "no mutations detected"), name = alpha_scale_name)
}

plot_cis_activated_samples_per_gene_bar_plot <- function(processed_data,
                                                         color_palette,
                                                         copy_number_normalized = FALSE,
                                                         mutation_type = "all",
                                                         matched_TAD = "anywhere",
                                                         sort_by = "occurrence",
                                                         filter_occurrence = 1,
                                                         filter_mutations = 1,
                                                         pc_only = false,
                                                         blacklist_genes = character(0),
                                                         axis.text.y.size = 1.7) {
    # for manual plot generation
    if (copy_number_normalized == TRUE) {
        temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
            dplyr::filter(cis_activated_gene_copy_number_normalized == TRUE)
    } else {
        temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
            dplyr::filter(cis_activated_gene == TRUE)
    }

    temp_df2 <- temp_df1 %>%
        dplyr::mutate(n_muts_of_mut_type = dplyr::case_when(
            (matched_TAD == "anywhere") & (mutation_type == "all") ~ n_muts,
            (matched_TAD == "anywhere") & (mutation_type == "SV_only") ~ n_SVs,
            (matched_TAD == "anywhere") & (mutation_type == "CNA_only") ~ n_CNAs,
            (matched_TAD == "anywhere") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs,
            (matched_TAD == "anywhere") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_valid_for_cis_activated_genes,
            (matched_TAD == "gene") & (mutation_type == "all") ~ n_muts_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "SV_only") ~ n_SVs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "CNA_only") ~ n_CNAs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_gene_TAD,
            (matched_TAD == "gene") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes,
            (matched_TAD == "genehancer") & (mutation_type == "all") ~ n_muts_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "SV_only") ~ n_SVs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "CNA_only") ~ n_CNAs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_genehancer_TAD,
            (matched_TAD == "genehancer") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes,
            (matched_TAD == "chipseq") & (mutation_type == "all") ~ n_muts_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "SV_only") ~ n_SVs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "CNA_only") ~ n_CNAs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_with_chipseq,
            (matched_TAD == "chipseq") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes,
            TRUE ~ n_muts
        )) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(n_cis_activated_samples_across_cohort = dplyr::n()) %>%
        dplyr::mutate(n_cis_activated_samples_with_mut_across_cohort = sum(n_muts_of_mut_type > 0)) %>%
        dplyr::ungroup()

    alpha_scale_name <-
        generate_alpha_scale_label_for_cis_activated_samples_per_gene_bar_plot(matched_TAD, mutation_type)

    if (sort_by == "occurrence") {
        temp_df3 <- temp_df2 %>%
            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_across_cohort))
    }
    if (sort_by == "mutations") {
        temp_df3 <- temp_df2 %>%
            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_with_mut_across_cohort))
    }
    if (!(sort_by %in% c("occurrence", "mutations"))) {
        cat("WRONG sort_by parameter\n")
    }

    temp_df4 <- temp_df3 %>%
        dplyr::filter(n_cis_activated_samples_across_cohort >= filter_occurrence)

    temp_df5 <- temp_df4 %>%
        dplyr::filter(n_cis_activated_samples_with_mut_across_cohort >= filter_mutations)

    final_df <- temp_df5 %>%
        dplyr::filter(gene_type == "protein_coding" | (!pc_only)) %>%
        dplyr::filter(!(gene_name %in% blacklist_genes))

    return(plot_cis_activated_samples_per_gene_bar_plot_integrated(final_df, color_palette, alpha_scale_name, axis.text.y.size = axis.text.y.size))
}

create_all_cis_activated_samples_per_gene_bar_plots <- function(processed_data, color_palette, HTML_report_figure_directory) {
    cat("Creating togglable plots...\n")
    progressr::with_progress({
        p <- progressr::progressor(steps = 960)

        # a.k. toggle plot
        for (copy_number_normalized in c(FALSE, TRUE)) {
            if (copy_number_normalized == TRUE) {
                temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
                    dplyr::filter(cis_activated_gene_copy_number_normalized == TRUE)
            } else {
                temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
                    dplyr::filter(cis_activated_gene == TRUE)
            }

            for (mutation_type in c("all", "SV_only", "CNA_only", "somatic_SNV_only", "relevant_somatic_SNV_only")) {
                for (matched_TAD in c("anywhere", "gene", "genehancer", "chipseq")) {
                    temp_df2 <- temp_df1 %>%
                        dplyr::mutate(n_muts_of_mut_type = dplyr::case_when(
                            (matched_TAD == "anywhere") & (mutation_type == "all") ~ n_muts,
                            (matched_TAD == "anywhere") & (mutation_type == "SV_only") ~ n_SVs,
                            (matched_TAD == "anywhere") & (mutation_type == "CNA_only") ~ n_CNAs,
                            (matched_TAD == "anywhere") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs,
                            (matched_TAD == "anywhere") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_valid_for_cis_activated_genes,
                            (matched_TAD == "gene") & (mutation_type == "all") ~ n_muts_from_gene_TAD,
                            (matched_TAD == "gene") & (mutation_type == "SV_only") ~ n_SVs_from_gene_TAD,
                            (matched_TAD == "gene") & (mutation_type == "CNA_only") ~ n_CNAs_from_gene_TAD,
                            (matched_TAD == "gene") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_gene_TAD,
                            (matched_TAD == "gene") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes,
                            (matched_TAD == "genehancer") & (mutation_type == "all") ~ n_muts_from_genehancer_TAD,
                            (matched_TAD == "genehancer") & (mutation_type == "SV_only") ~ n_SVs_from_genehancer_TAD,
                            (matched_TAD == "genehancer") & (mutation_type == "CNA_only") ~ n_CNAs_from_genehancer_TAD,
                            (matched_TAD == "genehancer") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_genehancer_TAD,
                            (matched_TAD == "genehancer") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes,
                            (matched_TAD == "chipseq") & (mutation_type == "all") ~ n_muts_with_chipseq,
                            (matched_TAD == "chipseq") & (mutation_type == "SV_only") ~ n_SVs_with_chipseq,
                            (matched_TAD == "chipseq") & (mutation_type == "CNA_only") ~ n_CNAs_with_chipseq,
                            (matched_TAD == "chipseq") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_with_chipseq,
                            (matched_TAD == "chipseq") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes,
                            TRUE ~ n_muts
                        )) %>%
                        dplyr::group_by(gene_name) %>%
                        dplyr::mutate(n_cis_activated_samples_across_cohort = dplyr::n()) %>%
                        dplyr::mutate(n_cis_activated_samples_with_mut_across_cohort = sum(n_muts_of_mut_type > 0)) %>%
                        dplyr::ungroup()

                    alpha_scale_name <-
                        generate_alpha_scale_label_for_cis_activated_samples_per_gene_bar_plot(matched_TAD, mutation_type)

                    for (sort_by in c("occurrence", "mutations")) {
                        if (sort_by == "occurrence") {
                            temp_df3 <- temp_df2 %>%
                                dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_across_cohort))
                        }
                        if (sort_by == "mutations") {
                            temp_df3 <- temp_df2 %>%
                                dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_with_mut_across_cohort))
                        }

                        for (filter_occurrence in c(1, 2, 3)) {
                            temp_df4 <- temp_df3 %>%
                                dplyr::filter(n_cis_activated_samples_across_cohort >= filter_occurrence)

                            for (filter_mutations in c(0, 1)) {
                                temp_df5 <- temp_df4 %>%
                                    dplyr::filter(n_cis_activated_samples_with_mut_across_cohort >= filter_mutations)



                                for (pc_only in c(TRUE, FALSE)) {
                                    final_df <- temp_df5 %>%
                                        dplyr::filter(gene_type == "protein_coding" | (!pc_only))

                                    cis_activated_samples_per_gene_across_cohort_bar_plot_path <- file.path(HTML_report_figure_directory, paste0("cis_activated_samples", ifelse(copy_number_normalized, "__copy_number_normalized_", ""), "_per_gene_across_cohort_bar_plot___matched_TAD_", matched_TAD, "__mutation_type_", mutation_type, "__filter_occurrence_", filter_occurrence, "__filter_mutations_", filter_mutations, "__sort_", sort_by, ifelse(pc_only, "__pc_only", ""), ".svg"))

                                    plot_cis_activated_samples_per_gene_bar_plot_integrated(final_df, color_palette, alpha_scale_name) %>%
                                        save_ggplot(cis_activated_samples_per_gene_across_cohort_bar_plot_path)

                                    p()
                                }
                            }
                        }
                    }
                }
            }
        }
    })
}

create_all_cis_activated_samples_per_gene_bar_plots_PARALLELIZATION <- function(processed_data, color_palette, HTML_report_figure_directory) {
    # a.k. toggle plot
    future::plan(future::multisession)

    cat("Creating togglable plots... (on ")
    cat(future::nbrOfWorkers())
    cat(" workers)\n")

    progressr::with_progress({
        p <- progressr::progressor(steps = 960)
    
        for (copy_number_normalized in c(FALSE, TRUE)) {
            if (copy_number_normalized == TRUE) {
                temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
                    dplyr::filter(cis_activated_gene_copy_number_normalized == TRUE)
            } else {
                temp_df1 <- processed_data$cis_activation_summary_table_mut_info %>%
                    dplyr::filter(cis_activated_gene == TRUE)
            }

            # OPTIMAL PLACE TO INITIALIZE PARALLELIZATION

            mutation_type_array <- character(20)
            matched_TAD_array <- character(20)
            iterators <- numeric(20)

            iterator <- 1

            for (mutation_type in c("all", "SV_only", "CNA_only", "somatic_SNV_only", "relevant_somatic_SNV_only")) {
                for (matched_TAD in c("anywhere", "gene", "genehancer", "chipseq")) {
                    mutation_type_array[iterator] <- mutation_type
                    matched_TAD_array[iterator] <- matched_TAD
                    iterators[iterator] <- iterator
                    iterator <- iterator + 1
                }
            }

            furrr::future_walk(iterators, function(i) {
                mutation_type <- mutation_type_array[i]
                matched_TAD <- matched_TAD_array[i]

                temp_df2 <- temp_df1 %>%
                    dplyr::mutate(n_muts_of_mut_type = dplyr::case_when(
                        (matched_TAD == "anywhere") & (mutation_type == "all") ~ n_muts,
                        (matched_TAD == "anywhere") & (mutation_type == "SV_only") ~ n_SVs,
                        (matched_TAD == "anywhere") & (mutation_type == "CNA_only") ~ n_CNAs,
                        (matched_TAD == "anywhere") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs,
                        (matched_TAD == "anywhere") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_valid_for_cis_activated_genes,
                        (matched_TAD == "gene") & (mutation_type == "all") ~ n_muts_from_gene_TAD,
                        (matched_TAD == "gene") & (mutation_type == "SV_only") ~ n_SVs_from_gene_TAD,
                        (matched_TAD == "gene") & (mutation_type == "CNA_only") ~ n_CNAs_from_gene_TAD,
                        (matched_TAD == "gene") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_gene_TAD,
                        (matched_TAD == "gene") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes,
                        (matched_TAD == "genehancer") & (mutation_type == "all") ~ n_muts_from_genehancer_TAD,
                        (matched_TAD == "genehancer") & (mutation_type == "SV_only") ~ n_SVs_from_genehancer_TAD,
                        (matched_TAD == "genehancer") & (mutation_type == "CNA_only") ~ n_CNAs_from_genehancer_TAD,
                        (matched_TAD == "genehancer") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_from_genehancer_TAD,
                        (matched_TAD == "genehancer") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes,
                        (matched_TAD == "chipseq") & (mutation_type == "all") ~ n_muts_with_chipseq,
                        (matched_TAD == "chipseq") & (mutation_type == "SV_only") ~ n_SVs_with_chipseq,
                        (matched_TAD == "chipseq") & (mutation_type == "CNA_only") ~ n_CNAs_with_chipseq,
                        (matched_TAD == "chipseq") & (mutation_type == "somatic_SNV_only") ~ n_somatic_SNVs_with_chipseq,
                        (matched_TAD == "chipseq") & (mutation_type == "relevant_somatic_SNV_only") ~ n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes,
                        TRUE ~ n_muts
                    )) %>%
                    dplyr::group_by(gene_name) %>%
                    dplyr::mutate(n_cis_activated_samples_across_cohort = dplyr::n()) %>%
                    dplyr::mutate(n_cis_activated_samples_with_mut_across_cohort = sum(n_muts_of_mut_type > 0)) %>%
                    dplyr::ungroup()

                alpha_scale_name <-
                    generate_alpha_scale_label_for_cis_activated_samples_per_gene_bar_plot(matched_TAD, mutation_type)

                for (sort_by in c("occurrence", "mutations")) {
                    if (sort_by == "occurrence") {
                        temp_df3 <- temp_df2 %>%
                            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_across_cohort))
                    }
                    if (sort_by == "mutations") {
                        temp_df3 <- temp_df2 %>%
                            dplyr::mutate(gene_name = forcats::fct_reorder(gene_name, n_cis_activated_samples_with_mut_across_cohort))
                    }

                    for (filter_occurrence in c(1, 2, 3)) {
                        temp_df4 <- temp_df3 %>%
                            dplyr::filter(n_cis_activated_samples_across_cohort >= filter_occurrence)

                        for (filter_mutations in c(0, 1)) {
                            temp_df5 <- temp_df4 %>%
                                dplyr::filter(n_cis_activated_samples_with_mut_across_cohort >= filter_mutations)



                            for (pc_only in c(TRUE, FALSE)) {
                                final_df <- temp_df5 %>%
                                    dplyr::filter(gene_type == "protein_coding" | (!pc_only))

                                cis_activated_samples_per_gene_across_cohort_bar_plot_path <- file.path(HTML_report_figure_directory, paste0("cis_activated_samples", ifelse(copy_number_normalized, "__copy_number_normalized_", ""), "_per_gene_across_cohort_bar_plot___matched_TAD_", matched_TAD, "__mutation_type_", mutation_type, "__filter_occurrence_", filter_occurrence, "__filter_mutations_", filter_mutations, "__sort_", sort_by, ifelse(pc_only, "__pc_only", ""), ".svg"))

                                plot_cis_activated_samples_per_gene_bar_plot_integrated(final_df, color_palette, alpha_scale_name) %>%
                                    save_ggplot(cis_activated_samples_per_gene_across_cohort_bar_plot_path)

                                p()
                            }
                        }
                    }
                }
            })
        }
    })
}

create_HTML_for_cis_activated_samples_per_gene_bar_plots <- function() {
    # a.k. toggle plot
    toggle_plot_HTML <- ""
    for (mutation_type in c("all", "SV_only", "CNA_only", "somatic_SNV_only", "relevant_somatic_SNV_only")) {
        plot_card_title <- dplyr::case_when(
            (mutation_type == "all") ~ "Mutations highlighted",
            (mutation_type == "SV_only") ~ "SVs highlighted",
            (mutation_type == "CNA_only") ~ "CNAs highlighted",
            (mutation_type == "somatic_SNV_only") ~ "somatic SNVs highlighted",
            (mutation_type == "relevant_somatic_SNV_only") ~ "somatic SNVs with relevant TF binding site highlighted",
            TRUE ~ ""
        )

        toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
            <div class = "card my-4">
                        <h5 class="card-header">Detected Cis-Activated Samples Per Gene <br/>', plot_card_title, '
                            <span class="badge bg-danger float-end">Results</span>
                        </h5>
                        <div class="card-body p-3 border-bottom">
                            <p>Applied Outlier High Expression</p>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-copy_number_normalized_false"
                            data-hide-class="toggleplot_', mutation_type, '-copy_number_normalized_true"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_true"
                            id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_false"
                            disabled>Standard</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-copy_number_normalized_true"
                            data-hide-class="toggleplot_', mutation_type, '-copy_number_normalized_false"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_false"
                            id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_true"
                            >Copy Number Normalized</button>
                        </div>
                        <div class="card-body p-3 border-bottom">
                            <p>Mutations matched via</p>
                            <button
                                class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-toggleplot_', mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, "-matched_TAD_NOT_anywhere btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                                data-show-class="toggleplot_', mutation_type, '-matched_TAD_gene"
                                data-hide-class="toggleplot_', mutation_type, '-matched_TAD_NOT_gene"
                                data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_gene"
                                id="btn-show-toggleplot_', mutation_type, '-matched_TAD_gene"
                                >Gene TAD</button>
                            <button
                                class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-toggleplot_', mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_anywhere btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                                data-show-class="toggleplot_', mutation_type, '-matched_TAD_genehancer"
                                data-hide-class="toggleplot_', mutation_type, '-matched_TAD_NOT_genehancer"
                                data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_genehancer"
                                id="btn-show-toggleplot_', mutation_type, '-matched_TAD_genehancer"
                                >Genehancer (TAD)</button>
                            <button
                                class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-toggleplot_', mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                                data-show-class="toggleplot_', mutation_type, '-matched_TAD_anywhere"
                                data-hide-class="toggleplot_', mutation_type, '-matched_TAD_NOT_anywhere"
                                data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_anywhere"
                                id="btn-show-toggleplot_', mutation_type, '-matched_TAD_anywhere"
                                disabled>BOTH</button>
                            <button
                                class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-toggleplot_', mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, '-matched_TAD_NOT_anywhere"
                                data-show-class="toggleplot_', mutation_type, '-matched_TAD_chipseq"
                                data-hide-class="toggleplot_', mutation_type, '-matched_TAD_NOT_chipseq"
                                data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_chipseq"
                                id="btn-show-toggleplot_', mutation_type, '-matched_TAD_chipseq"
                                >Gene TAD + Mutation affects ChIP-Seq region</button>
                        </div>
                        <div class="card-body p-3 border-bottom">
                            <p>Sort by</p>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-sort_by_occurrence"
                            data-hide-class="toggleplot_', mutation_type, '-sort_by_mutations"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-sort_by_mutations"
                            id="btn-show-toggleplot_', mutation_type, '-sort_by_occurrence"
                            disabled>Cis-activated samples</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-sort_by_mutations"
                            data-hide-class="toggleplot_', mutation_type, '-sort_by_occurrence"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-sort_by_occurrence"
                            id="btn-show-toggleplot_', mutation_type, '-sort_by_mutations"
                            >Detected mutations</button>
                        </div>
                        <div class="card-body p-3 border-bottom">
                            <p>Filter by</p>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-show-toggleplot_', mutation_type, "-filter_occurrence_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_NOT_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_NOT_3"
                            data-show-class="toggleplot_', mutation_type, '-filter_occurrence_1"
                            data-hide-class="toggleplot_', mutation_type, '-filter_occurrence_NOT_1"
                            data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_1"
                            id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_1"
                            >&gt;= 1 cis-activated samples</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-show-toggleplot_', mutation_type, "-filter_occurrence_NOT_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_NOT_3"
                            data-show-class="toggleplot_', mutation_type, '-filter_occurrence_2"
                            data-hide-class="toggleplot_', mutation_type, '-filter_occurrence_NOT_2"
                            data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_2"
                            id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_2"
                            disabled>&gt;= 2 cis-activated samples</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable btn-show-toggleplot_', mutation_type, "-filter_occurrence__NOT_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_NOT_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_3"
                            data-show-class="toggleplot_', mutation_type, '-filter_occurrence_3"
                            data-hide-class="toggleplot_', mutation_type, '-filter_occurrence_NOT_3"
                            data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_3"
                            id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_3"
                            >&gt;= 3 cis-activated samples</button>
                        </div>
                        <div class="card-body p-3 border-bottom">
                            <p>Filter by</p>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-filter_mutations_0"
                            data-hide-class="toggleplot_', mutation_type, '-filter_mutations_1"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-filter_mutations_1"
                            id="btn-show-toggleplot_', mutation_type, '-filter_mutations_0"
                            disabled>&gt;= 0 detected mutations</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-filter_mutations_1"
                            data-hide-class="toggleplot_', mutation_type, '-filter_mutations_0"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-filter_mutations_0"
                            id="btn-show-toggleplot_', mutation_type, '-filter_mutations_1"
                            >&gt;= 1 detected mutations</button>
                        </div>
                        <div class="card-body p-3 border-bottom">
                            <p>Filter by</p>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-pc_only_false"
                            data-hide-class="toggleplot_', mutation_type, '-pc_only_true"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-pc_only_true"
                            id="btn-show-toggleplot_', mutation_type, '-pc_only_false"
                            disabled>All Gene Types</button>
                            <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="toggleplot_', mutation_type, '-pc_only_true"
                            data-hide-class="toggleplot_', mutation_type, '-pc_only_false"
                            data-enable-id="btn-show-toggleplot_', mutation_type, '-pc_only_false"
                            id="btn-show-toggleplot_', mutation_type, '-pc_only_true"
                            >Protein coding genes only</button>
                        </div>
                        <div class="p-3 text-center">
        ')
        for (copy_number_normalized in c(TRUE, FALSE)) {
            copy_number_normalized_string <- ifelse(copy_number_normalized, "true", "false")
            toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
            <div class="toggleplot_', mutation_type, "-copy_number_normalized_", copy_number_normalized_string, '" ', ifelse(!(copy_number_normalized_string == "false"), 'style="display: none;" ', ""), ">")

            for (matched_TAD in c("anywhere", "gene", "genehancer", "chipseq")) {
                toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "
                <div
                    ", ifelse(matched_TAD == "anywhere", paste0('class="toggleplot_', mutation_type, "-matched_TAD_anywhere toggleplot_", mutation_type, "-matched_TAD_NOT_gene toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq" '), ""), "
                    ", ifelse(matched_TAD == "gene", paste0('class="toggleplot_', mutation_type, "-matched_TAD_NOT_anywhere toggleplot_", mutation_type, "-matched_TAD_gene toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq" '), ""), "
                    ", ifelse(matched_TAD == "genehancer", paste0('class="toggleplot_', mutation_type, "-matched_TAD_NOT_anywhere toggleplot_", mutation_type, "-matched_TAD_NOT_gene toggleplot_", mutation_type, "-matched_TAD_genehancer toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq" '), ""), "
                    ", ifelse(matched_TAD == "chipseq", paste0('class="toggleplot_', mutation_type, "-matched_TAD_NOT_anywhere toggleplot_", mutation_type, "-matched_TAD_NOT_gene toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer toggleplot_", mutation_type, '-matched_TAD_chipseq" '), ""), "
                    ", ifelse(!(matched_TAD == "anywhere"), 'style="display: none;" ', ""), "
                    >")

                for (filter_occurrence in c(1, 2, 3)) {
                    toggle_plot_HTML <- paste0(
                        toggle_plot_HTML, "\n", '
                    <div class="',
                        ifelse(filter_occurrence == 1, paste0("toggleplot_", mutation_type, "-filter_occurrence_1"), paste0("toggleplot_", mutation_type, "-filter_occurrence_NOT_1")), " ",
                        ifelse(filter_occurrence == 2, paste0("toggleplot_", mutation_type, "-filter_occurrence_2"), paste0("toggleplot_", mutation_type, "-filter_occurrence_NOT_2")), " ",
                        ifelse(filter_occurrence == 3, paste0("toggleplot_", mutation_type, "-filter_occurrence_3"), paste0("toggleplot_", mutation_type, "-filter_occurrence_NOT_3")),
                        '" ',
                        ifelse(!(filter_occurrence == 2), 'style="display: none;" ', ""),
                        ">"
                    )

                    for (filter_mutations in c(0, 1)) {
                        toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
                        <div class="toggleplot_', mutation_type, "-filter_mutations_", filter_mutations, '" ', ifelse(!(filter_mutations == 0), 'style="display: none;" ', ""), ">")

                        for (sort_by in c("occurrence", "mutations")) {
                            toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
                            <div class="toggleplot_', mutation_type, "-sort_by_", sort_by, '" ', ifelse(!(sort_by == "occurrence"), 'style="display: none;" ', ""), ">")

                            for (pc_only in c(TRUE, FALSE)) {
                                pc_only_string <- ifelse(pc_only, "true", "false")
                                toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
                                <div class="toggleplot_', mutation_type, "-pc_only_", pc_only_string, '" ', ifelse(!(pc_only_string == "false"), 'style="display: none;" ', ""), ">")

                                # IMAGE BODY
                                toggleplot_this_img_file_name <- paste0("cis_activated_samples", ifelse(copy_number_normalized, "__copy_number_normalized_", ""), "_per_gene_across_cohort_bar_plot___matched_TAD_", matched_TAD, "__mutation_type_", mutation_type, "__filter_occurrence_", filter_occurrence, "__filter_mutations_", filter_mutations, "__sort_", sort_by, ifelse(pc_only, "__pc_only", ""), ".svg")
                                toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '<img loading="lazy" class ="w-100" height="672" width="672" src="figures/', toggleplot_this_img_file_name, '" />')

                                toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
                            }
                            toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
                        }
                        toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
                    }
                    toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
                }
                toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
            }
            toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div>")
        }
        toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", "</div></div>")
    }

    return(toggle_plot_HTML)
}

create_HTML_for_cis_activated_samples_per_gene_bar_plots_NEW <- function() {
    toggle_plot_HTML <- '
    <script>
        class ToggleImage{
            constructor(targetElement, mutationType){
                this.targetElement = targetElement
                this.mutation_type = mutationType
                this.copy_number_normalized = "false"
                this.matched_TAD = "anywhere"
                this.filter_occurrence = "2"
                this.filter_mutations = "0"
                this.sort_by = "occurrence"
                this.pc_only = "false"
                this.targetId = targetElement.id
                this.updateImgSource()
            }

            updateImgSource(){
                this.targetElement.src = this.getImgSource()
            }

            getImgSource(){
                const source = `figures/cis_activated_samples${this.copy_number_normalized =="true" ? "__copy_number_normalized_" : ""}_per_gene_across_cohort_bar_plot___matched_TAD_${this.matched_TAD}__mutation_type_${this.mutation_type}__filter_occurrence_${this.filter_occurrence}__filter_mutations_${this.filter_mutations}__sort_${this.sort_by}${this.pc_only =="true" ? "__pc_only" : ""}.svg`
                return source
            }

            setCopy_number_normalized(copy_number_normalized){
                this.copy_number_normalized = copy_number_normalized
                this.updateImgSource()
            }
            setMatched_TAD(matched_TAD){
                this.matched_TAD = matched_TAD
                this.updateImgSource()
            }
            setFilter_occurrence(filter_occurrence){
                this.filter_occurrence = filter_occurrence
                this.updateImgSource()
            }
            setFilter_mutations(filter_mutations){
                this.filter_mutations = filter_mutations
                this.updateImgSource()
            }
            setSort_by(sort_by){
                this.sort_by = sort_by
                this.updateImgSource()
            }
            setPc_only(pc_only){
                this.pc_only = pc_only
                this.updateImgSource()
            }


            setAttribute(attribute, value){
                if(attribute==="copy_number_normalized"){
                    this.setCopy_number_normalized(value);
                    return;
                }
                if(attribute==="matched_TAD"){
                    this.setMatched_TAD(value);
                    return;
                }
                if(attribute==="filter_occurrence"){
                    this.setFilter_occurrence(value);
                    return;
                }
                if(attribute==="filter_mutations"){
                    this.setFilter_mutations(value);
                    return;
                }
                if(attribute==="sort_by"){
                    this.setSort_by(value);
                    return;
                }
                if(attribute==="pc_only"){
                    this.setPc_only(value);
                    return;
                }
                console.error("ERROR WRONG ATTRIBUTE")
            }
        }
    </script>

    '
    for (mutation_type in c("all", "SV_only", "CNA_only", "somatic_SNV_only", "relevant_somatic_SNV_only")) {
        plot_card_title <- dplyr::case_when(
            (mutation_type == "all") ~ "Mutations highlighted",
            (mutation_type == "SV_only") ~ "SVs highlighted",
            (mutation_type == "CNA_only") ~ "CNAs highlighted",
            (mutation_type == "somatic_SNV_only") ~ "somatic SNVs highlighted",
            (mutation_type == "relevant_somatic_SNV_only") ~ "somatic SNVs with relevant TF binding site highlighted",
            TRUE ~ ""
        )

        toggle_plot_HTML <- paste0(toggle_plot_HTML, "\n", '
            <div class = "card my-4">
                <h5 class="card-header">Detected Cis-Activated Samples Per Gene <br/>', plot_card_title, '
                    <span class="badge bg-primary float-end">Input data</span>
                </h5>
                <div class="card-body p-3 border-bottom">
                    <p>Applied Outlier High Expression</p>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="-copy_number_normalized"
                    data-targetvalue="false"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_true"
                    id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_false"
                    disabled>Standard</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="copy_number_normalized"
                    data-targetvalue="true"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_false"
                    id="btn-show-toggleplot_', mutation_type, '-copy_number_normalized_true"
                    >Copy Number Normalized</button>
                </div>
                <div class="card-body p-3 border-bottom">
                    <p>Mutations matched via</p>
                    <button
                        class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, "-matched_TAD_NOT_anywhere btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                        data-targetattribute="matched_TAD"
                        data-targetvalue="gene"
                        data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_gene"
                        id="btn-show-toggleplot_', mutation_type, '-matched_TAD_gene"
                        >Gene TAD</button>
                    <button
                        class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-toggleplot_", mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_anywhere btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                        data-targetattribute="matched_TAD"
                        data-targetvalue="genehancer"
                        data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_genehancer"
                        id="btn-show-toggleplot_', mutation_type, '-matched_TAD_genehancer"
                        >Genehancer (TAD)</button>
                    <button
                        class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-toggleplot_", mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, '-matched_TAD_NOT_chipseq"
                        data-targetattribute="matched_TAD"
                        data-targetvalue="anywhere"
                        data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_anywhere"
                        id="btn-show-toggleplot_', mutation_type, '-matched_TAD_anywhere"
                        disabled>BOTH</button>
                    <button
                        class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-toggleplot_", mutation_type, "-matched_TAD_NOT_gene btn-toggleplot_", mutation_type, "-matched_TAD_NOT_genehancer btn-toggleplot_", mutation_type, '-matched_TAD_NOT_anywhere"
                        data-targetattribute="matched_TAD"
                        data-targetvalue="chipseq"
                        data-enable-class="btn-toggleplot_', mutation_type, '-matched_TAD_NOT_chipseq"
                        id="btn-show-toggleplot_', mutation_type, '-matched_TAD_chipseq"
                        >Gene TAD + Mutation affects ChIP-Seq region</button>
                </div>
                <div class="card-body p-3 border-bottom">
                    <p>Sort by</p>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="sort_by"
                    data-targetvalue="occurrence"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-sort_by_mutations"
                    id="btn-show-toggleplot_', mutation_type, '-sort_by_occurrence"
                    disabled>Cis-activated samples</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="sort_by"
                    data-targetvalue="mutations"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-sort_by_occurrence"
                    id="btn-show-toggleplot_', mutation_type, '-sort_by_mutations"
                    >Detected mutations</button>
                </div>
                <div class="card-body p-3 border-bottom">
                    <p>Filter by</p>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-show-toggleplot_", mutation_type, "-filter_occurrence_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_NOT_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_NOT_3"
                    data-targetattribute="filter_occurrence"
                    data-targetvalue="1"
                    data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_1"
                    id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_1"
                    >&gt;= 1 cis-activated samples</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-show-toggleplot_", mutation_type, "-filter_occurrence_NOT_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_NOT_3"
                    data-targetattribute="filter_occurrence"
                    data-targetvalue="2"
                    data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_2"
                    id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_2"
                    disabled>&gt;= 2 cis-activated samples</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, " btn-js-disable-self btn-js-enable btn-show-toggleplot_", mutation_type, "-filter_occurrence__NOT_1 btn-show-toggleplot_", mutation_type, "-filter_occurrence_NOT_2 btn-show-toggleplot_", mutation_type, '-filter_occurrence_3"
                    data-targetattribute="filter_occurrence"
                    data-targetvalue="3"
                    data-enable-class="btn-show-toggleplot_', mutation_type, '-filter_occurrence_NOT_3"
                    id="btn-show-toggleplot_', mutation_type, '-filter_occurrence_3"
                    >&gt;= 3 cis-activated samples</button>
                </div>
                <div class="card-body p-3 border-bottom">
                    <p>Filter by</p>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="filter_mutations"
                    data-targetvalue="0"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-filter_mutations_1"
                    id="btn-show-toggleplot_', mutation_type, '-filter_mutations_0"
                    disabled>&gt;= 0 detected mutations</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="filter_mutations"
                    data-targetvalue="1"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-filter_mutations_0"
                    id="btn-show-toggleplot_', mutation_type, '-filter_mutations_1"
                    >&gt;= 1 detected mutations</button>
                </div>
                <div class="card-body p-3 border-bottom">
                    <p>Filter by</p>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="pc_only"
                    data-targetvalue="false"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-pc_only_true"
                    id="btn-show-toggleplot_', mutation_type, '-pc_only_false"
                    disabled>All Gene Types</button>
                    <button
                    class="btn btn-secondary btn-js-setsource-', mutation_type, ' btn-js-disable-self btn-js-enable"
                    data-targetattribute="pc_only"
                    data-targetvalue="true"
                    data-enable-id="btn-show-toggleplot_', mutation_type, '-pc_only_false"
                    id="btn-show-toggleplot_', mutation_type, '-pc_only_true"
                    >Protein coding genes only</button>
                </div>
                <div class="p-3 text-center">
                    <img id="toggle-image-', mutation_type, '" loading="lazy" class ="w-100" src="about:blank" />
                </div>
                <script>
                const imgToggler', mutation_type, ' = new ToggleImage(
                    targetElement = document.getElementById("toggle-image-', mutation_type, '"),
                    mutationType = "', mutation_type, '"
                )
                for (let btnElem of document.getElementsByClassName("btn-js-setsource-', mutation_type, '")){
                    const targetattribute = btnElem.dataset.targetattribute
                    const targetvalue = btnElem.dataset.targetvalue
                    const clickHandler = () => imgToggler', mutation_type, '.setAttribute(targetattribute, targetvalue)
                    btnElem.addEventListener("click",clickHandler)
                    console.log("LISTENER ADDED")
                }
                </script>
            </div>
        ')
    }

    return(toggle_plot_HTML)
}