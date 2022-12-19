create_plots_part_1_and_2 <- function(processed_data,
                                      HTML_report_figure_directory,
                                      skip_creating_toggle_images = FALSE,
                                      color_palette,
                                      use_parallelization = FALSE) {

    
    # Input data -----------------------------------
    cat("Creating input data plots...\n")
    progressr::with_progress({
        p <- progressr::progressor(steps = 8)

        # subgroups bar chart
        plot_subgroups_bar_chart(processed_data, color_palette) %>%
            save_ggplot(file.path(HTML_report_figure_directory, "subgroups_bar_chart.png"))

        p()


        # markers bar chart
        markers_bar_chart_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("marker_bar_chart_bygroup_", subgroup, ".png"))
                plot_markers_bar_chart_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # somatic SNV bar chart
        somatic_SNV_bar_chart_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("somatic_SNV_bar_chart_bygroup_", subgroup, ".png"))
                plot_somatic_SNV_bar_chart_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # SV bar chart
        SV_bar_chart_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("SV_bar_chart_bygroup_", subgroup, ".png"))
                plot_SV_bar_chart_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # CNA bar chart
        CNA_bar_chart_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("CNA_bar_chart_bygroup_", subgroup, ".png"))
                plot_CNA_bar_chart_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # expression box plot
        expression_box_plot_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("expression_box_plot_bygroup_", subgroup, ".png"))
                plot_expression_box_plot_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # PCA expression plot
        expression_PCA_plot_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("expression_PCA_plot_bygroup_", subgroup, ".png"))
                plot_expression_PCA_plot_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

        # cis-activated genes box plot list
        cis_activated_genes_bar_chart_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("cis_activated_genes_bar_chart_bygroup_", subgroup, ".png"))
                plot_cis_activated_genes_bar_chart_by_subgroup(processed_data, subgroup, color_palette) %>%
                    save_ggplot(path)
                return(path)
            })

        p()

    })


    # Results overview -----------------------------------
    
    cat("Creating results overview plots...\n")
    progressr::with_progress({
        p <- progressr::progressor(steps = 5)

        # ASE detection - upset plot - by group
        ase_detection_upset_plot_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("ase_detection_upset_plot_bygroup_", subgroup, ".svg"))
                plot <- plot_ase_detection_upset_plot_by_subgroup(processed_data, subgroup)

                # cat("Saving UpSetR image\n")
                svglite::svglite(path)
                print(plot)
                dev.off()
                return(path)
            })

        p()

        # ASE detection - upset plot - whole cohort
        ase_detection_upset_plot_whole_cohort_path <- file.path(HTML_report_figure_directory, "ase_detection_upset_plot_whole_cohort.svg")
        ase_detection_upset_plot_whole_cohort_plot <- plot_ase_detection_upset_plot_whole_cohort(processed_data)

        # cat("Saving UpSetR image\n")
        svglite::svglite(ase_detection_upset_plot_whole_cohort_path)
        print(ase_detection_upset_plot_whole_cohort_plot)
        dev.off()

        p()

        # applied filters - upset plot - by group
        applied_filters_upset_plot_bygroup_path <- purrr::set_names(processed_data$subgroups) %>%
            purrr::map(function(subgroup) {
                path <- file.path(HTML_report_figure_directory, paste0("applied_filters_upset_plot_bygroup_", subgroup, ".svg"))
                plot <- plot_applied_filters_upset_plot_by_subgroup(processed_data, subgroup)

                # cat("Saving UpSetR image\n")
                svglite::svglite(path)
                print(plot)
                dev.off()
                return(path)
            })

        p()

        # applied filters - upset plot - whole cohort
        applied_filters_upset_plot_whole_cohort_path <- file.path(HTML_report_figure_directory, "applied_filters_upset_plot_whole_cohort.svg")
        applied_filters_upset_plot_whole_cohort_plot <- plot_applied_filters_upset_plot_whole_cohort(processed_data)

        # cat("Saving UpSetR image\n")
        svglite::svglite(applied_filters_upset_plot_whole_cohort_path)
        print(applied_filters_upset_plot_whole_cohort_plot)
        dev.off()

        p()

        # OHE and ASE q-values
        plot_OHE_ASE_q_value_dot_plot(processed_data, color_palette) %>%
            save_ggplot(file.path(HTML_report_figure_directory, "OHE_ASE_q_value_dot_plot.png"))

        p()


    })

    if (!skip_creating_toggle_images) {
        # detected cis-activated samples per gene across cohort
        time_before <- Sys.time()

        if(use_parallelization == TRUE){
            create_all_cis_activated_samples_per_gene_bar_plots_PARALLELIZATION(processed_data, color_palette, HTML_report_figure_directory)
        }else{
            create_all_cis_activated_samples_per_gene_bar_plots(processed_data, color_palette, HTML_report_figure_directory)
        }
        
        # time_after <- Sys.time()
        # print(time_after- time_before)
    }
}

create_plots_and_HTML_part_4 <- function(processed_data,
                                         color_palette,
                                         HTML_report_figure_directory,
                                         HTML_report_subdocument_directory,
                                         max_n_genes_to_analyze_further,
                                         additional_genes_to_analyze_further,
                                         add_custom_HTML_to_gene_subdocument = NULL,
                                         add_custom_HTML_to_gene_subdocument_sample_data = NULL,
                                         genome_name = "hg19",
                                         use_parallelization = FALSE) {
    genes_to_analyze_further_ordered <- processed_data$cis_activated_genes_by_gene %>%
        dplyr::filter(
            ((n_cis_activated_samples >= 3) & (gene_type == "protein_coding")) |
                ((n_CA_plus_SV + n_CA_plus_CNA) > 0)
        ) %>%
        # dplyr::filter(gene_type == 'protein_coding') %>%
        dplyr::arrange(dplyr::desc(n_cis_activated_samples)) %>%
        dplyr::pull(gene_name)

    if (max_n_genes_to_analyze_further > length(genes_to_analyze_further_ordered)) {
        n_genes_to_analyze_further <- length(genes_to_analyze_further_ordered)
    } else {
        n_genes_to_analyze_further <- max_n_genes_to_analyze_further
    }

    genes_to_analyze_further <- c(additional_genes_to_analyze_further, genes_to_analyze_further_ordered[seq_len(n_genes_to_analyze_further)])

    # results by gene LOOP -------------------------------

    n_genes_to_actually_analyze_further <- length(genes_to_analyze_further)
    
    if(use_parallelization){
        cat("Creating by-gene analysis... (on ")
        cat(future::nbrOfWorkers())
        cat(" workers)\n")
        future::plan(future::multisession)

        progressr::with_progress({
            p <- progressr::progressor(steps = n_genes_to_actually_analyze_further)

            part4_html_by_gene <- genes_to_analyze_further %>% furrr::future_map(function(GENE_NAME) {
                html_this_gene <- create_report_gene_subdocument(
                    processed_data,
                    GENE_NAME,
                    color_palette,
                    HTML_report_subdocument_directory,
                    add_custom_HTML_to_gene_subdocument,
                    add_custom_HTML_to_gene_subdocument_sample_data,
                    genome_name
                )
                p()
                return(html_this_gene)
            },
            # adresses: <RngFutureWarning: UNRELIABLE VALUE: Future ('<none>') unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore".>
            .options = furrr::furrr_options(seed = TRUE))
        })
    }else{
        cat("Creating by-gene analysis ...\n")

        progressr::with_progress({
            p <- progressr::progressor(steps = n_genes_to_actually_analyze_further)

            part4_html_by_gene <- genes_to_analyze_further %>% purrr::map(function(GENE_NAME) {
                html_this_gene <- create_report_gene_subdocument(
                    processed_data,
                    GENE_NAME,
                    color_palette,
                    HTML_report_subdocument_directory,
                    add_custom_HTML_to_gene_subdocument,
                    add_custom_HTML_to_gene_subdocument_sample_data,
                    genome_name
                )
                p()
                return(html_this_gene)
            })
        })
    }
    

    # write gene HTML to subdocument files -----------------------------------
    part4_html_by_gene %>% purrr::walk(function(GENE_HTML_DATA) {
        subdocument_path <- GENE_HTML_DATA$subdocument_path %>% create_file_if_missing()
        write_HTML_to_file(
            html_content = GENE_HTML_DATA$subdocument_HTML,
            write_path = subdocument_path,
            print_before = paste0("    ",GENE_HTML_DATA$gene_name, " - Writing HTML subdocument to file..."),
            print_after = "DONE \n"
        )
    })

    # wrap gene overview HTML -----------------------------------
    analysis_by_gene_HTML <- part4_html_by_gene %>%
        purrr::map(~ .x[["overview_HTML"]]) %>%
        unlist() %>%
        paste0(collapse = "\n") %>%
        {
            paste0('<div class="list-group">', ., "</div>")
        }

    return(analysis_by_gene_HTML)
}

#' Create Interactive HTML report for Revana results
#'
#' @param HTML_report_output_dir_path path to the directory where the HTML report and its figures should be stored
#' @param output_paths_file_path path to the results path file - created by revana::run() . If several subgroups should be summarized the respective results paths files need to be merged in advance
#' @param use_cache debugging option to use cache for loaded data
#' @param create_cache debugging option to create cache for loaded data, if not already in place
#' @param max_n_genes_to_analyze_further max number of genes that are analyzed in 'by gene' analysis
#' @param additional_genes_to_analyze_further additional genes that are to be included in 'by gene' analysis
#' @param skip_creating_toggle_images should the cis-activated samples per gene - alias toggle plots - be created? defaults to TRUE
#' @param has_run_tf_binding_site_analysis Was TF binding site analysis conducted?
#' @param add_custom_HTML_to_gene_subdocument plugin function for custom HTML in by gene subdocument ()
#' @param add_custom_HTML_to_gene_subdocument_sample_data plugin function for custom HTML in by gene subdocument - sample data section
#' @param genome_name name of the reference genome - tested with "hg19" and "hg38"
#' @param use_parallelization should parallelization be used to speed up report generation (TRUE / FALSE)
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#'
#'
#' @export

create_results_HTML_report <- function(HTML_report_output_dir_path, output_paths_file_path, use_cache = F, create_cache = F, max_n_genes_to_analyze_further = 500, additional_genes_to_analyze_further = character(0), skip_creating_toggle_images = F, has_run_tf_binding_site_analysis = F, add_custom_HTML_to_gene_subdocument = NULL, add_custom_HTML_to_gene_subdocument_sample_data = NULL, genome_name = "hg19", use_parallelization = FALSE, verbose = FALSE) {
    # INITIALIZATION #############################################

    cat("Generating HTML report...\n")

    # check paths and create files -------------------------------
    create_dir_if_missing(HTML_report_output_dir_path)
    HTML_report_figure_directory <- create_dir_if_missing(file.path(HTML_report_output_dir_path, "figures"))
    HTML_report_subdocument_directory <- create_dir_if_missing(file.path(HTML_report_output_dir_path, "genes"))
    HTML_report_output_path <- create_file_if_missing(file.path(HTML_report_output_dir_path, "report.html"))


    # DATA IMPORT ###################################################

    # use cache if available --------------------------------------
    cache_available <- FALSE
    cache_file_path <- file.path(HTML_report_output_dir_path, "cached_data.Rds")
    if (use_cache) {
        if (file.exists(cache_file_path)) {
            cache_available <- TRUE

            cat("using cached data from ")
            cat(cache_file_path)
            cat("\n")

            data <- readRDS(cache_file_path)
        }
    }

    # import data from results if cache not available --------------------------------------
    if (!cache_available) {
        data <- load_data_for_HTML_report_from_result_files(output_paths_file_path, has_run_tf_binding_site_analysis)

        if (create_cache) {
            # save cache
            saveRDS(data, file = cache_file_path)
        }
    }

    # DATA PROCESSING ###########################################
    processed_data <- process_data_for_HTML_report(data, has_run_tf_binding_site_analysis)


    # CREATE PLOTS AND HTML ####################################
    color_palette <- create_color_palette(processed_data$subgroups)

    # part 1/2
    create_plots_part_1_and_2(processed_data, HTML_report_figure_directory, skip_creating_toggle_images, color_palette, use_parallelization = use_parallelization)
    toggle_plot_HTML <- create_HTML_for_cis_activated_samples_per_gene_bar_plots_NEW()
    SV_recurrent_TAD_combinations_HTML <- create_HTML_for_SV_recurrent_TAD_combinations(processed_data)

    # part 3
    HTML_mut_analysis <- create_mut_analysis_HTML(processed_data, has_run_tf_binding_site_analysis, HTML_report_figure_directory)

    # part 4
    analysis_by_gene_HTML <- create_plots_and_HTML_part_4(
        processed_data,
        color_palette,
        HTML_report_figure_directory,
        HTML_report_subdocument_directory,
        max_n_genes_to_analyze_further,
        additional_genes_to_analyze_further,
        add_custom_HTML_to_gene_subdocument,
        add_custom_HTML_to_gene_subdocument_sample_data,
        genome_name,
        use_parallelization = use_parallelization
    )

    # MAIN DOCUMENT HTML ####################################
    main_document_HTML <- create_report_HTML(
        processed_data,
        toggle_plot_HTML,
        SV_recurrent_TAD_combinations_HTML,
        HTML_mut_analysis,
        analysis_by_gene_HTML
    )

    #  write HTML to file -------------------------------------
    write_HTML_to_file(main_document_HTML, HTML_report_output_path)
}