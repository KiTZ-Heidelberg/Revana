tf_binding_site_step4_workflow <- function(sample_dir_path,
                                           sample_id,
                                           somatic_SNV_tf_binding_data_file_path,
                                           chipseq,
                                           cis_activation_summary_TADs,
                                           genehancer_with_cis_activation_summary_TADs,
                                           store_somatic_SNV_tf_binding_data_genes_all = FALSE,
                                           store_somatic_SNV_tf_binding_data_genehancer_all = FALSE) {
    # LINK SNVs TF BINDING DATA TO GENES AND GENEHANCERS ######################
    cat("LINK SNVs TF BINDING DATA TO GENES AND GENEHANCERS\n")

    # import required data -----------------------------------
    cat("reading and processing chunks...\n")
    somatic_SNV_tf_binding_data <- readRDS(somatic_SNV_tf_binding_data_file_path)

    # process data 1 - CHIPSEQ -------------------------------
    data.table::setDT(somatic_SNV_tf_binding_data)
    add_chipseq_info_to_SNV_tf_binding_data(somatic_SNV_tf_binding_data, chipseq)


    # chunk data ---------------------------------------------
    SNV_chunk_size <- 500000
    n_somatic_SNV_tf_binding_data_rows <- nrow(somatic_SNV_tf_binding_data)
    # somatic_SNV_tf_binding_data_chunked_list <- split(somatic_SNV_tf_binding_data, ((seq(nrow(somatic_SNV_tf_binding_data))-1) %/% SNV_chunk_size)+1)
    max_chunk_index <- ((n_somatic_SNV_tf_binding_data_rows - 1) %/% SNV_chunk_size) + 1

    # initialize lists for summaries -------------------------
    somatic_SNV_tf_binding_data_gene_summary_max_score_diff_chunk_list <- vector(mode = "list", length = max_chunk_index)
    somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_chunk_list <- vector(mode = "list", length = max_chunk_index)
    somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_chunk_list <- vector(mode = "list", length = max_chunk_index)
    somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs_chunk_list <- vector(mode = "list", length = max_chunk_index)
    somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs_chunk_list <- vector(mode = "list", length = max_chunk_index)

    # define output paths ------------------------------------
    somatic_SNV_tf_binding_data_genes_all_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.all.txt"))
    somatic_SNV_tf_binding_data_genes_cis_activated_only_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.txt"))
    somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.txt"))
    somatic_SNV_tf_binding_data_genehancer_all_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.all.txt"))
    somatic_SNV_tf_binding_data_genehancer_cis_activated_only_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.txt"))
    somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only_txt <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.txt"))

    # process data 2 - chunkwise workflow ----------------------
    for (chunk_i in seq_len(max_chunk_index)) {
        cat(paste0(chunk_i, " / ", max_chunk_index, "\n"))

        # run garbage control
        gc()


        if (chunk_i < max_chunk_index) {
            rows_of_chunk <- seq(from = ((chunk_i - 1) * SNV_chunk_size + 1), to = (chunk_i * SNV_chunk_size), by = 1)
        } else {
            # only include existing rows for last chunk
            rows_of_chunk <- seq(from = ((chunk_i - 1) * SNV_chunk_size + 1), to = n_somatic_SNV_tf_binding_data_rows, by = 1)
        }

        # link SNV TF data to genes ------------------------------------

        # via gene TAD
        somatic_SNV_tf_binding_data_gene_combinations <- link_SNV_tf_binding_data_to_genes(somatic_SNV_tf_binding_data[rows_of_chunk], cis_activation_summary_TADs)

        # via genehancer
        somatic_SNV_tf_binding_data_genehancer_combinations <- link_SNV_tf_binding_data_to_genehancer(somatic_SNV_tf_binding_data[rows_of_chunk], genehancer_with_cis_activation_summary_TADs)

        # summarize ----------------------------------------------------
        # keep only entries with max score diff => CHANGE DIMENSIONS: SNV * TF * gene -> gene
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff_chunk_list[[chunk_i]] <- somatic_SNV_tf_binding_data_gene_combinations %>%
            dplyr::group_by(gene_name) %>%
            dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup()

        # keep only entries with max score diff => CHANGE DIMENSIONS: SNV * TF * gene * genehancer -> gene
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_chunk_list[[chunk_i]] <- somatic_SNV_tf_binding_data_genehancer_combinations %>%
            dplyr::group_by(connected_gene) %>%
            dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(gene_name = connected_gene)

        # keep only entries with max score diff => CHANGE DIMENSIONS: SNV * TF * gene * chipseq -> gene
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_chunk_list[[chunk_i]] <- somatic_SNV_tf_binding_data_gene_combinations %>%
            dplyr::filter(overlaps_with_chipseq == TRUE) %>%
            dplyr::group_by(gene_name) %>%
            dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup()

        # count n relevant somatic SNVs
        somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs_chunk_list[[chunk_i]] <- somatic_SNV_tf_binding_data_gene_combinations %>%
            dplyr::group_by(gene_name) %>%
            dplyr::summarize(
                n_relevant_somatic_SNVs_via_gene_TAD = sum(relevant_tf_binding_site, na.rm = TRUE),
                n_relevant_somatic_SNVs_via_chipseq = sum((relevant_tf_binding_site & overlaps_with_chipseq), na.rm = TRUE)
            )

        somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs_chunk_list[[chunk_i]] <- somatic_SNV_tf_binding_data_genehancer_combinations %>%
            dplyr::rename(gene_name = connected_gene) %>%
            dplyr::group_by(gene_name) %>%
            dplyr::summarize(n_relevant_somatic_SNVs_via_genehancer = sum(relevant_tf_binding_site, na.rm = TRUE))

        # store results - chunkwise ------------------------------------

        # check for existing files to prevent appending to old data
        if (chunk_i == 1) {
            if (store_somatic_SNV_tf_binding_data_genes_all & file.exists(somatic_SNV_tf_binding_data_genes_all_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genes.all.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genes.all.txt - File already exists. Appending will lead to incorrect data")
            }

            if (file.exists(somatic_SNV_tf_binding_data_genes_cis_activated_only_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genes.cis_activated_only.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genes.cis_activated_only.txt - File already exists. Appending will lead to incorrect data")
            }

            if (file.exists(somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.txt - File already exists. Appending will lead to incorrect data")
            }

            if (store_somatic_SNV_tf_binding_data_genehancer_all & file.exists(somatic_SNV_tf_binding_data_genehancer_all_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genehancer.all.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genehancer.all.txt - File already exists. Appending will lead to incorrect data")
            }

            if (file.exists(somatic_SNV_tf_binding_data_genehancer_cis_activated_only_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genehancer.cis_activated_only.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.txt - File already exists. Appending will lead to incorrect data")
            }

            if (file.exists(somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only_txt)) {
                cat("ERROR: .somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.txt - File already exists. Appending will lead to incorrect data\n")
                stop(".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.txt - File already exists. Appending will lead to incorrect data")
            }
        }

        # save somatic_SNV_tf_binding_data_genes_all
        if (store_somatic_SNV_tf_binding_data_genes_all) {
            saveRDS(
                somatic_SNV_tf_binding_data_gene_combinations,
                file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.all.__chunk_", chunk_i, ".Rds"))
            )

            data.table::fwrite(
                somatic_SNV_tf_binding_data_gene_combinations,
                file = somatic_SNV_tf_binding_data_genes_all_txt,
                sep = "\t",
                append = TRUE
            )
        }

        # save somatic_SNV_tf_binding_data_genes_cis_activated_only
        saveRDS(
            somatic_SNV_tf_binding_data_gene_combinations[cis_activated_gene == TRUE],
            file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.__chunk_", chunk_i, ".Rds"))
        )
        data.table::fwrite(
            somatic_SNV_tf_binding_data_gene_combinations[cis_activated_gene == TRUE],
            file = somatic_SNV_tf_binding_data_genes_cis_activated_only_txt,
            sep = "\t",
            append = TRUE
        )

        # save somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only
        saveRDS(
            somatic_SNV_tf_binding_data_gene_combinations[(cis_activated_gene == TRUE) & (relevant_tf_binding_site == TRUE)],
            file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genes.cis_activated_only.relevant_only.__chunk_", chunk_i, ".Rds"))
        )
        data.table::fwrite(
            somatic_SNV_tf_binding_data_gene_combinations[(cis_activated_gene == TRUE) & (relevant_tf_binding_site == TRUE)],
            file = somatic_SNV_tf_binding_data_genes_cis_activated_only_relevant_only_txt,
            sep = "\t",
            append = TRUE
        )
        rm(somatic_SNV_tf_binding_data_gene_combinations)
        gc() # free memory


        # save somatic_SNV_tf_binding_data_genehancer_all
        if (store_somatic_SNV_tf_binding_data_genehancer_all) {
            saveRDS(
                somatic_SNV_tf_binding_data_genehancer_combinations,
                file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.all.__chunk_", chunk_i, ".Rds"))
            )

            data.table::fwrite(
                somatic_SNV_tf_binding_data_genehancer_combinations,
                file = somatic_SNV_tf_binding_data_genehancer_all_txt,
                sep = "\t",
                append = TRUE
            )
        }

        # save somatic_SNV_tf_binding_data_genehancer_cis_activated_only
        saveRDS(
            somatic_SNV_tf_binding_data_genehancer_combinations[cis_activated_gene == TRUE],
            file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.__chunk_", chunk_i, ".Rds"))
        )
        data.table::fwrite(
            somatic_SNV_tf_binding_data_genehancer_combinations[cis_activated_gene == TRUE],
            file = somatic_SNV_tf_binding_data_genehancer_cis_activated_only_txt,
            sep = "\t",
            append = TRUE
        )
        # save somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only
        saveRDS(
            somatic_SNV_tf_binding_data_genehancer_combinations[(cis_activated_gene == TRUE) & (relevant_tf_binding_site == TRUE)],
            file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.genehancer.cis_activated_only.relevant_only.__chunk_", chunk_i, ".Rds"))
        )
        data.table::fwrite(
            somatic_SNV_tf_binding_data_genehancer_combinations[(cis_activated_gene == TRUE) & (relevant_tf_binding_site == TRUE)],
            file = somatic_SNV_tf_binding_data_genehancer_cis_activated_only_relevant_only_txt,
            sep = "\t",
            append = TRUE
        )
        rm(somatic_SNV_tf_binding_data_genehancer_combinations)
        gc() # free memory
    }

    # summarize chunks ----------------------------
    somatic_SNV_tf_binding_data_gene_summary_max_score_diff <- data.table::rbindlist(somatic_SNV_tf_binding_data_gene_summary_max_score_diff_chunk_list) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff <- data.table::rbindlist(somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_chunk_list) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff <- data.table::rbindlist(somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_chunk_list) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::slice_max(order_by = score_diff, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    # summarize n of relevant somatic SNVs
    somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs <- data.table::rbindlist(somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs_chunk_list) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarize(
            n_relevant_somatic_SNVs_via_gene_TAD = sum(n_relevant_somatic_SNVs_via_gene_TAD, na.rm = TRUE),
            n_relevant_somatic_SNVs_via_chipseq = sum(n_relevant_somatic_SNVs_via_chipseq, na.rm = TRUE),
        )

    somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs <- data.table::rbindlist(somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs_chunk_list) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarize(n_relevant_somatic_SNVs_via_genehancer = sum(n_relevant_somatic_SNVs_via_genehancer, na.rm = TRUE))

    cis_activation_summary_with_n_relevant_somatic_SNVs <- cis_activation_summary_TADs %>%
        dplyr::left_join(somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs, by = "gene_name") %>%
        dplyr::left_join(somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs, by = "gene_name") %>%
        # EDIT 3 Mai 14:21 UNTESTED
        tidyr::replace_na(list(n_relevant_somatic_SNVs_via_gene_TAD = 0, n_relevant_somatic_SNVs_via_chipseq = 0, n_relevant_somatic_SNVs_via_genehancer = 0))

    # write summarized chunks ----------------------------
    #
    saveRDS(
        cis_activation_summary_with_n_relevant_somatic_SNVs,
        file = file.path(sample_dir_path, paste0(sample_id, ".cis_activation_summary_with_n_relevant_somatic_SNVs.Rds"))
    )
    readr::write_tsv(
        cis_activation_summary_with_n_relevant_somatic_SNVs,
        file = file.path(sample_dir_path, paste0(sample_id, ".cis_activation_summary_with_n_relevant_somatic_SNVs.txt"))
    )

    #
    saveRDS(
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_gene_summary_max_score_diff.Rds"))
    )

    readr::write_tsv(
        somatic_SNV_tf_binding_data_gene_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_gene_summary_max_score_diff.txt")),
    )

    #
    saveRDS(
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff.Rds"))
    )

    readr::write_tsv(
        somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff.txt")),
    )

    #
    saveRDS(
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff.Rds"))
    )

    readr::write_tsv(
        somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff.txt")),
    )

    rm(somatic_SNV_tf_binding_data, n_somatic_SNV_tf_binding_data_rows, max_chunk_index, cis_activation_summary_with_n_relevant_somatic_SNVs, somatic_SNV_tf_binding_data_gene_summary_max_score_diff, somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff, somatic_SNV_tf_binding_data_gene_summary_max_score_diff_chunk_list, somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_chunk_list, somatic_SNV_tf_binding_data_gene_summary_n_relevant_somatic_SNVs_chunk_list, somatic_SNV_tf_binding_data_genehancer_summary_n_relevant_somatic_SNVs_chunk_list)
    gc()
}