process_data_for_HTML_report <- function(data, has_run_tf_binding_site_analysis) {
    cat("Processing data...\n")

    # util functions -----------------------------
    binarify <- function(x) {
        return(
            ifelse(
                x & (!is.na(x)),
                1,
                0
            )
        )
    }



    # subgroups -----------------------------------
    subgroups_table <- data$result_paths %>%
        dplyr::group_by(subgroup) %>%
        dplyr::summarize(n = dplyr::n())

    subgroups <- subgroups_table$subgroup
    samples_by_subgroup <- purrr::set_names(subgroups) %>%
        purrr::map(~ (
            data$result_paths %>%
                dplyr::filter(subgroup == .x) %>%
                dplyr::pull(sample_id)
        ))


    # add groups to cis-activation-summary -----------------------------------
    cis_activation_summary_table_with_groups <- data$cis_activation_summary_table %>%
        dplyr::left_join(
            (data$result_paths %>%
                dplyr::select(sample_ID = sample_id, subgroup) %>%
                dplyr::distinct(sample_ID, .keep_all = TRUE)
            ),
            by = c("sample_ID")
        )

    # add gene_name col to gene genehancer combinations -----------------------------------
    data$SV_genehancer_all_table[, gene_name := connected_gene]
    data$CNA_genehancer_all_table[, gene_name := connected_gene]
    data$somatic_SNV_genehancer_all_table[, gene_name := connected_gene]
    if (has_run_tf_binding_site_analysis) {
        data$somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table[, gene_name := connected_gene]
    }

    # set DT for imported data -----------------------------
    data.table::setDT(data$SV_genehancer_all_table)
    data.table::setDT(data$SV_genehancer_all_table)
    data.table::setDT(data$SV_genehancer_all_table)


    # mutation counts by gene ------------------------------------------------------
    # merge tables
    SV_all_merged_table <- data.table::rbindlist(
        list(gene_TAD = data$SV_genes_all_table, genehancer_TAD = data$SV_genehancer_all_table),
        idcol = "source_table",
        fill = T
    ) %>%
        dplyr::select(sample_ID, gene_name, cis_activated_gene, sv_break1_chrom, sv_break1_pos, sv_break2_chrom, sv_break2_pos, sv_type, eventInversion, TAD_sv_break1_chrom, TAD_sv_break2_chrom, TAD_sv_break1_start, TAD_sv_break2_start, TAD_sv_break1_end, TAD_sv_break2_end, TAD_combination) %>%
        dplyr::distinct()

    if(nrow(SV_all_merged_table) > 0){
        GenomeInfoDb::seqlevelsStyle(SV_all_merged_table$sv_break1_chrom) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(SV_all_merged_table$sv_break2_chrom) <- "UCSC"
    }
    SV_all_merged_table <- SV_all_merged_table %>% dplyr::distinct()

    CNA_all_merged_table <- data.table::rbindlist(
        list(gene_TAD = data$CNA_genes_all_table, genehancer_TAD = data$CNA_genehancer_all_table),
        idcol = "source_table",
        fill = T
    ) %>%
        dplyr::select(sample_ID, gene_name, cis_activated_gene, cna_chrom, cna_start, cna_end, copy_number, CNA_type, log2) %>%
        dplyr::distinct()

    if(nrow(CNA_all_merged_table) > 0){
        GenomeInfoDb::seqlevelsStyle(CNA_all_merged_table$cna_chrom) <- "UCSC"
    }
    CNA_all_merged_table <- CNA_all_merged_table %>% dplyr::distinct()

    somatic_SNV_all_merged_table <- data.table::rbindlist(
        list(gene_TAD = data$somatic_SNV_genes_all_table, genehancer_TAD = data$somatic_SNV_genehancer_all_table),
        idcol = "source_table",
        fill = T
    ) %>%
        dplyr::select(sample_ID, gene_name, cis_activated_gene, snv_chrom, snv_pos, ref, alt) %>%
        dplyr::distinct()

    if(nrow(somatic_SNV_all_merged_table) > 0){
        GenomeInfoDb::seqlevelsStyle(somatic_SNV_all_merged_table$snv_chrom) <- "UCSC"
    }

    somatic_SNV_all_merged_table <- somatic_SNV_all_merged_table %>% dplyr::distinct()

    if (has_run_tf_binding_site_analysis) {
        somatic_SNV_tf_binding_data_cis_activated_only_merged_table <- data.table::rbindlist(list(gene_TAD = data$somatic_SNV_tf_binding_data_genes_cis_activated_only_table, genehancer_TAD = data$somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table), idcol = "source_table", fill = T)
    }

    # from gene TAD
    genes_n_SVs_from_gene_TAD_all_table <- data$SV_genes_all_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_SVs_from_gene_TAD = dplyr::n())

    genes_n_CNAs_from_gene_TAD_all_table <- data$CNA_genes_all_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_CNAs_from_gene_TAD = dplyr::n())

    genes_n_somatic_SNVs_from_gene_TAD_all_table <- data$somatic_SNV_genes_all_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_somatic_SNVs_from_gene_TAD = dplyr::n())

    if (has_run_tf_binding_site_analysis) {
        genes_n_relevant_somatic_SNVs_from_gene_TAD_cis_activated_only_table <- data$somatic_SNV_tf_binding_data_genes_cis_activated_only_table %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
            dplyr::summarize(n_relevant_tf_binding_site = sum(relevant_tf_binding_site, na.rm = TRUE)) %>%
            dplyr::group_by(sample_ID, gene_name) %>%
            dplyr::summarize(n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes = sum(n_relevant_tf_binding_site > 0))
    } else {
        genes_n_relevant_somatic_SNVs_from_gene_TAD_cis_activated_only_table <- data.frame(
            sample_ID = character(),
            gene_name = character(),
            n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes = integer()
        )
    }

    # from genehancer TAD
    genes_n_SVs_from_genehancer_TAD_all_table <- data$SV_genehancer_all_table %>%
        dplyr::group_by(sample_ID, gene_name, sv_break1_chrom, sv_break1_pos, sv_break2_chrom, sv_break2_pos) %>%
        dplyr::summarize(n_genehancers_for_SV_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_SVs_from_genehancer_TAD = dplyr::n())

    genes_n_CNAs_from_genehancer_TAD_all_table <- data$CNA_genehancer_all_table %>%
        dplyr::group_by(sample_ID, gene_name, cna_chrom, cna_start, cna_end) %>%
        dplyr::summarize(n_genehancers_for_CNA_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_CNAs_from_genehancer_TAD = dplyr::n())

    genes_n_somatic_SNVs_from_genehancer_TAD_all_table <- data$somatic_SNV_genehancer_all_table %>%
        dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
        dplyr::summarize(n_genehancers_for_somatic_SNV_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_somatic_SNVs_from_genehancer_TAD = dplyr::n())

    if (has_run_tf_binding_site_analysis) {
        genes_n_relevant_somatic_SNVs_from_genehancer_TAD_cis_activated_only_table <- data$somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos, genehancer_id) %>%
            dplyr::summarize(n_relevant_tf_binding_site = sum(relevant_tf_binding_site, na.rm = TRUE)) %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
            # NEXT LINE: n_relevant_tf_binding_site should be the same for all genehancers of 1 SNV
            dplyr::summarize(n_genehancers_for_SNV_gene_combination = dplyr::n(), n_relevant_tf_binding_site = max(n_relevant_tf_binding_site)) %>% # gets discarded right now later a relevant count can be implemented here
            dplyr::group_by(sample_ID, gene_name) %>%
            dplyr::summarize(n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes = sum(n_relevant_tf_binding_site > 0))
    } else {
        genes_n_relevant_somatic_SNVs_from_genehancer_TAD_cis_activated_only_table <- data.frame(
            sample_ID = character(),
            gene_name = character(),
            n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes = integer()
        )
    }

    # from anywhere
    genes_n_SVs_from_anywhere_all_table <- SV_all_merged_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_SVs = dplyr::n())

    genes_n_CNAs_from_anywhere_all_table <- CNA_all_merged_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_CNAs = dplyr::n())

    genes_n_somatic_SNVs_from_anywhere_all_table <- somatic_SNV_all_merged_table %>%
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_somatic_SNVs = dplyr::n())

    if (has_run_tf_binding_site_analysis) {
        genes_n_relevant_somatic_SNVs_from_anywhere_cis_activated_only_table <- somatic_SNV_tf_binding_data_cis_activated_only_merged_table %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos, genehancer_id) %>%
            dplyr::summarize(n_relevant_tf_binding_site = sum(relevant_tf_binding_site, na.rm = TRUE)) %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
            # NEXT LINE: n_relevant_tf_binding_site should be the same for all genehancers of 1 SNV
            dplyr::summarize(n_genehancers_for_SNV_gene_combination = dplyr::n(), n_relevant_tf_binding_site = max(n_relevant_tf_binding_site)) %>% # gets discarded right now later a relevant count can be implemented here
            dplyr::group_by(sample_ID, gene_name) %>%
            dplyr::summarize(n_relevant_somatic_SNVs_valid_for_cis_activated_genes = sum(n_relevant_tf_binding_site > 0))
    } else {
        genes_n_relevant_somatic_SNVs_from_anywhere_cis_activated_only_table <- data.frame(
            sample_ID = character(),
            gene_name = character(),
            n_relevant_somatic_SNVs_valid_for_cis_activated_genes = integer()
        )
    }

    # with chipseq region affected
    genes_n_SVs_with_chipseq_all_table <- data$SV_chipseq_all_table %>%
        dplyr::group_by(sample_ID, gene_name, sv_break1_chrom, sv_break1_pos, sv_break2_chrom, sv_break2_pos) %>%
        dplyr::summarize(n_chipseqs_for_SV_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_SVs_with_chipseq = dplyr::n())

    genes_n_CNAs_with_chipseq_all_table <- data$CNA_chipseq_all_table %>%
        dplyr::group_by(sample_ID, gene_name, cna_chrom, cna_start, cna_end) %>%
        dplyr::summarize(n_chipseqs_for_CNA_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_CNAs_with_chipseq = dplyr::n())

    genes_n_somatic_SNVs_with_chipseq_all_table <- data$somatic_SNV_chipseq_all_table %>%
        dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
        dplyr::summarize(n_chipseqs_for_somatic_SNV_gene_combination = dplyr::n()) %>% # gets discarded right now later a relevant count can be implemented here
        dplyr::group_by(sample_ID, gene_name) %>%
        dplyr::summarize(n_somatic_SNVs_with_chipseq = dplyr::n())


    if (has_run_tf_binding_site_analysis) {
        genes_n_relevant_somatic_SNVs_with_chipseq_cis_activated_only_table <- somatic_SNV_tf_binding_data_cis_activated_only_merged_table %>%
            dplyr::filter(overlaps_with_chipseq == TRUE) %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos, genehancer_id) %>%
            dplyr::summarize(n_relevant_tf_binding_site = sum(relevant_tf_binding_site, na.rm = TRUE)) %>%
            dplyr::group_by(sample_ID, gene_name, snv_chrom, snv_pos) %>%
            # NEXT LINE: n_relevant_tf_binding_site should be the same for all genehancers of 1 SNV
            dplyr::summarize(n_genehancers_for_SNV_gene_combination = dplyr::n(), n_relevant_tf_binding_site = max(n_relevant_tf_binding_site)) %>% # gets discarded right now later a relevant count can be implemented here
            dplyr::group_by(sample_ID, gene_name) %>%
            dplyr::summarize(n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes = sum(n_relevant_tf_binding_site > 0))
    } else {
        genes_n_relevant_somatic_SNVs_with_chipseq_cis_activated_only_table <- data.frame(
            sample_ID = character(),
            gene_name = character(),
            n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes = integer()
        )
    }

    # add mut counts to cis-activation-summary
    cis_activation_summary_table_mut_info <- cis_activation_summary_table_with_groups %>%
        dplyr::left_join(genes_n_SVs_from_gene_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_CNAs_from_gene_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_somatic_SNVs_from_gene_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_relevant_somatic_SNVs_from_gene_TAD_cis_activated_only_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_SVs_from_genehancer_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_CNAs_from_genehancer_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_somatic_SNVs_from_genehancer_TAD_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_relevant_somatic_SNVs_from_genehancer_TAD_cis_activated_only_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_SVs_from_anywhere_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_CNAs_from_anywhere_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_somatic_SNVs_from_anywhere_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_relevant_somatic_SNVs_from_anywhere_cis_activated_only_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_SVs_with_chipseq_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_CNAs_with_chipseq_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_somatic_SNVs_with_chipseq_all_table, by = c("sample_ID", "gene_name")) %>%
        dplyr::left_join(genes_n_relevant_somatic_SNVs_with_chipseq_cis_activated_only_table, by = c("sample_ID", "gene_name")) %>%
        tidyr::replace_na(
            list(
                n_SVs_from_gene_TAD = 0,
                n_CNAs_from_gene_TAD = 0,
                n_somatic_SNVs_from_gene_TAD = 0,
                n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes = 0,
                n_SVs_from_genehancer_TAD = 0,
                n_CNAs_from_genehancer_TAD = 0,
                n_somatic_SNVs_from_genehancer_TAD = 0,
                n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes = 0,
                n_SVs = 0,
                n_CNAs = 0,
                n_somatic_SNVs = 0,
                n_relevant_somatic_SNVs_valid_for_cis_activated_genes = 0,
                n_SVs_with_chipseq = 0,
                n_CNAs_with_chipseq = 0,
                n_somatic_SNVs_with_chipseq = 0,
                n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes = 0
            )
        ) %>%
        dplyr::mutate(
            n_muts_from_gene_TAD = n_CNAs_from_gene_TAD + n_SVs_from_gene_TAD + n_somatic_SNVs_from_gene_TAD,
            n_muts_from_genehancer_TAD = n_CNAs_from_genehancer_TAD + n_SVs_from_genehancer_TAD + n_somatic_SNVs_from_genehancer_TAD,
            n_muts_with_chipseq = n_CNAs_with_chipseq + n_SVs_with_chipseq + n_somatic_SNVs_with_chipseq,
            n_muts = n_CNAs + n_SVs + n_somatic_SNVs
        )

    # cis activation counts by gene/sample -----------------------------------
    # by sample
    cis_activated_genes_by_sample <- cis_activation_summary_table_mut_info %>%
        # dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::group_by(sample_ID) %>%
        dplyr::summarize(n_cis_activated_genes = sum(cis_activated_gene, na.rm = T))

    cis_activated_genes_by_gene <- cis_activation_summary_table_mut_info %>%
        dplyr::group_by(gene_name, start, end, chrom, gene_type) %>%
        dplyr::summarize(
            n_cis_activated_samples = sum((cis_activated_gene == TRUE), na.rm = T),
            n_CA_plus_SV = sum(cis_activated_gene & (n_SVs > 0), na.rm = T),
            n_CA_plus_CNA = sum(cis_activated_gene & (n_CNAs > 0), na.rm = T),
            n_CA_plus_mut = sum(cis_activated_gene & (n_muts > 0), na.rm = T)
        )

    # convert cis activation summary for upsetr -----------------------------------
    cis_activation_summary_table_for_upsetr <- cis_activation_summary_table_mut_info %>%
        dplyr::mutate(is_imprinted_gene = (!(imprinting_status == "no_imprinting"))) %>%
        dplyr::mutate(
            any_muts = (n_muts > 0),
            any_CNAs = (n_CNAs > 0),
            any_SVs = (n_SVs > 0),
            any_somatic_SNVs = (n_somatic_SNVs > 0),
            any_relevant_somatic_SNVs = (n_relevant_somatic_SNVs_valid_for_cis_activated_genes > 0),
        ) %>%
        dplyr::mutate(passed_any_ASE_filter = (passed_ASE_filter | ASE_filter_rescued_oncogene | ASE_marker_run_overlap)) %>%
        dplyr::mutate_at(
            c("is_cancer_gene", "is_imprinted_gene", "passed_FPKM_filter", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene", "any_muts", "any_CNAs", "any_SVs", "any_somatic_SNVs", "any_relevant_somatic_SNVs"),
            binarify
        ) %>%
        dplyr::mutate(set_id = paste0(sample_ID, "___", gene_name)) %>%
        dplyr::select(c("set_id", "sample_ID", "gene_name", "is_cancer_gene", "is_imprinted_gene", "passed_FPKM_filter", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap", "passed_any_ASE_filter", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene", "any_muts", "any_CNAs", "any_SVs", "any_somatic_SNVs", "any_relevant_somatic_SNVs"))


    # convert cis activation summary for mut analysis -----------------------------------
    cis_activation_summary_table_ANY_mut_info <- cis_activation_summary_table_mut_info %>%
        dplyr::mutate(
            any_muts = (n_muts > 0),
            any_muts_from_gene_TAD = (n_muts_from_gene_TAD > 0),
            any_muts_from_genehancer_TAD = (n_muts_from_genehancer_TAD > 0),
            any_muts_with_chipseq = (n_muts_with_chipseq > 0),
            any_CNAs = (n_CNAs > 0),
            any_CNAs_from_gene_TAD = (n_CNAs_from_gene_TAD > 0),
            any_CNAs_from_genehancer_TAD = (n_CNAs_from_genehancer_TAD > 0),
            any_CNAs_with_chipseq = (n_CNAs_with_chipseq > 0),
            any_SVs = (n_SVs > 0),
            any_SVs_from_gene_TAD = (n_SVs_from_gene_TAD > 0),
            any_SVs_from_genehancer_TAD = (n_SVs_from_genehancer_TAD > 0),
            any_SVs_with_chipseq = (n_SVs_with_chipseq > 0),
            any_somatic_SNVs = (n_somatic_SNVs > 0),
            any_somatic_SNVs_from_gene_TAD = (n_somatic_SNVs_from_gene_TAD > 0),
            any_somatic_SNVs_from_genehancer_TAD = (n_somatic_SNVs_from_genehancer_TAD > 0),
            any_somatic_SNVs_with_chipseq = (n_somatic_SNVs_with_chipseq > 0),
            any_relevant_somatic_SNVs = (n_relevant_somatic_SNVs_valid_for_cis_activated_genes > 0),
            any_relevant_somatic_SNVs_from_gene_TAD = (n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes > 0),
            any_relevant_somatic_SNVs_from_genehancer_TAD = (n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes > 0),
            any_relevant_somatic_SNVs_with_chipseq = (n_relevant_somatic_SNVs_with_chipseq_valid_for_cis_activated_genes > 0),
        ) %>%
        # convert NA -> FALSE
        dplyr::mutate(cis_activated_gene = (cis_activated_gene & (!(is.na(cis_activated_gene)))))


    # recurrent SV TAD combinations -----------------------------------
    # sample_ID, gene_name, cis_activated_gene, sv_break1_chrom, sv_break1_pos, sv_break2_chrom, sv_break2_pos, sv_type, eventInversion
    recurrent_SV_TAD_combinations <- SV_all_merged_table %>%
        # only cross TAD combinations
        dplyr::filter(!(
            (TAD_sv_break1_chrom == TAD_sv_break1_chrom) &
                (TAD_sv_break1_start == TAD_sv_break2_start) &
                (TAD_sv_break1_end == TAD_sv_break2_end))) %>%
        dplyr::group_by(sample_ID, gene_name, cis_activated_gene, TAD_combination) %>%
        dplyr::summarise(n_SVs = dplyr::n()) %>%
        dplyr::group_by(sample_ID, TAD_combination) %>%
        dplyr::summarise(
            affected_genes = paste(gene_name, collapse = ", "),
            n_affected_genes = dplyr::n(),
            affected_cis_activated_genes = paste(gene_name[(!is.na(cis_activated_gene)) & cis_activated_gene], collapse = ", "),
            n_affected_cis_activated_genes = sum(cis_activated_gene, na.rm = TRUE)
        )

    recurrent_SV_TAD_combinations_summary <- recurrent_SV_TAD_combinations %>%
        dplyr::group_by(TAD_combination, affected_genes) %>%
        dplyr::summarise(
            n_samples_with_TAD_combination = dplyr::n(),
            samples_with_TAD_combination = paste(sample_ID, collapse = ", "),
            n_samples_with_TAD_combination_and_cis_activated_genes_affected = sum(affected_cis_activated_genes > 0),
            affected_cis_activated_genes = paste(paste0(sample_ID, ": ", affected_cis_activated_genes)[n_affected_cis_activated_genes > 0], collapse = "; ")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n_samples_with_TAD_combination > 1) %>%
        dplyr::arrange(desc(n_samples_with_TAD_combination_and_cis_activated_genes_affected))



    # RETURN LIST -----------------------------------
    if (has_run_tf_binding_site_analysis) {
        return(list(
            # imported data
            result_paths = data$result_paths,
            marker_data_table = data$marker_data_table,
            expression_data_table = data$expression_data_table,
            somatic_SNV_data_table = data$somatic_SNV_data_table,
            SV_data_table = data$SV_data_table,
            CNA_data_table = data$CNA_data_table,
            cis_activation_summary_table = data$cis_activation_summary_table,
            SV_genes_all_table = data$SV_genes_all_table,
            SV_genehancer_all_table = data$SV_genehancer_all_table,
            SV_chipseq_all_table = data$SV_chipseq_all_table,
            CNA_genes_all_table = data$CNA_genes_all_table,
            CNA_genehancer_all_table = data$CNA_genehancer_all_table,
            CNA_chipseq_all_table = data$CNA_chipseq_all_table,
            somatic_SNV_genes_all_table = data$somatic_SNV_genes_all_table,
            somatic_SNV_genehancer_all_table = data$somatic_SNV_genehancer_all_table,
            somatic_SNV_chipseq_all_table = data$somatic_SNV_chipseq_all_table,
            somatic_SNV_tf_binding_data_genes_cis_activated_only_table = data$somatic_SNV_tf_binding_data_genes_cis_activated_only_table,
            somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table = data$somatic_SNV_tf_binding_data_genehancer_cis_activated_only_table,
            somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table = data$somatic_SNV_tf_binding_data_gene_summary_max_score_diff_table,
            somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table = data$somatic_SNV_tf_binding_data_genehancer_summary_max_score_diff_table,
            somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table = data$somatic_SNV_tf_binding_data_chipseq_summary_max_score_diff_table,
            cis_activation_summary_with_n_relevant_somatic_SNVs_table = data$cis_activation_summary_with_n_relevant_somatic_SNVs_table,
            genehancer_with_cis_activation_summary_table = data$genehancer_with_cis_activation_summary_table,
            chipseq_by_genes_table = data$chipseq_by_genes_table,

            # processed data
            subgroups_table = subgroups_table,
            subgroups = subgroups,
            samples_by_subgroup = samples_by_subgroup,
            SV_all_merged_table = SV_all_merged_table,
            CNA_all_merged_table = CNA_all_merged_table,
            somatic_SNV_all_merged_table = somatic_SNV_all_merged_table,
            cis_activation_summary_table_mut_info = cis_activation_summary_table_mut_info,
            cis_activated_genes_by_sample = cis_activated_genes_by_sample,
            cis_activated_genes_by_gene = cis_activated_genes_by_gene,
            cis_activation_summary_table_for_upsetr = cis_activation_summary_table_for_upsetr,
            cis_activation_summary_table_ANY_mut_info = cis_activation_summary_table_ANY_mut_info,
            recurrent_SV_TAD_combinations = recurrent_SV_TAD_combinations,
            recurrent_SV_TAD_combinations_summary = recurrent_SV_TAD_combinations_summary
        ))
    } else {
        return(list(
            # imported data
            result_paths = data$result_paths,
            marker_data_table = data$marker_data_table,
            expression_data_table = data$expression_data_table,
            somatic_SNV_data_table = data$somatic_SNV_data_table,
            SV_data_table = data$SV_data_table,
            CNA_data_table = data$CNA_data_table,
            cis_activation_summary_table = data$cis_activation_summary_table,
            SV_genes_all_table = data$SV_genes_all_table,
            SV_genehancer_all_table = data$SV_genehancer_all_table,
            SV_chipseq_all_table = data$SV_chipseq_all_table,
            CNA_genes_all_table = data$CNA_genes_all_table,
            CNA_genehancer_all_table = data$CNA_genehancer_all_table,
            CNA_chipseq_all_table = data$CNA_chipseq_all_table,
            somatic_SNV_genes_all_table = data$somatic_SNV_genes_all_table,
            somatic_SNV_genehancer_all_table = data$somatic_SNV_genehancer_all_table,
            somatic_SNV_chipseq_all_table = data$somatic_SNV_chipseq_all_table,
            genehancer_with_cis_activation_summary_table = data$genehancer_with_cis_activation_summary_table,
            chipseq_by_genes_table = data$chipseq_by_genes_table,

            # processed data
            subgroups_table = subgroups_table,
            subgroups = subgroups,
            samples_by_subgroup = samples_by_subgroup,
            SV_all_merged_table = SV_all_merged_table,
            CNA_all_merged_table = CNA_all_merged_table,
            somatic_SNV_all_merged_table = somatic_SNV_all_merged_table,
            cis_activation_summary_table_mut_info = cis_activation_summary_table_mut_info,
            cis_activated_genes_by_sample = cis_activated_genes_by_sample,
            cis_activated_genes_by_gene = cis_activated_genes_by_gene,
            cis_activation_summary_table_for_upsetr = cis_activation_summary_table_for_upsetr,
            cis_activation_summary_table_ANY_mut_info = cis_activation_summary_table_ANY_mut_info,
            recurrent_SV_TAD_combinations = recurrent_SV_TAD_combinations,
            recurrent_SV_TAD_combinations_summary = recurrent_SV_TAD_combinations_summary
        ))
    }
}