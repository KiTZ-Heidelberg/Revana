create_report_gene_subdocument <- function(processed_data, GENE_NAME, color_palette, HTML_report_subdocument_directory, add_custom_HTML_to_gene_subdocument, add_custom_HTML_to_gene_subdocument_sample_data, genome_name = "hg19") {
    HTML_gene_directory <- create_dir_if_missing(file.path(HTML_report_subdocument_directory,GENE_NAME))
    HTML_gene_figure_directory <- create_dir_if_missing(file.path(HTML_gene_directory,"figures"))

    subdocument_path <- file.path(HTML_gene_directory,"gene_overview.html")
    cat(paste0("    Generating plots for GENE: ",GENE_NAME," ...\n"))

    # process data -----------------------------------
    # general data

    cis_activation_summary_table_of_gene <- processed_data$cis_activation_summary_table_mut_info %>%
        dplyr::filter(gene_name == GENE_NAME)

    cis_activated_samples_of_gene <- cis_activation_summary_table_of_gene %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::pull(sample_ID)

    cis_activated_samples_of_gene_with_any_SVs <- cis_activation_summary_table_of_gene %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::filter(n_SVs > 0) %>%
        dplyr::pull(sample_ID)

    cis_activated_samples_of_gene_with_any_CNAs <- cis_activation_summary_table_of_gene %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::filter(n_CNAs > 0) %>%
        dplyr::pull(sample_ID)

    cis_activated_samples_of_gene_with_any_somatic_SNVs <- cis_activation_summary_table_of_gene %>%
        dplyr::filter(cis_activated_gene == TRUE) %>%
        dplyr::filter(n_somatic_SNVs > 0) %>%
        dplyr::pull(sample_ID)

    samples_of_gene_to_display_sample_data <- cis_activation_summary_table_of_gene %>%
        dplyr::filter((cis_activated_gene == TRUE) | (n_SVs > 0) | (n_CNAs > 0)) %>%
        dplyr::arrange(dplyr::desc(cis_activated_gene)) %>% 
        dplyr::pull(sample_ID)

    n_cis_activated_samples_of_gene <- length(cis_activated_samples_of_gene)
    n_cis_activated_samples_of_gene_with_any_SVs <- length(cis_activated_samples_of_gene_with_any_SVs)
    n_cis_activated_samples_of_gene_with_any_CNAs <- length(cis_activated_samples_of_gene_with_any_CNAs)
    n_cis_activated_samples_of_gene_with_any_somatic_SNVs <- length(cis_activated_samples_of_gene_with_any_somatic_SNVs)

    gene_info <- cis_activation_summary_table_of_gene %>%
        dplyr::slice(1)
    

    # mutation subset of GENE_NAME
    SVs_of_gene <- processed_data$SV_genes_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    CNAs_of_gene <- processed_data$CNA_genes_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    somatic_SNVs_of_gene <- processed_data$somatic_SNV_genes_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    somatic_SNVs_of_gene_selected_cols <- somatic_SNVs_of_gene %>%
        dplyr::select(chr = snv_chrom, pos = snv_pos, sample_ID, cis_activated_gene) %>% #
        dplyr::mutate(label = "somatic_SNV")

    SV_genehancer_of_gene <- processed_data$SV_genehancer_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    CNA_genehancer_of_gene <- processed_data$CNA_genehancer_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    somatic_SNV_genehancer_of_gene <- processed_data$somatic_SNV_genehancer_all_table %>%
        dplyr::filter(gene_name == GENE_NAME)

    # somatic_SNV_tf_binding_data_of_gene <- processed_data$somatic_SNV_tf_binding_data_genes_cis_activated_only_table %>%
    #     dplyr::filter(gene_name == GENE_NAME)
    
    SV_left_break <- SVs_of_gene %>%
        dplyr::select(chr = sv_break1_chrom, pos = sv_break1_pos, sample_ID, cis_activated_gene) %>%
        dplyr::mutate(label = "SV")

    SV_right_break <- SVs_of_gene %>%
        dplyr::select(chr = sv_break2_chrom, pos = sv_break2_pos, sample_ID, cis_activated_gene) %>%
        dplyr::mutate(label = "SV")

    SV_all_breaks <- rbind(SV_left_break, SV_right_break) %>%
        dplyr::distinct()

    if(length(somatic_SNVs_of_gene_selected_cols$chr)>0){
        GenomeInfoDb::seqlevelsStyle(somatic_SNVs_of_gene_selected_cols$chr) <- "UCSC"
    }
    if(length(CNAs_of_gene$cna_chrom)>0){
        GenomeInfoDb::seqlevelsStyle(CNAs_of_gene$cna_chrom) <- "UCSC"
    }
    if(length(SV_all_breaks$chr)>0){
        GenomeInfoDb::seqlevelsStyle(SV_all_breaks$chr) <- "UCSC"
    }
    
    SVs_and_somatic_SNVs_of_gene <- rbind(somatic_SNVs_of_gene_selected_cols, SV_all_breaks, fill = T) %>%
        dplyr::filter(chr == gene_info[['chrom']]) %>% # to filter out SV breaks from other chroms -> as SVs_and_somatic_SNVs_of_gene is used in plotting the locus
        dplyr::mutate(pos = as.numeric(pos))


    # TAD boundaries
    min_TAD_start <- gene_info$variant_overlap_window_start
    max_TAD_end <- gene_info$variant_overlap_window_end

    # in case no TAD was matched for gene
        if(is.na(min_TAD_start)|is.na(min_TAD_start)){
            min_TAD_start <- gene_info[["start"]] - 1000000
            max_TAD_end <- gene_info[["end"]] + 1000000
        }

    # genes within TAD
    cis_activated_genes_by_gene_within_TAD <- processed_data$cis_activated_genes_by_gene %>% 
        dplyr::filter(
            (start <= !!max_TAD_end) &
            (end >= !!min_TAD_start) &
            (chrom == gene_info[["chrom"]])
        ) %>%
        dplyr::mutate(type = "cds")

    genes_within_TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        cis_activated_genes_by_gene_within_TAD,
        keep.extra.columns=TRUE
        )

    # chipseq within TAD
    chipseq_by_genes_within_TAD <- processed_data$chipseq_by_genes_table %>%
        dplyr::filter(gene_name == GENE_NAME) %>%
        dplyr::filter(
                (start <= !!max_TAD_end) &
                (end >= !!min_TAD_start) &
                (chrom == gene_info[["chrom"]])
            )

    # affected chipseq regions
    SV_chipseq_of_gene <- processed_data$SV_chipseq_all_table[gene_name==GENE_NAME]
    CNA_chipseq_of_gene <- processed_data$CNA_chipseq_all_table[gene_name==GENE_NAME]
    somatic_SNV_chipseq_of_gene <- processed_data$somatic_SNV_chipseq_all_table[gene_name==GENE_NAME]

    SV_chipseq_of_gene[, chipseq_id:=paste0(chipseq_chrom ,":" ,chipseq_start ,"-" ,chipseq_end)]
    CNA_chipseq_of_gene[, chipseq_id:=paste0(chipseq_chrom ,":" ,chipseq_start ,"-" ,chipseq_end)]
    somatic_SNV_chipseq_of_gene[, chipseq_id:=paste0(chipseq_chrom ,":" ,chipseq_start ,"-" ,chipseq_end)]

    SV_chipseq_of_gene_only_chipseq_columns <- SV_chipseq_of_gene[,c("chipseq_chrom","chipseq_start","chipseq_end", "chipseq_id")]
    CNA_chipseq_of_gene_only_chipseq_columns <- CNA_chipseq_of_gene[,c("chipseq_chrom","chipseq_start","chipseq_end", "chipseq_id")]
    somatic_SNV_chipseq_of_gene_only_chipseq_columns <- somatic_SNV_chipseq_of_gene[,c("chipseq_chrom","chipseq_start","chipseq_end", "chipseq_id")]

    affected_chipseq_regions_of_gene <- unique(rbind(SV_chipseq_of_gene_only_chipseq_columns ,CNA_chipseq_of_gene_only_chipseq_columns ,somatic_SNV_chipseq_of_gene_only_chipseq_columns ))


    # genehancer subset of gene_name
    genehancer_with_cis_activation_summary_of_gene <- processed_data$genehancer_with_cis_activation_summary_table %>%
        dplyr::filter(connected_gene == GENE_NAME)

    genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer <- genehancer_with_cis_activation_summary_of_gene %>%
        dplyr::group_by(genehancer_id) %>% 
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()

    genehancer_with_cis_activation_summary_of_gene_within_TAD <- genehancer_with_cis_activation_summary_of_gene %>%
        dplyr::filter(
            (start_genehancer <= !!max_TAD_end) &
            (end_genehancer >= !!min_TAD_start) &
            (chrom_genehancer == gene_info[["chrom"]])
        )

    # create plots -----------------------------------
    # mutation_type matrix
        mutation_type_matrix_plot_data <- cis_activation_summary_table_of_gene %>%
            dplyr::filter((cis_activated_gene == TRUE)| (n_muts>0) | (cis_activated_gene_copy_number_normalized == TRUE)) %>%
            dplyr::mutate(cis_activation = ifelse(cis_activated_gene == TRUE,1,0)) %>%
            dplyr::mutate(cis_activation_copy_number_normalized = ifelse(cis_activated_gene_copy_number_normalized == TRUE,1,0)) %>%
            dplyr::mutate(sample_ID = forcats::fct_reorder(sample_ID, cis_activation)) %>%
            tidyr::gather(n_somatic_SNVs, n_relevant_somatic_SNVs_valid_for_cis_activated_genes, n_SVs, n_CNAs, cis_activation, cis_activation_copy_number_normalized, key="mutation_type",value="n_of_detections", factor_key = T) %>%
            dplyr::filter(n_of_detections>0)

        if(nrow(mutation_type_matrix_plot_data) > 0){
             mutation_type_matrix_plot <- ggplot2::ggplot(mutation_type_matrix_plot_data) +
            ggplot2::geom_point(ggplot2::aes(y=sample_ID, x = mutation_type, color = subgroup), size=5) +
            ggplot2::geom_text(ggplot2::aes(y=sample_ID, x = mutation_type, label=ifelse(!(mutation_type %in% c("cis_activation", "cis_activation_copy_number_normalized")),n_of_detections,"")), color = "white") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) + 
            color_palette$scale_color_subgroup +
            ggplot2::scale_x_discrete(breaks=c("cis_activation","cis_activation_copy_number_normalized", "n_somatic_SNVs", "n_relevant_somatic_SNVs_valid_for_cis_activated_genes", "n_SVs", "n_CNAs"), labels = c("cis-activated", "cis-act.\n(CN norm.)", "SNV/Indels", "relevant SNV/Indels","SVs", "CNAs"),position = "top", name="Mutation Type", drop=F)
        }else{
            mutation_type_matrix_plot <-  ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0,0,label='No mutated or cis-activated samples')) +
                ggplot2::xlab(NULL)
        }
    

        mutation_type_matrix_plot_path <- file.path(HTML_gene_figure_directory, paste0("mutation_type_matrix_plot___",GENE_NAME,".svg"))
        save_ggplot(mutation_type_matrix_plot, mutation_type_matrix_plot_path)


        # passed_filter UpSetR all samples
        passed_filter_all_samples_upsetr_path <- file.path(HTML_gene_figure_directory, paste0("passed_filter_all_samples_upsetr___",GENE_NAME,".svg"))
        passed_filter_all_samples_upsetr_temp_df <- processed_data$cis_activation_summary_table_for_upsetr %>%
            dplyr::filter(gene_name == GENE_NAME) %>%
            # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
            dplyr::mutate(ones1 = 1, ones2 = 1) %>%
            dplyr::select("set_id","sample_ID","ones1","is_cancer_gene", "is_imprinted_gene", "passed_FPKM_filter", "passed_ASE_filter", "ASE_filter_rescued_oncogene", "ASE_marker_run_overlap", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene","ones2")
        
        are_all_sets_empty_in_passed_filter_all_samples_upsetr_data <- passed_filter_all_samples_upsetr_temp_df %>%
            { sum(c(
                .$is_cancer_gene,
                .$is_imprinted_gene,
                .$passed_FPKM_filter,
                .$passed_ASE_filter,
                .$ASE_filter_rescued_oncogene,
                .$ASE_marker_run_overlap,
                .$passed_mean_delta_abs_filter,
                .$passed_OHE_filter,
                .$cis_activated_gene))} %>%
                {. == 0 }

        if(are_all_sets_empty_in_passed_filter_all_samples_upsetr_data == TRUE){
            passed_filter_all_samples_upsetr_plot <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0,0,label='No Upset block available\nas no sample passed any filter')) +
                ggplot2::xlab(NULL)
        }else{
            passed_filter_all_samples_upsetr_plot <- UpSetR::upset(
                passed_filter_all_samples_upsetr_temp_df,
                sets = c("is_cancer_gene", "is_imprinted_gene", "passed_FPKM_filter", "passed_ASE_filter", "ASE_filter_rescued_oncogene","ASE_marker_run_overlap", "passed_mean_delta_abs_filter", "passed_OHE_filter", "cis_activated_gene")
            )
        }

        
        # cat("Saving UpSetR image\n")
        svglite::svglite(passed_filter_all_samples_upsetr_path)
        print(passed_filter_all_samples_upsetr_plot)
        dev.off()

        # mutation_type UpSetR all samples
        mutation_type_all_samples_upsetr_path <- file.path(HTML_gene_figure_directory, paste0("mutation_type_all_samples_upsetr___",GENE_NAME,".svg"))
        # CHANGE THIS LINE TO SELF GENERATED MUT COUNT DATA
        mutation_type_all_samples_upsetr_temp_df <- processed_data$cis_activation_summary_table_for_upsetr %>%
            dplyr::filter(gene_name == GENE_NAME) %>%
            # FIX FOR: currently UpSetR bug of not including first column of dataset if only made up of 0s
            dplyr::mutate(ones1 = 1, ones2 = 1) %>%
            dplyr::select("set_id","sample_ID","ones1","cis_activated_gene","any_muts", "any_SVs", "any_CNAs","any_somatic_SNVs", "any_relevant_somatic_SNVs","ones2")
        
        are_all_sets_empty_in_mutation_type_all_samples_upsetr_data <- mutation_type_all_samples_upsetr_temp_df %>%
            { sum(c(
                .$cis_activated_gene,
                .$any_muts,
                .$any_SVs,
                .$any_CNAs,
                .$any_somatic_SNVs,
                .$any_relevant_somatic_SNVs
            ))} %>% { . == 0 }

        if(are_all_sets_empty_in_mutation_type_all_samples_upsetr_data == TRUE){
            mutation_type_all_samples_upsetr_plot <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0,0,label='No Upset block available\nas all no sample passed any filter')) +
                ggplot2::xlab(NULL)
        }else{
            mutation_type_all_samples_upsetr_plot <- UpSetR::upset(
                mutation_type_all_samples_upsetr_temp_df,
                sets = c("cis_activated_gene", "any_muts", "any_SVs", "any_CNAs","any_somatic_SNVs", "any_relevant_somatic_SNVs")
            )
        }

        # cat("Saving UpSetR image\n")
        svglite::svglite(mutation_type_all_samples_upsetr_path)
        print(mutation_type_all_samples_upsetr_plot)
        dev.off()

        # Expression by sample
        expression_by_sample_path <- file.path(HTML_gene_figure_directory, paste0("expression_by_sample_plot___",GENE_NAME,".svg"))
        expression_by_sample_plot <- cis_activation_summary_table_of_gene %>%
            dplyr::mutate(sample_ID = forcats::fct_reorder(sample_ID, FPKM)) %>%
            ggplot2::ggplot() +
                ggplot2::geom_bar(ggplot2::aes(x=sample_ID, y = FPKM, fill = subgroup, alpha = as.factor(cis_activated_gene)), stat="identity") +
                ggplot2::geom_text(ggplot2::aes(x=sample_ID, y = FPKM, label=FPKM), color = "white", position = ggplot2::position_stack(vjust = 0.5), size = 1.7) +
                # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0) + 
                color_palette$scale_fill_subgroup +
                ggplot2::scale_alpha_manual(values = c(1,0.5), breaks = c(TRUE,FALSE), name = "cis-activated", labels =c("yes", "no"), na.value=0.5) +
                ggplot2::scale_y_discrete(position = "top", name="FPKM") +
                ggplot2::coord_flip() +
                ggplot2::labs(y="Sample ID",x="Expression (FPKM)") +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 2))

        save_ggplot(expression_by_sample_plot, expression_by_sample_path)

        # Expression by group
        expression_by_group_path <- file.path(HTML_gene_figure_directory, paste0("expression_by_group_plot___",GENE_NAME,".svg"))
        expression_by_group_plot <- cis_activation_summary_table_of_gene %>%

            # to adress warning: Removed X rows containing non-finite values (`stat_boxplot()`).
            tidyr::drop_na(FPKM) %>%

            ggplot2::ggplot() +
                ggplot2::geom_boxplot(ggplot2::aes(x=subgroup, y = FPKM)) +
                ggplot2::geom_point(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=FPKM, x = subgroup, color = subgroup), size=3) +
                ggrepel::geom_text_repel(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=FPKM, x = subgroup, label = sample_ID), hjust = 0, nudge_x = 0.05, size=2) +
                ggplot2::geom_label((cis_activation_summary_table_of_gene%>%tidyr::drop_na(FPKM)%>%dplyr::group_by(subgroup)%>%dplyr::summarise(n=dplyr::n(), mean = mean(FPKM))), mapping = ggplot2::aes(x=subgroup, y=mean, label =paste0("n=",n))) +
                color_palette$scale_color_subgroup +
                ggplot2::labs(x="Subgroup",y="Expression (FPKM)")

        save_ggplot(expression_by_group_plot, expression_by_group_path)

        # Coverage Ratio by sample
        cov_ratio_by_sample_path <- file.path(HTML_gene_figure_directory, paste0("cov_ratio_by_sample_plot___",GENE_NAME,".svg"))
        cov_ratio_by_sample_plot <- cis_activation_summary_table_of_gene %>%
            dplyr::mutate(sample_ID = forcats::fct_reorder(sample_ID, FPKM)) %>%
            ggplot2::ggplot() +
                ggplot2::geom_bar(ggplot2::aes(x=sample_ID, y = avg_cov_ratio, fill = subgroup, alpha = as.factor(cis_activated_gene)), stat="identity") +
                ggplot2::geom_text(ggplot2::aes(x=sample_ID, y = avg_cov_ratio, label=avg_cov_ratio), color = "white", position = ggplot2::position_stack(vjust = 0.5), size = 1.7) +
                # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0) + 
                color_palette$scale_fill_subgroup +
                ggplot2::scale_alpha_manual(values = c(1,0.5), breaks = c(TRUE,FALSE), name = "cis-activated", labels =c("yes", "no"), na.value=0.5) +
                ggplot2::scale_y_discrete(position = "top", name="Cov Ratio") +
                ggplot2::coord_flip() +
                ggplot2::labs(y="Sample ID",x="Cov Ratio") +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 2))

        save_ggplot(cov_ratio_by_sample_plot, cov_ratio_by_sample_path)

        # Expression by group - copy_number_normalized
        expression_copy_number_normalized_by_group_path <- file.path(HTML_gene_figure_directory, paste0("expression_copy_number_normalized_by_group_plot___",GENE_NAME,".svg"))
        expression_copy_number_normalized_by_group_plot <- cis_activation_summary_table_of_gene %>%

            # to adress warning: Removed X rows containing non-finite values (`stat_boxplot()`).
            tidyr::drop_na(FPKM_copy_number_normalized) %>%

            ggplot2::ggplot() +
                ggplot2::geom_boxplot(ggplot2::aes(x=subgroup, y = FPKM_copy_number_normalized)) +
                ggplot2::geom_point(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=FPKM_copy_number_normalized, x = subgroup, color = subgroup), size=3) +
                ggrepel::geom_text_repel(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=FPKM_copy_number_normalized, x = subgroup, label = sample_ID), hjust = 0, nudge_x = 0.05, size=2) +
                ggplot2::geom_label((cis_activation_summary_table_of_gene%>%tidyr::drop_na(FPKM_copy_number_normalized)%>%dplyr::group_by(subgroup)%>%dplyr::summarise(n=dplyr::n(), mean = mean(FPKM_copy_number_normalized))), mapping = ggplot2::aes(x=subgroup, y=mean, label =paste0("n=",n))) +
                color_palette$scale_color_subgroup +
                ggplot2::labs(x="Subgroup",y="Expression (FPKM) - copy number normalized")

        save_ggplot(expression_copy_number_normalized_by_group_plot, expression_copy_number_normalized_by_group_path)


        # ASE detection type
        if(n_cis_activated_samples_of_gene > 0) {
            ASE_detection_type_matrix_plot <- cis_activation_summary_table_of_gene %>%
                dplyr::filter((cis_activated_gene == TRUE)) %>%
                dplyr::mutate(passed_ASE_filter_binary = ifelse(passed_ASE_filter == TRUE,1,0)) %>%
                dplyr::mutate(ASE_filter_rescued_oncogene_binary = ifelse(ASE_filter_rescued_oncogene == TRUE,1,0)) %>%
                dplyr::mutate(ASE_marker_run_overlap_binary = ifelse(ASE_marker_run_overlap == TRUE,1,0)) %>%
                tidyr::gather(passed_ASE_filter_binary, ASE_filter_rescued_oncogene_binary, ASE_marker_run_overlap_binary, key="ase_detection_type",value="detected_0_or_1", factor_key = T) %>%
                dplyr::filter(detected_0_or_1>0) %>%
                ggplot2::ggplot() +
                    ggplot2::geom_point(ggplot2::aes(y=sample_ID, x = ase_detection_type, color = subgroup), size=5) +
                    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) + 
                    color_palette$scale_color_subgroup +
                    ggplot2::scale_x_discrete(breaks=c("passed_ASE_filter_binary", "ASE_filter_rescued_oncogene_binary", "ASE_marker_run_overlap_binary"), labels = c("normal", "rescued oncogene", "marker run overlap"),position = "top", name="ASE Detection Type", drop=F)

            ASE_detection_type_matrix_plot_path <- file.path(HTML_gene_figure_directory, paste0("ASE_detection_type_matrix_plot___",GENE_NAME,".svg"))
            save_ggplot(ASE_detection_type_matrix_plot, ASE_detection_type_matrix_plot_path)
        }


        # Allelic imbalance by group
        allelic_imbalance_by_group_path <- file.path(HTML_gene_figure_directory, paste0("allelic_imbalance_by_group_plot___",GENE_NAME,".svg"))
        
        has_any_sample_mean_delta_abs_values <- sum(!is.na(cis_activation_summary_table_of_gene$mean_delta_abs)) > 0
        
        if(has_any_sample_mean_delta_abs_values){
            allelic_imbalance_by_group_plot <- cis_activation_summary_table_of_gene %>%

            # to adress warning: Removed X rows containing non-finite values (`stat_boxplot()`).
            tidyr::drop_na(mean_delta_abs) %>%

            ggplot2::ggplot() +
                ggplot2::geom_boxplot(ggplot2::aes(y=mean_delta_abs, x=subgroup)) +
                ggplot2::geom_point(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=mean_delta_abs, x = subgroup, color = subgroup), size=3) +
                ggrepel::geom_text_repel(dplyr::filter(cis_activation_summary_table_of_gene, cis_activated_gene == TRUE), mapping = ggplot2::aes(y=mean_delta_abs, x = subgroup, label = sample_ID), hjust = 0, nudge_x = 0.05, size=2) +
                ggplot2::geom_label((cis_activation_summary_table_of_gene%>%tidyr::drop_na(mean_delta_abs)%>%dplyr::group_by(subgroup)%>%dplyr::summarise(n=dplyr::n(), mean = median(mean_delta_abs))), mapping = ggplot2::aes(x=subgroup, y=mean, label =paste0("n=",n))) +
                color_palette$scale_color_subgroup +
                ggplot2::labs(y="Allelic Imbalance (mean_delta_abs)",x="Subgroup")
        }else{
            allelic_imbalance_by_group_plot <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0,0,label='Allelic imbalance could not be calculated for any sample\ne.g. due to lack of SNP markers')) +
                ggplot2::xlab(NULL)
        }
        
        save_ggplot(allelic_imbalance_by_group_plot, allelic_imbalance_by_group_path)

        # Allelic imbalance & Expression Dot Plot
        allelic_imbalance_expression_dot_plot_path <- file.path(HTML_gene_figure_directory, paste0("allelic_imbalance_expression_dot_plot___",GENE_NAME,".svg"))
        allelic_imbalance_expression_dot_plot <- cis_activation_summary_table_of_gene %>%
        tidyr::replace_na(list(mean_delta_abs= 0)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(y=FPKM, x = mean_delta_abs, color = subgroup, shape=as.factor(ifelse(cis_activated_gene & (!is.na(cis_activated_gene)),ifelse((n_SVs > 0) | (n_CNAs > 0), 5, ifelse(n_muts > 0, 4,3)),ifelse((n_SVs > 0) | (n_CNAs > 0), 2, ifelse(n_muts > 0, 1,0)))))) +
        ggrepel::geom_text_repel(dplyr::filter(tidyr::replace_na(cis_activation_summary_table_of_gene, list(mean_delta_abs= 0)), cis_activated_gene == TRUE), mapping = ggplot2::aes(y=FPKM, x = mean_delta_abs, label = sample_ID), size=2) + 
        color_palette$scale_color_subgroup +
        ggplot2::scale_shape_manual(values = c(1,0,2,19,15,17), breaks = c(0,1,2,3,4,5), labels=c("not cis-activated","not c.a. with detected somatic SNV","not c.a. with SV/CNA","cis-activated", "c.a. with detected somatic SNV", "c.a. with detected SV/CNA"), name = "Cis-activation") +
        ggplot2::labs(y="Expression (FPKM)",x="Assumed Allelic Imbalance (mean_delta_abs)*", caption = "*assumed that AI = 0 for genes without any SNP markers")

        save_ggplot(allelic_imbalance_expression_dot_plot, allelic_imbalance_expression_dot_plot_path)

        # Allelic imbalance & Expression Dot Plot - SV/CNA highlighted
        allelic_imbalance_expression_dot_plot_SV_CNA_highlighted_path <- file.path(HTML_gene_figure_directory, paste0("allelic_imbalance_expression_dot_plot___",GENE_NAME,"_SV_CNA_highlighted.svg"))
        allelic_imbalance_expression_dot_plot_SV_CNA_highlighted <- cis_activation_summary_table_of_gene %>%
        tidyr::replace_na(list(mean_delta_abs= 0)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(y=FPKM, x = mean_delta_abs, color = subgroup, shape=as.factor(ifelse(cis_activated_gene & (!is.na(cis_activated_gene)),ifelse((n_SVs > 0) | (n_CNAs > 0), 5, ifelse(n_muts > 0, 4,3)),ifelse((n_SVs > 0) | (n_CNAs > 0), 2, ifelse(n_muts > 0, 1,0)))), alpha=as.factor(ifelse((n_SVs > 0) | (n_CNAs > 0), 1, 0)))) +
        ggrepel::geom_text_repel(dplyr::filter(tidyr::replace_na(cis_activation_summary_table_of_gene, list(mean_delta_abs= 0)), (n_SVs > 0) | (n_CNAs > 0)), mapping = ggplot2::aes(y=FPKM, x = mean_delta_abs, label = sample_ID), size=2) + 
        color_palette$scale_color_subgroup +
        ggplot2::scale_alpha_manual(values = c(0.3,1), breaks = c(0,1), labels=c("no SVs/CNAs","with SVs/CNAs"), name = "SVs/CNAs") +
        ggplot2::scale_shape_manual(values = c(1,0,2,19,15,17), breaks = c(0,1,2,3,4,5), labels=c("not cis-activated","not c.a. with detected somatic SNV","not c.a. with SV/CNA","cis-activated", "c.a. with detected somatic SNV", "c.a. with detected SV/CNA"), name = "Cis-activation") +
        ggplot2::labs(y="Expression (FPKM)",x="Assumed Allelic Imbalance (mean_delta_abs)*", caption = "*assumed that AI = 0 for genes without any SNP markers")

        save_ggplot(allelic_imbalance_expression_dot_plot_SV_CNA_highlighted, allelic_imbalance_expression_dot_plot_SV_CNA_highlighted_path)



        # Allelic imbalance & Expression Dot Plot copy_number_normalized
        allelic_imbalance_expression_copy_number_normalized_dot_plot_path <- file.path(HTML_gene_figure_directory, paste0("allelic_imbalance_expression_copy_number_normalized_dot_plot___",GENE_NAME,".svg"))
        allelic_imbalance_expression_copy_number_normalized_dot_plot <- cis_activation_summary_table_of_gene %>%
        tidyr::replace_na(list(mean_delta_abs= 0)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(y=FPKM_copy_number_normalized, x = mean_delta_abs, color = subgroup, shape=as.factor(ifelse(cis_activated_gene_copy_number_normalized & (!is.na(cis_activated_gene_copy_number_normalized)),ifelse((n_SVs > 0) | (n_CNAs > 0), 5, ifelse(n_muts > 0, 4,3)),ifelse((n_SVs > 0) | (n_CNAs > 0), 2, ifelse(n_muts > 0, 1,0)))))) +
        ggrepel::geom_text_repel(dplyr::filter(tidyr::replace_na(cis_activation_summary_table_of_gene, list(mean_delta_abs= 0)), cis_activated_gene_copy_number_normalized == TRUE), mapping = ggplot2::aes(y=FPKM_copy_number_normalized, x = mean_delta_abs, label = sample_ID), size=2) + 
        color_palette$scale_color_subgroup +
        ggplot2::scale_shape_manual(values = c(1,0,2,19,15,17), breaks = c(0,1,2,3,4,5), labels=c("not cis-activated","not c.a. with detected somatic SNV","not c.a. with SV/CNA","cis-activated", "c.a. with detected somatic SNV", "c.a. with detected SV/CNA"), name = "Cis-activation") +
        ggplot2::labs(y="Expression (FPKM) - copy number normalized",x="Assumed Allelic Imbalance (mean_delta_abs)*", caption = "*assumed that AI = 0 for genes without any SNP markers")

        save_ggplot(allelic_imbalance_expression_copy_number_normalized_dot_plot, allelic_imbalance_expression_copy_number_normalized_dot_plot_path)

        # Allelic imbalance & Expression Dot Plot - SV/CNA highlighted - copy_number_normalized
        allelic_imbalance_expression_copy_number_normalized_dot_plot_SV_CNA_highlighted_path <- file.path(HTML_gene_figure_directory, paste0("allelic_imbalance_expression_copy_number_normalized_dot_plot___",GENE_NAME,"_SV_CNA_highlighted.svg"))
        allelic_imbalance_expression_copy_number_normalized_dot_plot_SV_CNA_highlighted <- cis_activation_summary_table_of_gene %>%
        tidyr::replace_na(list(mean_delta_abs= 0)) %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(y=FPKM_copy_number_normalized, x = mean_delta_abs, color = subgroup, shape=as.factor(ifelse(cis_activated_gene_copy_number_normalized & (!is.na(cis_activated_gene_copy_number_normalized)),ifelse((n_SVs > 0) | (n_CNAs > 0), 5, ifelse(n_muts > 0, 4,3)),ifelse((n_SVs > 0) | (n_CNAs > 0), 2, ifelse(n_muts > 0, 1,0)))), alpha=as.factor(ifelse((n_SVs > 0) | (n_CNAs > 0), 1, 0)))) +
        ggrepel::geom_text_repel(dplyr::filter(tidyr::replace_na(cis_activation_summary_table_of_gene, list(mean_delta_abs= 0)), (n_SVs > 0) | (n_CNAs > 0)), mapping = ggplot2::aes(y=FPKM_copy_number_normalized, x = mean_delta_abs, label = sample_ID), size=2) + 
        color_palette$scale_color_subgroup +
        ggplot2::scale_alpha_manual(values = c(0.3,1), breaks = c(0,1), labels=c("no SVs/CNAs","with SVs/CNAs"), name = "SVs/CNAs") +
        ggplot2::scale_shape_manual(values = c(1,0,2,19,15,17), breaks = c(0,1,2,3,4,5), labels=c("not cis-activated","not c.a. with detected somatic SNV","not c.a. with SV/CNA","cis-activated", "c.a. with detected somatic SNV", "c.a. with detected SV/CNA"), name = "Cis-activation") +
        ggplot2::labs(y="Expression (FPKM) - copy number normalized",x="Assumed Allelic Imbalance (mean_delta_abs)*", caption = "*assumed that AI = 0 for genes without any SNP markers")

        save_ggplot(allelic_imbalance_expression_copy_number_normalized_dot_plot_SV_CNA_highlighted, allelic_imbalance_expression_copy_number_normalized_dot_plot_SV_CNA_highlighted_path)

        # gene locus plot
        create_gene_locus_plot <- function (
            min_TAD_start,
            max_TAD_end,
            chipseq_by_genes_within_TAD,
            genehancer_with_cis_activation_summary_of_gene_within_TAD,
            SVs_and_somatic_SNVs_of_gene,
            CNAs_of_gene,
            scale_color_mutation_type,
            cis_activated_genes_by_gene_within_TAD,
            genome_name = "hg19"
        ){
            # silence MESSAGE: packageStartupMessage in packageStartupMessage(msg, domain = NA): Registered S3 method overwritten by 'GGally'
            # silence MESSAGE: Constructing Graphics...
            suppress_certain_message({
                TAD_range <- GenomicRanges::GRanges(seqnames = gene_info[["chrom"]], ranges = IRanges::IRanges(start = min_TAD_start, end = max_TAD_end))
                p_ideo <- ggbio::Ideogram(genome = genome_name, xlabel = TRUE, zoom.region = c(GenomicRanges::start(TAD_range),GenomicRanges::end(TAD_range)), subchr = gene_info[["chrom"]])

                if(nrow(chipseq_by_genes_within_TAD) > 0){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale
                    p_chipseq_y_numeric <- as.numeric(as.factor(chipseq_by_genes_within_TAD$subgroup))
                    p_chipseq_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(chipseq_by_genes_within_TAD$subgroup)))))
                    p_chipseq_y_labels <- levels(as.factor(chipseq_by_genes_within_TAD$subgroup))

                    p_chipseq <- chipseq_by_genes_within_TAD %>%
                        ggplot2::ggplot() +
                        ggplot2::geom_errorbar(mapping = ggplot2::aes(xmin = start, xmax = end, y = p_chipseq_y_numeric), color = "red") +
                        ggplot2::scale_y_discrete(breaks = p_chipseq_y_breaks, labels = p_chipseq_y_labels) +
                        ggplot2::labs(y="")
                }else{
                    p_chipseq <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                if(nrow(genehancer_with_cis_activation_summary_of_gene_within_TAD) > 0){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale
                    p_genehancer_y_numeric <- as.numeric(as.factor(genehancer_with_cis_activation_summary_of_gene_within_TAD$connected_gene))
                    p_genehancer_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(genehancer_with_cis_activation_summary_of_gene_within_TAD$connected_gene)))))
                    p_genehancer_y_labels <- levels(as.factor(genehancer_with_cis_activation_summary_of_gene_within_TAD$connected_gene))

                    p_genehancer <- genehancer_with_cis_activation_summary_of_gene_within_TAD %>%
                    ggplot2::ggplot() +
                    ggplot2::geom_errorbar(mapping = ggplot2::aes(xmin = start_genehancer, xmax = end_genehancer, y = p_genehancer_y_numeric), color = "red") +
                    ggplot2::scale_y_discrete(breaks = p_genehancer_y_breaks, labels = p_genehancer_y_labels) +
                    ggplot2::labs(y="") # +
                    # ggrepel::geom_label_repel(ggplot2::aes(x=start_genehancer, y=as.factor(connected_gene), label=paste0(genehancer_id, " - ", feature_name)), size = 2)

                }else{
                    p_genehancer <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                
                SVs_cis_activated_of_gene <- SVs_and_somatic_SNVs_of_gene %>%
                    dplyr::filter(cis_activated_gene == TRUE) %>%
                    dplyr::filter(label == "SV")

                CNAs_cis_activated_of_gene <- CNAs_of_gene %>% 
                    dplyr::filter(cis_activated_gene == TRUE, cna_chrom == gene_info[["chrom"]])

                if((nrow(SVs_cis_activated_of_gene) > 0)|(nrow(CNAs_cis_activated_of_gene) > 0)){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale
                    CNAs_and_SVs_y <- c(CNAs_cis_activated_of_gene$sample_ID, SVs_cis_activated_of_gene$sample_ID)

                    CNAs_y_numeric <- as.numeric(as.factor(c(CNAs_and_SVs_y)))[0:length(CNAs_cis_activated_of_gene$sample_ID)]
                    SVs_y_numeric <- {if(nrow(SVs_cis_activated_of_gene) > 0) as.numeric(as.factor(c(CNAs_and_SVs_y)))[(length(CNAs_cis_activated_of_gene$sample_ID)+1):length(CNAs_and_SVs_y)] else numeric(0)}
                    p_SVs_and_CNAs_cis_activated_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(CNAs_and_SVs_y)))))
                    p_SVs_and_CNAs_cis_activated_y_labels <- levels(as.factor(CNAs_and_SVs_y))

                    p_SVs_and_CNAs_cis_activated <- SVs_cis_activated_of_gene %>%
                        ggplot2::ggplot() +
                        ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                        ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                        ggplot2::geom_point(ggplot2::aes(x=pos, y = SVs_y_numeric, color=as.factor(label)), show.legend = F) +
                        ggplot2::geom_errorbar(CNAs_cis_activated_of_gene, mapping = ggplot2::aes(xmin = cna_start, xmax = cna_end, y = CNAs_y_numeric), color = "green", width=0.5) + 
                        scale_color_mutation_type +
                        ggplot2::scale_y_continuous(breaks = p_SVs_and_CNAs_cis_activated_y_breaks, labels = p_SVs_and_CNAs_cis_activated_y_labels) +
                        ggplot2::labs(y="")+
                        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
                }else{
                    p_SVs_and_CNAs_cis_activated <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                # p_SVs_and_CNAs_cis_activated <- SVs_and_somatic_SNVs_of_gene %>%
                #     dplyr::filter(cis_activated_gene == TRUE) %>%
                #     dplyr::filter(label == "SV") %>%
                #     ggplot2::ggplot() +
                #     ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                #     ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                #     ggplot2::geom_point(ggplot2::aes(x=pos, y = as.factor(sample_ID), color=as.factor(label)), show.legend = F) +
                #     ggplot2::geom_errorbar(dplyr::filter(CNAs_of_gene, cis_activated_gene == TRUE, cna_chrom == gene_info[["chrom"]]), mapping = ggplot2::aes(xmin = cna_start, xmax = cna_end, y = as.factor(sample_ID)), color = "green") + 
                #     scale_color_mutation_type +
                #     ggplot2::labs(y="")

                SNVs_cis_activated_of_gene <- SVs_and_somatic_SNVs_of_gene %>%
                    dplyr::filter(cis_activated_gene == TRUE) %>%
                    dplyr::filter(label == "somatic_SNV")

                if(nrow(SNVs_cis_activated_of_gene) > 0){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale

                    p_SNVs_cis_activated_y_numeric <- as.numeric(as.factor(SNVs_cis_activated_of_gene$sample_ID))
                    p_SNVs_cis_activated_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(SNVs_cis_activated_of_gene$sample_ID)))))
                    p_SNVs_cis_activated_y_labels <- levels(as.factor(SNVs_cis_activated_of_gene$sample_ID))

                    p_SNVs_cis_activated <- SNVs_cis_activated_of_gene %>%
                        ggplot2::ggplot() +
                        ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                        ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                        { if(length(p_SNVs_cis_activated_y_numeric)>0)ggplot2::geom_point(ggplot2::aes(x=pos, y = p_SNVs_cis_activated_y_numeric, color=as.factor(label)), show.legend = F) } +
                        scale_color_mutation_type +
                        ggplot2::scale_y_continuous(breaks = p_SNVs_cis_activated_y_breaks, labels = p_SNVs_cis_activated_y_labels) +
                        ggplot2::labs(y="")+
                        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

                }else{
                    p_SNVs_cis_activated <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                # p_SNVs_cis_activated <- SVs_and_somatic_SNVs_of_gene %>%
                #     dplyr::filter(cis_activated_gene == TRUE) %>%
                #     dplyr::filter(label == "somatic_SNV") %>%
                #     ggplot2::ggplot() +
                #     ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                #     ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                #     ggplot2::geom_point(ggplot2::aes(x=pos, y = as.factor(sample_ID), color=as.factor(label)), show.legend = F) +
                #     scale_color_mutation_type +
                #     ggplot2::labs(y="")

                SVs_others_of_gene <- SVs_and_somatic_SNVs_of_gene %>%
                    dplyr::filter(is.na(cis_activated_gene) | (cis_activated_gene == FALSE)) %>%
                    dplyr::filter(label == "SV")

                CNAs_others_of_gene <- CNAs_of_gene %>% 
                    dplyr::filter((is.na(cis_activated_gene) | (cis_activated_gene == FALSE)), cna_chrom == gene_info[["chrom"]])

                if((nrow(SVs_others_of_gene) > 0)|(nrow(CNAs_others_of_gene) > 0)){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale
                    CNAs_and_SVs_others_y <- c(CNAs_others_of_gene$sample_ID, SVs_others_of_gene$sample_ID)

                    CNAs_others_y_numeric <- as.numeric(as.factor(c(CNAs_and_SVs_others_y)))[0:length(CNAs_others_of_gene$sample_ID)]
                    SVs_others_y_numeric <- {if(nrow(SVs_others_of_gene) > 0) as.numeric(as.factor(c(CNAs_and_SVs_others_y)))[(length(CNAs_others_of_gene$sample_ID)+1):length(CNAs_and_SVs_others_y)] else numeric(0)}
                    p_SVs_and_CNAs_others_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(CNAs_and_SVs_others_y)))))
                    p_SVs_and_CNAs_others_y_labels <- levels(as.factor(CNAs_and_SVs_others_y))

                    p_SVs_and_CNAs_others <- SVs_others_of_gene %>%
                        ggplot2::ggplot() +
                        ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                        ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                        { if(length(SVs_others_y_numeric) > 0) ggplot2::geom_point(ggplot2::aes(x=pos, y = SVs_others_y_numeric, color=as.factor(label)), show.legend = F) } +
                        ggplot2::geom_errorbar(CNAs_others_of_gene, mapping = ggplot2::aes(xmin = cna_start, xmax = cna_end, y = CNAs_others_y_numeric), color = "green", width=0.5) + 
                        scale_color_mutation_type +
                        ggplot2::scale_y_continuous(breaks = p_SVs_and_CNAs_others_y_breaks, labels = p_SVs_and_CNAs_others_y_labels) +
                        ggplot2::labs(y="")+
                        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
                }else{
                    p_SVs_and_CNAs_others <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                # p_SVs_and_CNAs_others <- SVs_and_somatic_SNVs_of_gene %>%
                #     dplyr::filter(is.na(cis_activated_gene) | (cis_activated_gene == FALSE)) %>%
                #     dplyr::filter(label == "SV") %>%
                #     ggplot2::ggplot() +
                #     ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                #     ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                #     ggplot2::geom_point(ggplot2::aes(x=pos, y = as.factor(sample_ID), color=as.factor(label)), show.legend = F) +
                #     scale_color_mutation_type +
                #     ggplot2::geom_errorbar(dplyr::filter(CNAs_of_gene, is.na(cis_activated_gene) | (cis_activated_gene == FALSE), cna_chrom == gene_info[["chrom"]]), mapping = ggplot2::aes(xmin = cna_start, xmax = cna_end, y = as.factor(sample_ID)), color = "green") +
                #     ggplot2::labs(y="")

                SNVs_others_of_gene <- SVs_and_somatic_SNVs_of_gene %>%
                    dplyr::filter(is.na(cis_activated_gene) | (cis_activated_gene == FALSE)) %>%
                    dplyr::filter(label == "somatic_SNV")

                if(nrow(SNVs_others_of_gene) > 0){
                    # factor to numeric conversion is needed as ggbio::tracks throws error with factors in y scale

                    p_SNVs_others_y_numeric <- as.numeric(as.factor(SNVs_others_of_gene$sample_ID))
                    p_SNVs_others_y_breaks <- as.numeric(levels(as.factor(as.numeric(as.factor(SNVs_others_of_gene$sample_ID)))))
                    p_SNVs_others_y_labels <- levels(as.factor(SNVs_others_of_gene$sample_ID))

                    p_SNVs_others <- SNVs_others_of_gene %>%
                        ggplot2::ggplot() +
                        ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                        ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                        ggplot2::geom_point(ggplot2::aes(x=pos, y = p_SNVs_others_y_numeric, color=as.factor(label)), show.legend = F) +
                        scale_color_mutation_type +
                        ggplot2::scale_y_continuous(breaks = p_SNVs_others_y_breaks, labels = p_SNVs_others_y_labels) +
                        ggplot2::labs(y="")+
                        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

                }else{
                    p_SNVs_others <- ggplot2::ggplot() + ggplot2::geom_blank()
                }

                # p_SNVs_others <- SVs_and_somatic_SNVs_of_gene %>%
                #     dplyr::filter(is.na(cis_activated_gene) | (cis_activated_gene == FALSE)) %>%
                #     dplyr::filter(label == "somatic_SNV") %>%
                #     ggplot2::ggplot() +
                #     ggplot2::geom_vline(xintercept = gene_info[["start"]]) +
                #     ggplot2::geom_vline(xintercept = gene_info[["end"]]) + 
                #     ggplot2::geom_point(ggplot2::aes(x=pos, y = as.factor(sample_ID), color=as.factor(label)), show.legend = F) +
                #     scale_color_mutation_type +
                #     ggplot2::labs(y="")
                
                if(nrow(cis_activated_genes_by_gene_within_TAD) > 0){
                    genes_within_TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
                        cis_activated_genes_by_gene_within_TAD,
                        keep.extra.columns=TRUE
                    )
                    genes_GRangesList<- GenomicRanges::split(genes_within_TAD_GRanges, as.factor(genes_within_TAD_GRanges$gene_name))
                    
                    # apparently ggbio does not use type anymore -> let's see if this won't break
                    # old:
                    # p_genes <- ggplot2::ggplot() + ggbio::geom_alignment(genes_GRangesList, type = "type")
                    # new:
                    p_genes <- ggplot2::ggplot() + ggbio::geom_alignment(genes_GRangesList)

                    # browser();
                    gene_locus_plot <- ggbio::tracks(
                        "Chromosome" = p_ideo,
                        "ChIP-seq" = p_chipseq,
                        "Genehancer" = p_genehancer,
                        "SVs / CNAs\n(cis-act.)" = p_SVs_and_CNAs_cis_activated,
                        "SNVs\n(cis-act.)" = p_SNVs_cis_activated,
                        "SVs/CNAs\n(others)" = p_SVs_and_CNAs_others,
                        "SNVs\n(others)" = p_SNVs_others,
                        "Genes" = p_genes,
                        xlim = c(IRanges::start(TAD_range),IRanges::end(TAD_range)),
                        label.text.cex = 0.7
                    )
                }else{
                    gene_locus_plot <- ggbio::tracks(
                        "Chromosome" = p_ideo,
                        "ChIP-seq" = p_chipseq,
                        "Genehancer" = p_genehancer,
                        "SVs / CNAs\n(cis-act.)" = p_SVs_and_CNAs_cis_activated,
                        "SNVs\n(cis-act.)" = p_SNVs_cis_activated,
                        "SVs/CNAs\n(others)" = p_SVs_and_CNAs_others,
                        "SNVs\n(others)" = p_SNVs_others,
                        xlim = c(IRanges::start(TAD_range),IRanges::end(TAD_range)),
                        label.text.cex = 0.7
                    )
                }
            },
            c(
                "Registered S3 method overwritten by 'GGally'",
                "Constructing graphics"
            ))
            return(gene_locus_plot)
        }
        gene_locus_plot_path <- file.path(HTML_gene_figure_directory, paste0("gene_locus_plot___",GENE_NAME,".svg"))
        gene_locus_plot <- create_gene_locus_plot(
            min_TAD_start,
            max_TAD_end,
            chipseq_by_genes_within_TAD,
            genehancer_with_cis_activation_summary_of_gene_within_TAD,
            SVs_and_somatic_SNVs_of_gene,
            CNAs_of_gene,
            color_palette$scale_color_mutation_type,
            cis_activated_genes_by_gene_within_TAD,
            genome_name
        )
        # cat("Saving ggbio plot\n")
        svglite::svglite(gene_locus_plot_path)
        print(gene_locus_plot)
        dev.off()
    
    # create genehancer data tables -----------------------------------
    GENE_NAME_ESCAPED_DOTS <- stringr::str_replace_all(GENE_NAME, pattern="\\.", replacement="\\\\.") # to make bootstrap collapse work for gene_names with dots
    # cat("Creating Genehancer data tables\n")
    genehancer_locus_plot_ids <- character(length(genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$genehancer_id));
    genehancer_data_of_gene_HTML_list <- seq_len(length(genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$genehancer_id)) %>% 
        purrr::map(function(i) {
            gh_id <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$genehancer_id[i]
            gh_chrom <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$chrom_genehancer[i]
            gh_start <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$start_genehancer[i]
            gh_end <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$end_genehancer[i]
            gh_feature_name <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$feature_name[i]
            gh_min_TAD_start <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$variant_overlap_window_start_genehancer[i]
            gh_max_TAD_end <- genehancer_with_cis_activation_summary_of_gene_only_one_sample_per_genehancer$variant_overlap_window_end_genehancer[i]

            is_within_gene_TAD <- (
                (gh_start <= max_TAD_end) &
                (gh_end >= min_TAD_start) &
                (gh_chrom == gene_info[["chrom"]])
            )
            
            # CHANGE TO ALL - with col cis_activated
            SVs_of_this_genehancer <- SV_genehancer_of_gene %>% 
                dplyr::filter(genehancer_id == gh_id)

            CNAs_of_this_genehancer <- CNA_genehancer_of_gene %>% 
                dplyr::filter(genehancer_id == gh_id)

            somatic_SNVs_of_this_genehancer <- somatic_SNV_genehancer_of_gene %>% 
                dplyr::filter(genehancer_id == gh_id)

            genehancer_data_of_gene_SVs <- SVs_of_this_genehancer %>%
                dplyr::group_by(sample_ID, position_category_SV, cis_activated_gene) %>%
                dplyr::summarize(n = dplyr::n())

            genehancer_data_of_gene_CNAs <- CNAs_of_this_genehancer %>%
                dplyr::group_by(sample_ID, position_category_CNA, cis_activated_gene) %>%
                dplyr::summarize(n = dplyr::n())

            genehancer_data_of_gene_SNVs <- somatic_SNVs_of_this_genehancer %>%
                dplyr::group_by(sample_ID, cis_activated_gene) %>% 
                dplyr::summarize(n_within_genehancer = dplyr::n())

            # create genehancer gene locus plots
            TAD_as_string <- paste0(gh_chrom,"__",gh_min_TAD_start,"__",gh_max_TAD_end)
            if(!is_within_gene_TAD){
                if(TAD_as_string %in%  genehancer_locus_plot_ids){
                    cat("        Skipped creating Genehancer ggbio plot - it already exists\n")
                    genehancer_locus_plot_HTML <- paste0('
                        <h5>Genehancer Locus</h5>
                        <div class="text-center p-3">
                            <img class ="w-50" src="figures/genehancer_locus_plot___',GENE_NAME,'__genehancer__',TAD_as_string,'.svg"/>
                        </div>')
                }else{
                    # cat("Creating genehancer ggbioplot for ")
                    # cat(TAD_as_string)
                    # cat("\n")
                    genehancer_locus_plot_ids[i] <<- TAD_as_string #### <<- assigns to variable in parent scope
                    # create required data subsets for locus plot (like above for gene locus)
                    somatic_SNVs_of_this_genehancer_selected_cols <- somatic_SNVs_of_this_genehancer %>%
                        dplyr::select(chr = snv_chrom, pos = snv_pos, sample_ID, cis_activated_gene) %>% #
                        dplyr::mutate(label = "somatic_SNV")
                    
                    SVs_of_this_genehancer_left_break <- SVs_of_this_genehancer %>%
                        dplyr::select(chr = sv_break1_chrom, pos = sv_break1_pos, sample_ID, cis_activated_gene) %>%
                        dplyr::mutate(label = "SV")

                    SVs_of_this_genehancer_right_break <- SVs_of_this_genehancer %>%
                        dplyr::select(chr = sv_break2_chrom, pos = sv_break2_pos, sample_ID, cis_activated_gene) %>%
                        dplyr::mutate(label = "SV")

                    SVs_of_this_genehancer_all_breaks <- rbind(SVs_of_this_genehancer_left_break, SVs_of_this_genehancer_right_break) %>%
                        dplyr::distinct()

                    if(length(somatic_SNVs_of_this_genehancer_selected_cols$chr)>0){
                        GenomeInfoDb::seqlevelsStyle(somatic_SNVs_of_this_genehancer_selected_cols$chr) <- "UCSC"
                    }
                    if(length(CNAs_of_this_genehancer$cna_chrom)>0){
                        GenomeInfoDb::seqlevelsStyle(CNAs_of_this_genehancer$cna_chrom) <- "UCSC"
                    }
                    if(length(SVs_of_this_genehancer_all_breaks$chr)>0){
                        GenomeInfoDb::seqlevelsStyle(SVs_of_this_genehancer_all_breaks$chr) <- "UCSC"
                    }
                    
                    SVs_and_somatic_SNVs_of_this_genehancer <- rbind(somatic_SNVs_of_this_genehancer_selected_cols, SVs_of_this_genehancer_all_breaks, fill = T) %>%
                        dplyr::filter(chr == gh_chrom) %>% # to filter out SV breaks from other chroms -> as SVs_and_somatic_SNVs_of_gene is used in plotting the locus
                        dplyr::mutate(pos = as.numeric(pos))

                    genehancer_with_cis_activation_summary_of_gene_within_genehancer_TAD <- genehancer_with_cis_activation_summary_of_gene %>%
                        dplyr::filter(
                                (start_genehancer <= !!gh_max_TAD_end) &
                                (end_genehancer >= !!gh_min_TAD_start) &
                                (chrom_genehancer == gh_chrom)
                            )

                    chipseq_by_genes_within_genehancer_TAD <- processed_data$chipseq_by_genes_table %>%
                        dplyr::filter(gene_name == GENE_NAME) %>%
                        dplyr::filter(
                                (start <= !!gh_max_TAD_end) &
                                (end >= !!gh_min_TAD_start) &
                                (chrom == gh_chrom)
                            )
                    
                    cis_activated_genes_by_gene_within_genehancer_TAD <- processed_data$cis_activated_genes_by_gene %>% 
                        dplyr::filter(
                            (start <= !!gh_max_TAD_end) &
                            (end >= !!gh_min_TAD_start) &
                            (chrom == gh_chrom)
                        ) %>%
                        dplyr::mutate(type = "cds")

                    genehancer_locus_plot_path <- file.path(HTML_gene_figure_directory, paste0("genehancer_locus_plot___",GENE_NAME,"__genehancer__",TAD_as_string,".svg"))
                    genehancer_locus_plot <- create_gene_locus_plot (
                        min_TAD_start = gh_min_TAD_start,
                        max_TAD_end = gh_max_TAD_end,
                        chipseq_by_genes_within_TAD = chipseq_by_genes_within_genehancer_TAD,# maybe empty??
                        genehancer_with_cis_activation_summary_of_gene_within_TAD = genehancer_with_cis_activation_summary_of_gene_within_genehancer_TAD,
                        SVs_and_somatic_SNVs_of_gene = SVs_and_somatic_SNVs_of_this_genehancer,
                        CNAs_of_gene = CNAs_of_this_genehancer,
                        scale_color_mutation_type = color_palette$scale_color_mutation_type,
                        cis_activated_genes_by_gene_within_TAD = cis_activated_genes_by_gene_within_genehancer_TAD,
                        genome_name = genome_name
                    )
                    # cat("Saving ggbio plot\n")
                    svglite::svglite(genehancer_locus_plot_path)
                    print(genehancer_locus_plot)
                    dev.off()

                    genehancer_locus_plot_HTML <- paste0('
                        <h5>Genehancer Locus</h5>
                        <div class="text-center p-3">
                            <img class ="w-50" src="figures/genehancer_locus_plot___',GENE_NAME,'__genehancer__',TAD_as_string,'.svg"/>
                        </div>')
                }
            }else{
                genehancer_locus_plot_HTML <- ''
            }

            genehancer_data_of_gene_SV_HTML <- ifelse(nrow(genehancer_data_of_gene_SVs) > 0,
                paste0("
                <h5>SVs</h5>", 
                knitr::kable(genehancer_data_of_gene_SVs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ),
                "")

            genehancer_data_of_gene_CNA_HTML <- ifelse(nrow(genehancer_data_of_gene_CNAs) > 0,
                paste0("
                <h5>CNAs</h5>", 
                knitr::kable(genehancer_data_of_gene_CNAs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ),
                "")

            genehancer_data_of_gene_somatic_SNVs_HTML <- ifelse(nrow(genehancer_data_of_gene_SNVs) > 0,
                paste0("
                <h5>somatic SNVs</h5>", 
                knitr::kable(genehancer_data_of_gene_SNVs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ),
            "")

            genehancer_data_of_gene_title_HTML <- paste0("<h4>", gh_id, " - (", gh_chrom, ":", gh_start, "-", gh_end, ")</h4>")
            genehancer_data_of_gene_subtitle_HTML <- paste0("<p>",gh_feature_name, " - ",ifelse(is_within_gene_TAD, "within gene TAD", "outside gene TAD"), "</p>")
            genehancer_data_of_gene_HTML <- paste0(
                genehancer_data_of_gene_title_HTML, 
                genehancer_data_of_gene_subtitle_HTML,
                genehancer_locus_plot_HTML,
                genehancer_data_of_gene_SV_HTML, 
                genehancer_data_of_gene_CNA_HTML,
                genehancer_data_of_gene_somatic_SNVs_HTML)
        })

    complete_sample_data_of_genehancer_HTML <- paste0(
        unlist(genehancer_data_of_gene_HTML_list), 
        collapse = "\n")

    # create chipseq data tables -----------------------------------
    # input values to create
    #chipseq_with_cis_activation_summary_of_gene_only_one_sample_per_chipseq
    #SV_chipseq_of_gene
    #CNA_chipseq_of_gene
    #somatic_SNV_chipseq_of_gene
    #

    # cat("Creating Chipseq data tables\n")
    chipseq_data_of_gene_HTML_list <- seq_len(nrow(affected_chipseq_regions_of_gene)) %>% 
        purrr::map(function(i) {
            cs_id <- affected_chipseq_regions_of_gene$chipseq_id[i]
            cs_chrom <- affected_chipseq_regions_of_gene$chipseq_chrom[i]
            cs_start <- affected_chipseq_regions_of_gene$chipseq_start[i]
            cs_end <- affected_chipseq_regions_of_gene$chipseq_end[i]

            is_within_gene_TAD <- (
                (cs_start <= max_TAD_end) &
                (cs_end >= min_TAD_start) &
                (cs_chrom == gene_info[["chrom"]])
            )
            
            # CHANGE TO ALL - with col cis_activated
            chipseq_data_of_gene_SVs <- SV_chipseq_of_gene %>% 
                dplyr::filter(chipseq_id == cs_id) %>%
                dplyr::group_by(sample_ID, position_category_SV, cis_activated_gene) %>%
                dplyr::summarize(n = dplyr::n())

            chipseq_data_of_gene_CNAs <- CNA_chipseq_of_gene %>% 
                dplyr::filter(chipseq_id == cs_id) %>%
                dplyr::group_by(sample_ID, position_category_CNA, cis_activated_gene) %>%
                dplyr::summarize(n = dplyr::n())

            chipseq_data_of_gene_SNVs <- somatic_SNV_chipseq_of_gene %>% 
                dplyr::filter(chipseq_id == cs_id) %>%
                dplyr::group_by(sample_ID, cis_activated_gene) %>% 
                dplyr::summarize(n_within_chipseq = dplyr::n())

            chipseq_data_of_gene_SV_HTML <- ifelse(nrow(chipseq_data_of_gene_SVs) > 0,
                paste0("
                <h5>SVs</h5>", 
                knitr::kable(chipseq_data_of_gene_SVs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ),
                "")

            chipseq_data_of_gene_CNA_HTML <- ifelse(nrow(chipseq_data_of_gene_CNAs) > 0,
                paste0("
                <h5>CNAs</h5>", 
                knitr::kable(chipseq_data_of_gene_CNAs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ), "")

            chipseq_data_of_gene_somatic_SNVs_HTML <- ifelse(nrow(chipseq_data_of_gene_SNVs) > 0,
                paste0("
                <h5>somatic SNVs</h5>", 
                knitr::kable(chipseq_data_of_gene_SNVs, format = "html", table.attr = "class=\"table table-striped table-bordered\"")
                ), "")

            chipseq_data_of_gene_title_HTML <- paste0("<h4>", cs_id,"</h4>")
            chipseq_data_of_gene_subtitle_HTML <- paste0("<p>",ifelse(is_within_gene_TAD, "within gene TAD", "outside gene TAD"), "</p>")
            chipseq_data_of_gene_HTML <- paste0(
                chipseq_data_of_gene_title_HTML, 
                chipseq_data_of_gene_subtitle_HTML,
                chipseq_data_of_gene_SV_HTML, 
                chipseq_data_of_gene_CNA_HTML,
                chipseq_data_of_gene_somatic_SNVs_HTML)
        })

    complete_sample_data_of_chipseq_HTML <- paste0(
        unlist(chipseq_data_of_gene_HTML_list), 
        collapse = "\n")


    # create sample data tables -----------------------------------
    # cat("Creating sample data tables\n")
    sample_data_of_gene_HTML_list <- samples_of_gene_to_display_sample_data %>% purrr::map(
        function(s_id){
            n_CNAs <- cis_activation_summary_table_of_gene %>%
                dplyr::filter(sample_ID == s_id) %>%
                    dplyr::pull(n_CNAs)
            n_SVs <- cis_activation_summary_table_of_gene %>%
                dplyr::filter(sample_ID == s_id) %>%
                    dplyr::pull(n_SVs)
            
            is_sample_cis_activated <- cis_activation_summary_table_of_gene %>%
                dplyr::filter(sample_ID == s_id) %>%
                    dplyr::pull(cis_activated_gene) %>% { (. & (!is.na(.))) }

            is_sample_with_SVs_or_CNAs <- (n_CNAs > 0) | (n_SVs > 0)

            sample_data_of_gene_general <- cis_activation_summary_table_of_gene %>%
                dplyr::filter(sample_ID == s_id) %>%
                dplyr::select(
                    sample_id, subgroup, FPKM, n_markers, n_ASE_markers, n_cnv_markers, 
                    p_all, delta_abs_all, copynumber_tag_all, combined_p_value, 
                    mean_delta_abs, ABH, ASE_marker_run_overlap, 
                    passed_FPKM_filter, passed_ASE_filter, ASE_filter_rescued_oncogene, 
                    passed_mean_delta_abs_filter, OHE_used_reference, 
                    OHE_used_pvalue, passed_OHE_filter, cis_activated_gene)

            sample_data_of_gene_mutation_counts <- cis_activation_summary_table_of_gene %>% 
                dplyr::filter(sample_ID == s_id) %>%
                dplyr::select(
                sample_id, n_SVs_from_gene_TAD, n_CNAs_from_gene_TAD, 
                n_somatic_SNVs_from_gene_TAD, n_relevant_somatic_SNVs_from_gene_TAD_valid_for_cis_activated_genes, 
                n_SVs_from_genehancer_TAD, n_CNAs_from_genehancer_TAD, 
                n_somatic_SNVs_from_genehancer_TAD, n_relevant_somatic_SNVs_from_genehancer_TAD_valid_for_cis_activated_genes, 
                n_SVs, n_CNAs, n_somatic_SNVs, n_relevant_somatic_SNVs_valid_for_cis_activated_genes, 
                n_SVs_with_chipseq, n_CNAs_with_chipseq, n_somatic_SNVs_with_chipseq, 
                n_muts_from_gene_TAD, n_muts_from_genehancer_TAD, 
                n_muts)

            sample_data_of_gene_SVs <- processed_data$SV_all_merged_table %>%
                dplyr::filter(gene_name == GENE_NAME) %>%
                dplyr::filter(sample_ID == s_id) %>%
                dplyr::select(sv_break1_chrom, sv_break1_pos, sv_break2_chrom, sv_break2_pos, sv_type, eventInversion)

            sample_data_of_gene_CNAs <- processed_data$CNA_all_merged_table %>%
                dplyr::filter(gene_name == GENE_NAME) %>%
                dplyr::filter(sample_ID == s_id) %>%
                dplyr::select(cna_chrom, cna_start, cna_end, copy_number, CNA_type, log2)

            sample_data_of_gene_somatic_SNVs <- processed_data$somatic_SNV_all_merged_table %>%
                dplyr::filter(gene_name == GENE_NAME) %>%
                dplyr::filter(sample_ID == s_id) %>%
                dplyr::select(snv_chrom, snv_pos, ref, alt) # overlaps_chipseq

            # sample_data_of_gene_somatic_SNVs <- somatic_SNV_tf_binding_data_of_gene %>%
            #     dplyr::filter(sample_ID == s_id) %>%
            #     dplyr::group_by(snv_chrom, snv_pos, ref, alt) %>%
            #     dplyr::summarize(n_relevant_tf_binding_site = sum(relevant_tf_binding_site, na.rm = T))
            
            # iNSERT TF STUFF HERE

            sample_data_of_gene_general_HTML <- 
                    paste0('
                        <h5>General</h5>
                        <div class="overflow-scroll">',
                        knitr::kable(sample_data_of_gene_general, format = "html", table.attr = "class=\"table table-striped table-bordered\""),
                        '</div>'
                        )

            sample_data_of_gene_mut_counts_HTML <- 
                    paste0('
                        <h5>Mutation Counts</h5>
                        <table class="table table-striped table-bordered">
                            <thead>
                                <tr>
                                    <th>Variant Type</th>
                                    <th>n matched via gene TAD</th>
                                    <th>n matched via GeneHancer</th>
                                    <th>n affecting ChIPSeq</th>
                                    <th>total</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <th>SVs</th>
                                    <td>',sample_data_of_gene_mutation_counts$n_SVs_from_gene_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_SVs_from_genehancer_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_SVs_with_chipseq,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_SVs,'</td>
                                </tr>
                                <tr>
                                    <th>CNAs</th>
                                    <td>',sample_data_of_gene_mutation_counts$n_CNAs_from_gene_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_CNAs_from_genehancer_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_CNAs_with_chipseq,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_CNAs,'</td>
                                </tr>
                                <tr>
                                    <th>somatic SNVs/InDels</th>
                                    <td>',sample_data_of_gene_mutation_counts$n_somatic_SNVs_from_gene_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_somatic_SNVs_from_genehancer_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_somatic_SNVs_with_chipseq,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_somatic_SNVs,'</td>
                                </tr>
                                <tr>
                                    <th>relevant somatic SNVs/InDels</th>
                                    <td>',sample_data_of_gene_mutation_counts$n_relevant_somatic_SNVs_from_gene_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_relevant_somatic_SNVs_from_genehancer_TAD,'</td>
                                    <td>',
                                    # sample_data_of_gene_mutation_counts$n_relevant_somatic_SNVs_with_chipseq
                                    ' - '
                                    ,'</td>
                                    <td>',
                                    # sample_data_of_gene_mutation_counts$n_relevant_somatic_SNVs
                                    ' - '
                                    ,'</td>
                                </tr>
                                <tr>
                                    <th>all Mutations</th>
                                    <td>',sample_data_of_gene_mutation_counts$n_muts_from_gene_TAD,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_muts_from_genehancer_TAD,'</td>
                                    <td>',
                                    # sample_data_of_gene_mutation_counts$n_muts_with_chipseq
                                    ' - '
                                    ,'</td>
                                    <td>',sample_data_of_gene_mutation_counts$n_muts,'</td>
                                </tr>
                            </tbody>
                        </table>'
                    )

            sample_data_of_gene_SV_HTML <- ifelse(
                    nrow(sample_data_of_gene_SVs) > 0,
                    paste0('
                        <h5>SVs</h5>',
                        knitr::kable(sample_data_of_gene_SVs, format = "html", table.attr = "class=\"table table-striped table-bordered sv-table\"")
                    ),
                    ""
                )
            sample_data_of_gene_CNA_HTML <- ifelse(
                    nrow(sample_data_of_gene_CNAs) > 0,
                    paste0('
                        <h5>CNAs</h5>',
                        knitr::kable(sample_data_of_gene_CNAs, format = "html", table.attr = "class=\"table table-striped table-bordered cna-table\"")
                        ),
                    ""
                )
            sample_data_of_gene_somatic_SNVs_HTML <- ifelse(
                    nrow(sample_data_of_gene_somatic_SNVs) > 0,
                    paste0('<h5>Somatic SNVs</h5>',
                        knitr::kable(sample_data_of_gene_somatic_SNVs, format = "html", table.attr = "class=\"table table-striped table-bordered somatic-snv-table\"")
                        ),
                    ""
                )
            # INSERT TF binding stuff
            if(is.function(add_custom_HTML_to_gene_subdocument_sample_data)){
                cat(paste0("        running custom plugin: ", s_id, " - ", GENE_NAME,"\n"))
                sample_data_of_gene_custom_HTML <- add_custom_HTML_to_gene_subdocument_sample_data(s_id, GENE_NAME, gene_info[["chrom"]], sample_data_of_gene_SVs, sample_data_of_gene_CNAs, sample_data_of_gene_somatic_SNVs)
            }else{
                sample_data_of_gene_custom_HTML <- ""
            }
            

            sample_data_of_gene_HTML <- paste0('
                <h3>',
                s_id,
                ifelse(is_sample_cis_activated, '<span class="ms-3 badge bg-success fs-6 fw-light align-top">cis-activated</span>',''),
                ifelse(is_sample_with_SVs_or_CNAs, '<span class="ms-3 badge bg-primary fs-6 fw-light align-top">SVs/CNAs</span>',''),
                '
                <button
                    class="btn btn-secondary badge float-end ms-2 collapsed"
                    data-bs-toggle="collapse"
                    data-bs-target="#gene-card-', GENE_NAME_ESCAPED_DOTS, '-sample-',s_id,'-data"
                    aria-expanded="true"
                    aria-controls="gene-card-', GENE_NAME, '-sample-',s_id,'-data">show/hide</button></h3>
                <div class="border-bottom mb-3 collapse show" id="gene-card-', GENE_NAME, '-sample-',s_id,'-data">
                    <div class="igc-btns-container">
                        <h5>IGV</h5>
                        <div id="load-track-btn-container-',s_id,'" class="pb-3" style="display: none">
                            <p class="m-0">&#10003; Bam Files loaded</p>
                            <button class="btn btn-sm btn-primary load-track-btn" data-sample-id="',s_id,'">
                                Load Sample in IGV
                            </button>
                        </div>
                        <button class="btn btn-sm btn-secondary" onclick="document.getElementById(\'file-input-bam-',s_id,'\').click();">
                            Open Bam File
                        </button>
                        <button class="btn btn-sm btn-secondary" onclick="document.getElementById(\'file-input-index-',s_id,'\').click();">
                            Open Index File
                        </button>
                        <input
                            id="file-input-bam-',s_id,'"
                            class="file-input-bam"
                            type="file"
                            data-sample-id="',s_id,'"
                            data-load-track-btn-container-id="load-track-btn-container-',s_id,'"
                            style="display: none"
                        />
                        <input
                            id="file-input-index-',s_id,'"
                            class="file-input-index"
                            type="file"
                            data-sample-id="',s_id,'"
                            data-load-track-btn-container-id="load-track-btn-container-',s_id,'"
                            style="display: none"
                        />
                    </div>',
                    sample_data_of_gene_general_HTML,
                    sample_data_of_gene_mut_counts_HTML,
                    sample_data_of_gene_SV_HTML,
                    sample_data_of_gene_CNA_HTML,
                    sample_data_of_gene_somatic_SNVs_HTML,
                    sample_data_of_gene_custom_HTML,
                '</div>'
            )
            return(sample_data_of_gene_HTML)
        }
    )
    
    complete_sample_data_of_gene_HTML <- paste0(unlist(sample_data_of_gene_HTML_list), collapse = "\n")

    # html output -----------------------------------
    subdocument_HTML <- paste0('
        <html> 
            <head>
                <title>Revana Report - ', GENE_NAME ,'</title>
                <link
                    href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/css/bootstrap.min.css"
                    rel="stylesheet"
                    integrity="sha384-+0n0xVW2eSR5OomGNYDnhzAbDsOXxcvSN1TPprVMTNDbiYZCxYbOOl7+AMvyTG2x"
                    crossorigin="anonymous"
                />
                <script
                    src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/js/bootstrap.bundle.min.js"
                    integrity="sha384-gtEjrD/SeCtmISkJkNUaaKMoLD0//ElJ19smozuHV6z3Iehds+3Ulb9Bn9Plx0x4"
                    crossorigin="anonymous"
                ></script>
                <style>
                    .sv-pos-table,
                    .cna-pos-table {
                        table-layout: fixed;
                    }

                    .sv-pos-table td:first-child,
                    .sv-pos-table th:first-child  {
                        width: 15em;
                    }
                </style>
            </head>
            <body>
                <div class="container">
                    <div class = "card my-4 ', ifelse(gene_info[["gene_type"]] == "protein_coding", "analysis_by_gene_proteine_coding", "analysis_by_gene_not_proteine_coding"),'">
                        <h5 class="card-header">', GENE_NAME, '
                            <button class="btn btn-secondary badge float-end ms-2" data-bs-toggle="collapse" data-bs-target="#card-body-gene-', GENE_NAME_ESCAPED_DOTS, '" aria-expanded="true" aria-controls="card-body-gene-', GENE_NAME, '">show/hide</button>
                            <span class="badge bg-danger float-end">By gene</span>
                        </h5>
                        <div class="collapse show" id="card-body-gene-', GENE_NAME, '" style="">
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="p-4 w-50 border-end">
                                    <p>Chromosome: ',gene_info[["chrom"]],'</p>
                                    <p>Position: ',format_number(gene_info[["start"]]),' - ',format_number(gene_info[["end"]]),'</p>
                                    <p>Width: ',format_number(gene_info[["width"]]),'</p>
                                    <p>Strand: ',gene_info[["strand"]],'</p>
                                    <p>Gene Type: ',gene_info[["gene_type"]],'</p>
                                    <p>Imprinting Status: ',gene_info[["imprinting_status"]],'</p>
                                    <p>Cancer Gene: ',ifelse(gene_info[["is_cancer_gene"]],"yes","no"),'</p>
                                    <p>Role in Cancer: ',ifelse(is.na(gene_info[["cancer_gene_role_in_cancer"]]),"-",gene_info[["cancer_gene_role_in_cancer"]]),'</p>
                                    <p>TAD region: ', format_number(min_TAD_start),' - ', format_number(max_TAD_end),'</p>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/mutation_type_matrix_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/passed_filter_all_samples_upsetr___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/mutation_type_all_samples_upsetr___",GENE_NAME,".svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/expression_by_sample_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/expression_by_group_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/cov_ratio_by_sample_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/expression_copy_number_normalized_by_group_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    ',ifelse(n_cis_activated_samples_of_gene > 0,paste0(
                                    '<img class ="w-100 p-2" src="',paste0("figures/ASE_detection_type_matrix_plot___",GENE_NAME,".svg"),'"/>'),
                                    '<p>No cis-activated samples</p>'),'
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/allelic_imbalance_by_group_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/allelic_imbalance_expression_dot_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/allelic_imbalance_expression_dot_plot___",GENE_NAME,"_SV_CNA_highlighted.svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/allelic_imbalance_expression_copy_number_normalized_dot_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50">
                                    <img class ="w-100 p-2" src="',paste0("figures/allelic_imbalance_expression_copy_number_normalized_dot_plot___",GENE_NAME,"_SV_CNA_highlighted.svg"),'"/>
                                </div>
                            </div>
                            <div class="d-flex flex-wrap border-bottom">
                                <div class="w-50 border-end">
                                    <img class ="w-100 p-2" src="',paste0("figures/gene_locus_plot___",GENE_NAME,".svg"),'"/>
                                </div>
                                <div class="w-50 p-2">
                                    <div id="igv-window"></div>
                                    <div id="igv-btn-container" class="w-100 border border-primary h-100 d-flex align-items-center justify-content-center">
                                        <button class="btn btn-lg btn-primary" id="igv-start-btn" type="button">START IGV</button>
                                    </div>
                                </div>
                            </div>
                            <div class="p-3 border-bottom">
                                <h3>GeneHancer
                                    <button
                                        class="btn btn-secondary badge float-end ms-2 collapsed"
                                        data-bs-toggle="collapse"
                                        data-bs-target="#gene-card-', GENE_NAME_ESCAPED_DOTS, '-genehancer-data"
                                        aria-expanded="false"
                                        aria-controls="gene-card-', GENE_NAME, '-genehancer-data">show/hide</button>
                                </h3>
                                <div class="border-bottom mb-3 collapse" id="gene-card-', GENE_NAME, '-genehancer-data">',
                                complete_sample_data_of_genehancer_HTML,'
                                </div>
                            </div>
                            <div class="p-3 border-bottom">
                                <h3>ChIP-Seq
                                    <button
                                        class="btn btn-secondary badge float-end ms-2 collapsed"
                                        data-bs-toggle="collapse"
                                        data-bs-target="#gene-card-', GENE_NAME_ESCAPED_DOTS, '-chipseq-data"
                                        aria-expanded="false"
                                        aria-controls="gene-card-', GENE_NAME, '-chipseq-data">show/hide</button>
                                </h3>
                                <div class="border-bottom mb-3 collapse" id="gene-card-', GENE_NAME, '-chipseq-data">',
                                complete_sample_data_of_chipseq_HTML,'
                                </div>
                            </div>
                            <div class="p-3">
                            ',
                            complete_sample_data_of_gene_HTML,'
                            </div>
                        </div>
                    </div>
                </div>
                <script src="https://cdn.jsdelivr.net/npm/igv@2.10.5/dist/igv.min.js"></script>
                <script>
                    const igvDiv = document.getElementById("igv-window");
                    class IgvController {
                        constructor(igvOutElement) {
                            this.igvOutElement = igvOutElement;
                            this.browser = null;
                            this.isIgvOpen = false;
                            this.sampleBamFiles = [];
                            this.sampleBamIndexFiles = [];
                            this.onOpenIGVCompleteCallback = () => {}
                        }

                        setLocus(locus) {
                            console.log(`Going to ${locus}`)
                            if(this.isIgvOpen){
                                this.browser.search(locus)
                            }else{
                                this.locus = locus
                            }
                        }

                        setLocusBySpotAndPadding(chrom, location, padding) {
                            const start = location - padding;
                            const end = location + padding;
                            const locus = `${chrom}:${start}-${end}`
                            this.setLocus(locus)
                        }

                        setBrowser(browser){
                            this.browser = browser;
                        }

                        setOnOpenIGVCompleteCallback(callback){
                            this.onOpenIGVCompleteCallback = callback;
                        }

                        getOptions() {
                            const options = {
                                genome: "',genome_name,'",
                                locus: this.locus,
                            };
                            return options;
                        }

                        openIGV() {
                            const igvOut = this.igvOutElement;
                            const options = this.getOptions();
                            const onBrowserCreated = (browser) => {
                                console.log("Created IGV browser");
                                this.setBrowser(browser);
                                this.isIgvOpen = true;
                                this.onOpenIGVCompleteCallback();
                            };
                            igv.createBrowser(this.igvOutElement, options).then(onBrowserCreated);
                        }

                        loadTrack({label, file, indexFile, callback = () => {}}){
                            this.browser.loadTrack({
                                label: label,
                                url: file,
                                indexURL: indexFile,
                                format: "bam"
                            })
                            .then(callback)
                            .catch((e) => console.error(e))
                        }

                        addSampleBamFile(sampleId, file){
                            const sampleBamFiles = this.sampleBamFiles.filter((obj) => (!(obj.sampleId === sampleId)));
                            sampleBamFiles.push({
                                sampleId,
                                file
                            })
                            this.sampleBamFiles = sampleBamFiles
                        }

                        addSampleBamIndexFile(sampleId, file){
                            const sampleBamIndexFiles = this.sampleBamIndexFiles.filter((obj) => (!(obj.sampleId === sampleId)));
                            sampleBamIndexFiles.push({
                                sampleId,
                                file
                            })
                            this.sampleBamIndexFiles = sampleBamIndexFiles
                        }

                        getSampleBamFileBySampleId(sampleId) {
                            const sampleBamFileObj = this.sampleBamFiles.find(
                                (obj) => obj.sampleId === sampleId
                            );
                            return sampleBamFileObj.file;
                        }

                        getSampleBamIndexFileBySampleId(sampleId) {
                            const sampleBamIndexFileObj = this.sampleBamIndexFiles.find(
                                (obj) => obj.sampleId === sampleId
                            );
                            return sampleBamIndexFileObj.file;
                        }

                        areSampleFilesAvailable(sampleId){
                            const sampleBamFile = this.getSampleBamFileBySampleId(sampleId);
                            const sampleBamIndexFile = this.getSampleBamIndexFileBySampleId(sampleId);
                            return((!!sampleBamFile)&&(!!sampleBamIndexFile))
                        }

                        loadTrackBySampleId(sampleId){
                            if(this.areSampleFilesAvailable(sampleId) && this.isIgvOpen){
                                const bamFile = this.getSampleBamFileBySampleId(sampleId);
                                const bamIndexFile = this.getSampleBamIndexFileBySampleId(sampleId);

                                this.loadTrack({
                                    label: sampleId, 
                                    file: bamFile,
                                    indexFile: bamIndexFile,
                                    callback: () => {console.log("track loaded")}
                                })
                            }else{
                                if(!this.areSampleFilesAvailable(sampleId)){
                                    window.alert("sample files not available")
                                }
                                if(!this.isIgvOpen){
                                    window.alert("Please start IGV Plugin first!")
                                }
                            }
                        }
                    }

                    const igvController = new IgvController(igvDiv);
                    igvController.setLocus("',gene_info[["chrom"]],':',min_TAD_start,'-',max_TAD_end,'")



                    const igvBtn = document.getElementById("igv-start-btn");
                    const igvBtnContainer = document.getElementById("igv-btn-container");
                    const bamFileInputsElements = Array.from(
                        document.getElementsByClassName("file-input-bam")
                    );
                    const bamIndexFileInputsElements = Array.from(
                        document.getElementsByClassName("file-input-index")
                    );
                    const loadTrackBtnElements = Array.from(
                        document.getElementsByClassName("load-track-btn")
                    );
                    const svTableRowElements = Array.from(
                        document.querySelectorAll("table.sv-table > tbody > tr")
                    );
                    const cnaTableRowElements = Array.from(
                        document.querySelectorAll("table.cna-table > tbody > tr")
                    );
                    const somaticSnvTableRowElements = Array.from(
                        document.querySelectorAll("table.somatic-snv-table > tbody > tr")
                    );



                    const hideElement = (element) => {
                        element.style.display = "none";
                    }

                    const displayElement = (element) => {
                        element.style.display = "block";
                    }

                    const removeElement = (element) => {
                        element.remove();
                    }



                    const handleIgvBtnClick = () => {
                        igvController.openIGV();
                        removeElement(igvBtnContainer);
                    }

                    const handleBamFileInputChange = (e) => {
                        const file = e.target.files[0];
                        const sampleId = e.target.dataset.sampleId;
                        igvController.addSampleBamFile(sampleId, file)

                        if(igvController.areSampleFilesAvailable(sampleId)) {
                            const loadTrackBtnContainerElement = document.getElementById(e.target.dataset.loadTrackBtnContainerId)
                            displayElement(loadTrackBtnContainerElement)
                        }
                    }

                    const handleBamIndexFileInputChange = (e) => {
                        const file = e.target.files[0];
                        const sampleId = e.target.dataset.sampleId;
                        igvController.addSampleBamIndexFile(sampleId, file)

                        if(igvController.areSampleFilesAvailable(sampleId)) {
                            const loadTrackBtnContainerElement = document.getElementById(e.target.dataset.loadTrackBtnContainerId)
                            displayElement(loadTrackBtnContainerElement)
                        }
                    }

                    const handleLoadTrackBtnClick = (e) => {
                        const sampleId = e.target.dataset.sampleId;
                        igvController.loadTrackBySampleId(sampleId)
                    }



                    igvBtn.addEventListener("click", handleIgvBtnClick);

                    bamFileInputsElements.forEach((element) => {
                        element.addEventListener("change", (e) => handleBamFileInputChange(e));
                    });

                    bamIndexFileInputsElements.forEach((element) => {
                        element.addEventListener("change", (e) =>
                            handleBamIndexFileInputChange(e)
                        );
                    });

                    loadTrackBtnElements.forEach((element) => {
                        element.addEventListener("click", (e) =>
                            handleLoadTrackBtnClick(e)
                        );
                    });

                    const addViewInIgvButtons = () => {

                        svTableRowElements.forEach((element) => {
                            const tds = element.children;
                            const chrom1 = String(tds[0].innerHTML).trim();
                            const pos1 = parseInt(tds[1].innerHTML);
                            const chrom2 = String(tds[2].innerHTML).trim();
                            const pos2 = parseInt(tds[3].innerHTML);

                            const viewInIgvBtn1 = document.createElement("button");
                            viewInIgvBtn1.innerHTML = "View in IGV";
                            viewInIgvBtn1.className = "btn btn-sm btn-secondary py-0 px-1 my-0 mx-1"
                            viewInIgvBtn1.addEventListener("click", (e) => igvController.setLocusBySpotAndPadding(chrom1, pos1, 2000))

                            const viewInIgvBtn2 = document.createElement("button");
                            viewInIgvBtn2.innerHTML = "View in IGV";
                            viewInIgvBtn2.className = "btn btn-sm btn-secondary py-0 px-1 my-0 mx-1"
                            viewInIgvBtn2.addEventListener("click", (e) => igvController.setLocusBySpotAndPadding(chrom2, pos2, 2000))

                            tds[1].appendChild(viewInIgvBtn1)
                            tds[3].appendChild(viewInIgvBtn2)
                        })

                        cnaTableRowElements.forEach((element) => {
                            const tds = element.children;
                            const chrom = String(tds[0].innerHTML).trim();
                            const pos1 = parseInt(tds[1].innerHTML);
                            const pos2 = parseInt(tds[2].innerHTML);

                            const viewInIgvBtn1 = document.createElement("button");
                            viewInIgvBtn1.innerHTML = "View in IGV";
                            viewInIgvBtn1.className = "btn btn-sm btn-secondary py-0 px-1 my-0 mx-1"
                            viewInIgvBtn1.addEventListener("click", (e) => igvController.setLocusBySpotAndPadding(chrom, pos1, 2000))

                            const viewInIgvBtn2 = document.createElement("button");
                            viewInIgvBtn2.innerHTML = "View in IGV";
                            viewInIgvBtn2.className = "btn btn-sm btn-secondary py-0 px-1 my-0 mx-1"
                            viewInIgvBtn2.addEventListener("click", (e) => igvController.setLocusBySpotAndPadding(chrom, pos2, 2000))

                            tds[1].appendChild(viewInIgvBtn1)
                            tds[2].appendChild(viewInIgvBtn2)
                        })

                        somaticSnvTableRowElements.forEach((element) => {
                            const tds = element.children;
                            const chrom = String(tds[0].innerHTML).trim();
                            const pos = parseInt(tds[1].innerHTML);

                            const viewInIgvBtn = document.createElement("button");
                            viewInIgvBtn.innerHTML = "View in IGV";
                            viewInIgvBtn.className = "btn btn-sm btn-secondary py-0 px-1 my-0 mx-1"
                            viewInIgvBtn.addEventListener("click", (e) => igvController.setLocusBySpotAndPadding(chrom, pos, 2000))

                            tds[1].appendChild(viewInIgvBtn)
                        })
                    }

                    igvController.setOnOpenIGVCompleteCallback(addViewInIgvButtons)
                </script>
            </body>
        </html>
    ')

    overview_HTML <- paste0('
        <a href="genes/',file.path(GENE_NAME,"gene_overview.html"),'" target="_blank" class="list-group-item list-group-item-action ',ifelse(gene_info[["gene_type"]] == "protein_coding","analysis_by_gene_proteine_coding","analysis_by_gene_not_proteine_coding"),'" aria-current="true">
            <div class="d-flex w-100 justify-content-between">
            <h5 class="mb-1">',GENE_NAME,'</h5>
            <small>',gene_info[["gene_type"]],'</small>
            </div>
            <p class="mb-1">',gene_info[["chrom"]],' : ',format_number(gene_info[["start"]]),' - ',format_number(gene_info[["end"]]),'</p>
            <small>
            ',n_cis_activated_samples_of_gene,' cis-activated samples</br>
            ',n_cis_activated_samples_of_gene_with_any_SVs,' cis-activated samples with SVs</br>
            ',n_cis_activated_samples_of_gene_with_any_CNAs,' cis-activated samples with CNAs</br>
            ',n_cis_activated_samples_of_gene_with_any_somatic_SNVs,' cis-activated samples with somatic SNVs
            </small>
        </a>
    ')

    return(
        list(
            gene_name = GENE_NAME,
            subdocument_HTML = subdocument_HTML,
            subdocument_path = subdocument_path,
            overview_HTML = overview_HTML
        )
    )
}