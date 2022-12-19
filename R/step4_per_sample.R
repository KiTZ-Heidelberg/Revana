run_step4_per_sample_chunkwise_SNV_links <- function(output_dir,
                                                     sample_id,
                                                     TAD_file_path,
                                                     genehancer_Rds_path,
                                                     chipseq_file_path,
                                                     chipseq_by_genes_file_path,
                                                     check_paths = TRUE,
                                                     store_somatic_SNV_tf_binding_data_genes_all = FALSE,
                                                     store_somatic_SNV_tf_binding_data_genehancer_all = FALSE,
                                                     run_tf_binding_site_analysis = TRUE,
                                                     verbose = FALSE) {
    # get all required paths ---------------------------------
    cohort_expression_reference_path <- file.path(output_dir, "cohort.expression.reference.Rds")
    sample_dir_path <- file.path(output_dir, sample_id)
    marker_summary_file_path <- file.path(sample_dir_path, paste0(sample_id, ".geneanno_marker_summary.Rds"))
    SV_file_path <- file.path(sample_dir_path, paste0(sample_id, ".SV.Rds"))
    CNA_file_path <- file.path(sample_dir_path, paste0(sample_id, ".CNA.Rds"))
    somatic_SNV_file_path <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.Rds"))
    somatic_SNV_tf_binding_data_file_path <- file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV_tf_binding_data.Rds"))

    # check all required files/directories -------------------
    if (check_paths == TRUE) {
        check_sample_dir_existence(sample_dir_path, sample_id)
        check_file_existence(TAD_file_path, name_of_file_type = "TAD file")
        check_file_existence_auto_created(cohort_expression_reference_path, name_of_file_type = "cohort.expression.reference.Rds file")
        check_file_existence_auto_created(marker_summary_file_path, name_of_file_type = "geneanno_marker_summary.Rds file")
        check_file_existence_auto_created(SV_file_path, name_of_file_type = "SV.Rds file")
        check_file_existence_auto_created(CNA_file_path, name_of_file_type = "CNA.Rds file")

        if (run_tf_binding_site_analysis) {
            check_file_existence_auto_created(somatic_SNV_tf_binding_data_file_path, name_of_file_type = "somatic_SNV_tf_binding_data.Rds file")
        }
    }




    # CALCULATE OHE ##########################################

    # import required data -----------------------------------
    cohort_expression_reference <- readRDS(cohort_expression_reference_path)


    # process data  ------------------------------------------------
    OHE_table <- calculate_OHE(cohort_expression_reference, sample_id)


    # store results ------------------------------------------------
    saveRDS(
        OHE_table,
        file = file.path(sample_dir_path, paste0(sample_id, ".OHE.Rds"))
    )

    data.table::fwrite(
        OHE_table,
        file = file.path(sample_dir_path, paste0(sample_id, ".OHE.txt")),
        sep = "\t"
    )
    rm(cohort_expression_reference)
    gc() # free memory

    log_msg("OHE calculated", verbose = verbose, sample_id = sample_id)



    # CALL CIS-ACTIVATION #####################################

    # import required data -----------------------------------
    marker_summary <- readRDS(marker_summary_file_path)

    # process data  ------------------------------------------------
    cis_activation_summary <- call_cis_activation(OHE_table, marker_summary)

    log_msg("cis-activation called", verbose = verbose, sample_id = sample_id)


    # ADD TADs TO GENEs #####################################

    # import required data -----------------------------------
    TADs <- import_TAD_file(TAD_file_path)

    # process data  ------------------------------------------------
    cis_activation_summary_TADs <- add_TAD_boundaries_to_cis_activation_summary(cis_activation_summary, TADs)
    add_variant_overlap_window(cis_activation_summary_TADs)

    # store results ------------------------------------------------
    saveRDS(
        cis_activation_summary_TADs, # changed cis_activation_summary -> cis_activation_summary_TADs
        file = file.path(sample_dir_path, paste0(sample_id, ".cis.activation.summary.Rds"))
    )

    data.table::fwrite(
        cis_activation_summary_TADs,
        file = file.path(sample_dir_path, paste0(sample_id, ".cis.activation.summary.txt")),
        sep = "\t"
    )

    saveRDS(
        cis_activation_summary_TADs[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".cis.activated.genes.Rds"))
    )

    data.table::fwrite(
        cis_activation_summary_TADs[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".cis.activated.genes.txt")),
        sep = "\t"
    )

    rm(OHE_table, marker_summary)
    gc() # free memory

    log_msg("linked TADs to genes", verbose = verbose, sample_id = sample_id)

    # PREPARE GENEHANCER DATA ################################

    # import required data -----------------------------------
    genehancer <- readRDS(genehancer_Rds_path)

    # process data -------------------------------------------
    genehancer_with_cis_activation_summary <- merge_genehancer_with_cis_activation_summary(genehancer, cis_activation_summary_TADs)
    genehancer_with_cis_activation_summary_TADs <- add_TAD_boundaries_to_genehancer_with_cis_activation_summary(genehancer_with_cis_activation_summary, TADs)
    add_variant_overlap_window_to_genehancer_with_cis_activation_summary(genehancer_with_cis_activation_summary_TADs)

    # store results ------------------------------------------------
    saveRDS(
        genehancer_with_cis_activation_summary_TADs,
        file = file.path(sample_dir_path, paste0(sample_id, ".genehancer_with_cis_activation_summary.Rds"))
    )

    data.table::fwrite(
        genehancer_with_cis_activation_summary_TADs,
        file = file.path(sample_dir_path, paste0(sample_id, ".genehancer_with_cis_activation_summary.txt")),
        sep = "\t"
    )

    log_msg("GeneHancer data prepared", verbose = verbose, sample_id = sample_id)


    # PREPARE CHIPSEQ DATA ################################

    # import required data -----------------------------------
    chipseq <- readRDS(chipseq_file_path)
    data.table::setDT(chipseq)
    data.table::setnames(chipseq, c("chrom", "start", "end"), c("chipseq_chrom", "chipseq_start", "chipseq_end"))
    chipseq_by_genes <- readRDS(chipseq_by_genes_file_path)


    # process data -------------------------------------------
    # chipseq
    chipseq_with_TADs <- add_TAD_boundaries_to_chipseq(chipseq, TADs)
    add_variant_overlap_window_to_chipseq(chipseq_with_TADs)


    # chipseq_by_gene
    chipseq_with_cis_activation_summary <- merge_chipseq_gene_combination_with_cis_activation_summary(chipseq_by_genes, cis_activation_summary_TADs)
    chipseq_with_cis_activation_summary_TADs <- add_TAD_boundaries_to_chipseq_with_cis_activation_summary(chipseq_with_cis_activation_summary, TADs)
    add_variant_overlap_window_to_chipseq_with_cis_activation_summary(chipseq_with_cis_activation_summary_TADs)



    # store results ------------------------------------------------
    saveRDS(
        chipseq_with_cis_activation_summary_TADs,
        file = file.path(sample_dir_path, paste0(sample_id, ".chipseq_with_cis_activation_summary.Rds"))
    )

    data.table::fwrite(
        chipseq_with_cis_activation_summary_TADs,
        file = file.path(sample_dir_path, paste0(sample_id, ".chipseq_with_cis_activation_summary.txt")),
        sep = "\t"
    )

    log_msg("ChIP-Seq data prepared", verbose = verbose, sample_id = sample_id)

    # LINK SNVs TO GENES AND GENEHANCERS ######################
    # import required data -----------------------------------
    somatic_SNVs <- readRDS(somatic_SNV_file_path)
    data.table::setDT(somatic_SNVs)
    data.table::setnames(somatic_SNVs, c("chrom", "pos"), c("snv_chrom", "snv_pos"))

    # process data  ------------------------------------------------
    add_chipseq_info_to_SNVs(somatic_SNVs, chipseq)
    somatic_SNV_gene_combinations <- link_SNVs_to_genes(somatic_SNVs, cis_activation_summary_TADs)

    somatic_SNV_genehancer_combinations <- link_SNVs_to_genehancer(somatic_SNVs, genehancer_with_cis_activation_summary_TADs)

    somatic_SNV_chipseq_combinations <- link_SNVs_to_chipseq(somatic_SNVs, chipseq_with_cis_activation_summary_TADs)

    # store results ------------------------------------------------
    saveRDS(
        somatic_SNV_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genes.all.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genes.all.txt")),
        sep = "\t"
    )

    saveRDS(
        somatic_SNV_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genes.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genes.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        somatic_SNV_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genehancer.all.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genehancer.all.txt")),
        sep = "\t"
    )

    saveRDS(
        somatic_SNV_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genehancer.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.genehancer.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        somatic_SNV_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.chipseq.all.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.chipseq.all.txt")),
        sep = "\t"
    )

    saveRDS(
        somatic_SNV_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.chipseq.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        somatic_SNV_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".somatic_SNV.chipseq.cis_activated_only.txt")),
        sep = "\t"
    )

    rm(somatic_SNVs, somatic_SNV_gene_combinations, somatic_SNV_genehancer_combinations, somatic_SNV_chipseq_combinations)
    gc() # free memory

    log_msg("somatic SNVs linked to genes (and GeneHancers)", verbose = verbose, sample_id = sample_id)


    # LINK SNVs TF BINDING DATA TO GENES AND GENEHANCERS (optional) ######################
    if (run_tf_binding_site_analysis) {
        tf_binding_site_step4_workflow(
            sample_dir_path,
            sample_id,
            somatic_SNV_tf_binding_data_file_path,
            chipseq,
            cis_activation_summary_TADs,
            genehancer_with_cis_activation_summary_TADs,
            store_somatic_SNV_tf_binding_data_genes_all,
            store_somatic_SNV_tf_binding_data_genehancer_all
        )

        log_msg("somatic SNVs TF binding data linked to genes (and GeneHancers)", verbose = verbose, sample_id = sample_id)
    }

    # TODO: add new output files to step 5



    # LINK CNAs TO GENES #####################################

    # import required data -----------------------------------
    CNAs <- readRDS(CNA_file_path)
    data.table::setDT(CNAs)
    data.table::setnames(CNAs, c("chrom", "start", "end"), c("cna_chrom", "cna_start", "cna_end"))

    # process data  ------------------------------------------------
    CNA_gene_combinations <- link_CNAs_to_genes(CNAs, cis_activation_summary_TADs)

    CNA_genehancer_combinations <- link_CNAs_to_genehancer(CNAs, genehancer_with_cis_activation_summary_TADs)
    CNA_genehancer_combinations <- evaluate_CNA_genehancer_variant_position(CNA_genehancer_combinations)

    CNA_chipseq_combinations <- link_CNA_gene_combinations_to_chipseq(CNA_gene_combinations, chipseq_with_TADs)
    CNA_chipseq_combinations <- evaluate_CNA_chipseq_variant_position(CNA_chipseq_combinations)

    # store results ------------------------------------------------
    saveRDS(
        CNA_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genes.all.Rds"))
    )

    data.table::fwrite(
        CNA_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genes.all.txt")),
        sep = "\t"
    )

    saveRDS(
        CNA_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genes.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        CNA_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genes.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        CNA_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genehancer.all.Rds"))
    )

    data.table::fwrite(
        CNA_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genehancer.all.txt")),
        sep = "\t"
    )

    saveRDS(
        CNA_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genehancer.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        CNA_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.genehancer.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        CNA_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.chipseq.all.Rds"))
    )

    data.table::fwrite(
        CNA_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.chipseq.all.txt")),
        sep = "\t"
    )

    saveRDS(
        CNA_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.chipseq.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        CNA_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".CNA.chipseq.cis_activated_only.txt")),
        sep = "\t"
    )

    rm(CNAs, CNA_gene_combinations, CNA_genehancer_combinations, CNA_chipseq_combinations)
    gc() # free memory

    log_msg("CNAs linked to genes (and GeneHancers, ChIP-Seq)", verbose = verbose, sample_id = sample_id)


    # LINK SVs TO GENES #####################################

    # import required data -----------------------------------
    SVs <- readRDS(SV_file_path)
    data.table::setDT(SVs)
    data.table::setnames(
        SVs,
        c("chrom1", "pos1", "chrom2", "pos2"),
        c("sv_break1_chrom", "sv_break1_pos", "sv_break2_chrom", "sv_break2_pos")
    )

    # process data  ------------------------------------------------
    SVs_with_TAD_combinations <- add_TAD_combination_to_SVs(SVs, TADs)
    SV_gene_combinations <- link_SVs_to_genes(SVs_with_TAD_combinations, cis_activation_summary_TADs)

    SV_genehancer_combinations <- link_SVs_to_genehancer(SVs_with_TAD_combinations, genehancer_with_cis_activation_summary_TADs)
    SV_genehancer_combinations <- evaluate_SV_genehancer_variant_position(SV_genehancer_combinations)

    SV_chipseq_combinations <- link_SV_gene_combinations_to_chipseq(SV_gene_combinations, chipseq_with_TADs)
    SV_chipseq_combinations <- evaluate_SV_chipseq_variant_position(SV_chipseq_combinations)

    # store results ------------------------------------------------
    saveRDS(
        SV_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genes.all.Rds"))
    )

    data.table::fwrite(
        SV_gene_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genes.all.txt")),
        sep = "\t"
    )

    saveRDS(
        SV_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genes.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        SV_gene_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genes.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        SV_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genehancer.all.Rds"))
    )

    data.table::fwrite(
        SV_genehancer_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genehancer.all.txt")),
        sep = "\t"
    )

    saveRDS(
        SV_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genehancer.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        SV_genehancer_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.genehancer.cis_activated_only.txt")),
        sep = "\t"
    )

    saveRDS(
        SV_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.chipseq.all.Rds"))
    )

    data.table::fwrite(
        SV_chipseq_combinations,
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.chipseq.all.txt")),
        sep = "\t"
    )

    saveRDS(
        SV_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.chipseq.cis_activated_only.Rds"))
    )

    data.table::fwrite(
        SV_chipseq_combinations[cis_activated_gene == TRUE],
        file = file.path(sample_dir_path, paste0(sample_id, ".SV.chipseq.cis_activated_only.txt")),
        sep = "\t"
    )

    rm(SVs, SVs_with_TAD_combinations, SV_gene_combinations)
    gc() # free memory

    log_msg("SVs linked to genes (and GeneHancers, ChIP-Seq)", verbose = verbose, sample_id = sample_id)


    rm(cis_activation_summary_TADs)
    gc() # free memory
}


#' Run Step 4 of the Revana workflow on a cohort/subgroup of tumor sample data - Parallelized implementation
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param TAD_file_path Path to the TAD file. The TAD file contains genomic coordinates of Topologically Associated Domains (TADs). See the documentation for the exact format of this file.
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#' @param use_parallelization should parallelization be used to speed up Revana (TRUE / FALSE), defaults to TRUE
#'
#' @details
#'
#' This implementation of step 4 runs the samples in parallel. This is the optimized and default implementation.
#'
#' For each included sample step 4 of the revana workflow ...
#' * calls cis-activation
#' * adds the TADs to the genes
#' * associates genes with regulatory regions from GeneHancer and ChIP-Seq data
#' * matches somatic SNVs and InDels with genes and genehancers
#' * matches somatic SNVs and InDels TF binding data with genes and genehancers
#' * matches copy number aberrations (CNAs) with genes and genehancers
#' * matches structural variants (SVs) with genes and genehancers
#' @export
run_step4_per_cohort <- function(paths_file_path,
                                 output_dir,
                                 TAD_file_path,
                                 run_tf_binding_site_analysis = TRUE,
                                 verbose = FALSE,
                                 use_parallelization = TRUE) {
    # import paths file  -------------------------------------
    paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

    # get all required paths ---------------------------------
    sample_ids <- paths$sample_id
    sample_dir_paths <- file.path(output_dir, sample_ids)

    cohort_expression_reference_path <- file.path(output_dir, "cohort.expression.reference.Rds")
    genehancer_Rds_path <- file.path(output_dir, "genehancer_ref.Rds")
    chipseq_file_path <- file.path(output_dir, "chipseq.Rds")
    chipseq_by_genes_file_path <- file.path(output_dir, "chipseq_by_gene.Rds")

    marker_summary_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".geneanno_marker_summary.Rds"))
    SV_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".SV.Rds"))
    CNA_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".CNA.Rds"))
    somatic_SNV_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".somatic_SNV.Rds"))
    somatic_SNV_tf_binding_data_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".somatic_SNV_tf_binding_data.Rds"))

    # check all required files/directories -------------------
    purrr::walk2(sample_dir_paths, sample_ids, function(s_dir, s_id) {
        check_sample_dir_existence(s_dir, s_id)
    })

    check_file_existence(TAD_file_path, name_of_file_type = "TAD file")

    check_file_existence_auto_created(cohort_expression_reference_path, name_of_file_type = "cohort.expression.reference.Rds file")
    check_file_existence_auto_created(genehancer_Rds_path, name_of_file_type = "genehancer_ref.Rds file")
    check_file_existence_auto_created(chipseq_file_path, name_of_file_type = "chipseq.Rds file")
    check_file_existence_auto_created(chipseq_by_genes_file_path, name_of_file_type = "chipseq_by_gene.Rds file")

    purrr::walk(marker_summary_file_paths, function(path) {
        check_file_existence_auto_created(path, name_of_file_type = "ASE marker summary .Rds file")
    })
    purrr::walk(SV_file_paths, function(path) {
        check_file_existence_auto_created(path, name_of_file_type = "SV.Rds file")
    })
    purrr::walk(CNA_file_paths, function(path) {
        check_file_existence_auto_created(path, name_of_file_type = "CNA.Rds file")
    })
    purrr::walk(somatic_SNV_file_paths, function(path) {
        check_file_existence_auto_created(path, name_of_file_type = "somatic_SNV.Rds file")
    })
    if (run_tf_binding_site_analysis) {
        purrr::walk(somatic_SNV_tf_binding_data_file_paths, function(path) {
            check_file_existence_auto_created(path, name_of_file_type = "somatic_SNV_tf_binding_data.Rds file")
        })
    }

    # run step 4 on all samples ------------------------------
    n_samples <- nrow(paths)

    if(use_parallelization == TRUE){
        if(.Platform$OS.type == "unix") {
            future::plan(future::multicore)
        }else{
            future::plan(future::multisession)
        }

        cat("Step 4: running in parallel on ")
        cat(future::nbrOfWorkers())
        cat(" workers\n")


        progressr::with_progress({
            p <- progressr::progressor(steps = n_samples)
            unique(sample_ids) %>% furrr::future_walk(function(s_id) {
                ## for hypermutators
                # uncomment for easier debugging:
                # withCallingHandlers(message=handle_message, warning = handle_warning, {
                run_step4_per_sample_chunkwise_SNV_links(
                    output_dir,
                    s_id,
                    TAD_file_path,
                    genehancer_Rds_path,
                    chipseq_file_path,
                    chipseq_by_genes_file_path,
                    run_tf_binding_site_analysis = run_tf_binding_site_analysis,
                    verbose = verbose
                )
                # })

                # DEPRECATED:
                # normal
                # run_step4_per_sample(output_dir, s_id, TAD_file_path)

                p() # progress ++
            })
        })
    }else{
        cat("Step 4: running serially\n")


        progressr::with_progress({
            p <- progressr::progressor(steps = n_samples)
            for (s_id in unique(sample_ids)){
                ## for hypermutators
                run_step4_per_sample_chunkwise_SNV_links(
                    output_dir,
                    s_id,
                    TAD_file_path,
                    genehancer_Rds_path,
                    chipseq_file_path,
                    chipseq_by_genes_file_path,
                    run_tf_binding_site_analysis = run_tf_binding_site_analysis,
                    verbose = verbose
                )

                # DEPRECATED:
                # normal
                # run_step4_per_sample(output_dir, s_id, TAD_file_path)

                p() # progress ++
            }
        })
    }
}