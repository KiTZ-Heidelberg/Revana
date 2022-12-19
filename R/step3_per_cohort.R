#' Run Step 3 of the Revana workflow on a cohort/subgroup of tumor sample data
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param TAD_file_path Path to the TAD file. The TAD file contains genomic coordinates of Topologically Associated Domains (TADs). See the documentation for the exact format of this file.
#' @param genehancer_ref_file_path Path to the reference for GeneHancer. See the documentation for the exact format of this file and how to obtain it.
#' @param chipseq_file_path Path to the ChIP-Seq file. The ChIP-Seq file contains genomic coordinates of regulatory active regions and can be obtained from experimental data or online less specific online resources. This file is optional, although an empty dummy file has to be used, if no ChIP-Seq data is supplied. See the documentation for the exact format of this file.
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#' @details
#' Step 3 of the revana workflow ...
#' * creates the internal OHE reference
#' * processes GeneHancer reference data
#' * processes ChIPseq input data
#' * links the ChIPseq regions to their respective TADs
#'
#' @export
run_step3_per_cohort <- function(paths_file_path,
                                 output_dir,
                                 gene_annotation_ref_file_path,
                                 TAD_file_path,
                                 genehancer_ref_file_path = NULL,
                                 chipseq_file_path = NULL,
                                 verbose = FALSE) {
    # import paths file  --------------------------------------------
    paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

    # get all required paths ---------------------------------
    sample_ids <- paths$sample_id
    sample_dir_paths <- file.path(output_dir, sample_ids)
    marker_summary_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".geneanno_marker_summary.Rds"))
    expression_file_paths <- file.path(sample_dir_paths, paste0(sample_ids, ".expression.Rds"))

    # check all required files ---------------------------------
    purrr::walk2(sample_dir_paths, sample_ids, function(s_dir, s_id) {
        check_sample_dir_existence(s_dir, s_id)
    })
    purrr::walk(marker_summary_file_paths, function(path) {
        check_file_existence_auto_created(path, name_of_file_type = "ASE marker summary .Rds file")
    })
    purrr::walk(expression_file_paths, function(path) {
        check_file_existence(path, name_of_file_type = "expression file")
    })

    cohort_marker_summary_list <- vector("list", length(sample_ids))
    cohort_expression_list <- vector("list", length(sample_ids))
    names(cohort_marker_summary_list) <- sample_ids
    names(cohort_expression_list) <- sample_ids


    # CREATE OHE REFERENCE ####################################
    # import required data ------------------------------------
    for (i in seq_len(length(sample_ids))) {
        sample_id <- sample_ids[i]
        sample_dir_path <- sample_dir_paths[i]
        marker_summary_file_path <- marker_summary_file_paths[i]
        expression_file_path <- expression_file_paths[i]

        cohort_marker_summary_list[[sample_id]] <- readRDS(marker_summary_file_path)
        cohort_expression_list[[sample_id]] <- readRDS(expression_file_path)
    }

    cohort_expression_table <- convert_cohort_expression_list_to_table(cohort_expression_list)
    rm(cohort_expression_list)
    cohort_marker_summary_table <- convert_marker_summary_list_to_table(cohort_marker_summary_list)
    rm(cohort_marker_summary_list)

    # process data  ------------------------------------------------

    # this function updates the cohort_expression_table by reference
    # no assignment of the returned value is needed
    generate_OHE_reference(cohort_expression_table, cohort_marker_summary_table, threshold_min_biallelic_samples = 10)

    # store results ------------------------------------------------
    saveRDS(
        cohort_expression_table,
        file = file.path(output_dir, "cohort.expression.reference.Rds")
    )

    data.table::fwrite(
        cohort_expression_table,
        file = file.path(output_dir, "cohort.expression.reference.txt"),
        sep = "\t"
    )

    log_msg("OHE reference generated", verbose = verbose)


    # PROCESS GENEHANCER ############################################
    # import required data  ---------------------------------------
    genehancer_raw_import <- import_genehancer_file(genehancer_ref_file_path)
    gene_annotation <- import_gene_annotation_ref_file(gene_annotation_ref_file_path)

    # process data  ------------------------------------------------
    genehancer_ref <- process_genehancer(genehancer_raw_import, gene_annotation)

    # store results ------------------------------------------------
    saveRDS(
        genehancer_ref,
        file.path(output_dir, paste0("genehancer_ref.Rds"))
    )

    readr::write_tsv(
        genehancer_ref,
        file.path(output_dir, paste0("genehancer_ref.txt"))
    )

    log_msg("GeneHancer reference generated", verbose = verbose)

    # PROCESS CHIPSEQ ############################################
    # import required data  ---------------------------------------
    chipseq_raw_import <- import_chipseq_file(chipseq_file_path)

    # process data  ------------------------------------------------
    chipseq <- process_chipseq(chipseq_raw_import)

    # store results ------------------------------------------------
    saveRDS(
        chipseq,
        file.path(output_dir, paste0("chipseq.Rds"))
    )

    readr::write_tsv(
        chipseq,
        file.path(output_dir, paste0("chipseq.txt"))
    )

    log_msg("ChIP-Seq data processed", verbose = verbose)


    # LINK CHIPSEQ TO TADs #######################################
    # import required data -----------------------------------
    TADs <- import_TAD_file(TAD_file_path)

    # process data  ------------------------------------------------
    gene_annotation_TADs <- add_TAD_boundaries_to_gene_annotation(gene_annotation, TADs)
    gene_annotation_TADs <- add_variant_overlap_window_to_gene_annotation(gene_annotation_TADs)
    chipseq_by_gene <- link_chipseq_to_genes(chipseq, gene_annotation_TADs)

    # store results ------------------------------------------------
    saveRDS(
        chipseq_by_gene,
        file.path(output_dir, paste0("chipseq_by_gene.Rds"))
    )

    readr::write_tsv(
        chipseq_by_gene,
        file.path(output_dir, paste0("chipseq_by_gene.txt"))
    )

    log_msg("ChIP-Seq data linked to genes", verbose = verbose)
}