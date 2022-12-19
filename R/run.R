#' Run Revana on a cohort/subgroup of tumor sample data
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param gene_annotation_exons_ref_file_path Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param TAD_file_path Path to the TAD file. The TAD file contains genomic coordinates of Topologically Associated Domains (TADs). See the documentation for the exact format of this file.
#' @param fimo_motif_ref_path Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param motif_id_tf_gene_name_table_path Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param genehancer_ref_file_path Path to the reference for GeneHancer. See the documentation for the exact format of this file and how to obtain it.
#' @param chipseq_file_path Path to the ChIP-Seq file. The ChIP-Seq file contains genomic coordinates of regulatory active regions and can be obtained from experimental data or online less specific online resources. This file is optional, although an empty dummy file has to be used, if no ChIP-Seq data is supplied. See the documentation for the exact format of this file.
#' @param subgroup_name Name of the subgroup. The name of the subgroup is used as reference in later created HTML reports
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#' @param use_parallelization should parallelization be used to speed up Revana (TRUE / FALSE), defaults to TRUE
#'
#'
#' @export

run_analysis <- function(paths_file_path,
                output_dir,
                gene_annotation_ref_file_path,
                gene_annotation_exons_ref_file_path,
                TAD_file_path,
                fimo_motif_ref_path = NULL,
                motif_id_tf_gene_name_table_path = NULL,
                genehancer_ref_file_path = NULL,
                chipseq_file_path = NULL,
                subgroup_name = "DEFAULT_GROUP",
                run_tf_binding_site_analysis = FALSE,
                reference_genome = NULL,
                verbose = FALSE,
                use_parallelization = TRUE) {
    
    cat("\n")
    if(verbose == TRUE) {
        log_msg("ANALYSIS STARTED", log_time = verbose)
        cat("\n")
    }

    log_msg("Running Revana ...", log_time = FALSE)
    cat("\n")
    log_msg("Running Step 1 ...", log_time = FALSE)

    run_step1_per_cohort(
        paths_file_path,
        output_dir,
        gene_annotation_ref_file_path,
        gene_annotation_exons_ref_file_path,
        TAD_file_path,
        fimo_motif_ref_path,
        motif_id_tf_gene_name_table_path,
        genehancer_ref_file_path,
        chipseq_file_path,
        run_tf_binding_site_analysis,
        reference_genome
    )

    log_msg("Running Step 2 ...", log_time = FALSE)
    if(use_parallelization == TRUE){
        run_step2_per_cohort(
            paths_file_path,
            output_dir,
            gene_annotation_ref_file_path,
            gene_annotation_exons_ref_file_path,
            fimo_motif_ref_path,
            motif_id_tf_gene_name_table_path,
            run_tf_binding_site_analysis,
            reference_genome,
            verbose
        )
    }else{
        run_step2_per_cohort_serial(
            paths_file_path,
            output_dir,
            gene_annotation_ref_file_path,
            gene_annotation_exons_ref_file_path,
            fimo_motif_ref_path,
            motif_id_tf_gene_name_table_path,
            run_tf_binding_site_analysis,
            reference_genome,
            verbose
        )
    }
    

    log_msg("Running Step 3 ...", log_time = FALSE)
    run_step3_per_cohort(
        paths_file_path,
        output_dir,
        gene_annotation_ref_file_path,
        TAD_file_path,
        genehancer_ref_file_path,
        chipseq_file_path,
        verbose
    )

    log_msg("Running Step 4 ...", log_time = FALSE)
    run_step4_per_cohort(
        paths_file_path,
        output_dir,
        TAD_file_path,
        run_tf_binding_site_analysis,
        verbose,
        use_parallelization
    )

    log_msg("Running Step 5 ...", log_time = FALSE)
    run_step5_per_cohort(
        paths_file_path,
        output_dir,
        subgroup_name,
        run_tf_binding_site_analysis = run_tf_binding_site_analysis,
        verbose
    )

    cat("\n")
    log_msg("ANALYSIS COMPLETED", log_time = verbose)
    cat("\n")
}