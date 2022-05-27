#' Run Step 1 of the Revana workflow on a cohort/subgroup of tumor sample data
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param gene_annotation_exons_ref_file_path Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param fimo_motif_ref_path Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param motif_id_tf_gene_name_table_path Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param TAD_file_path Path to the TAD file. The TAD file contains genomic coordinates of Topologically Associated Domains (TADs). See the documentation for the exact format of this file.
#' @param genehancer_ref_file_path Path to the reference for GeneHancer. See the documentation for the exact format of this file and how to obtain it.
#' @param chipseq_file_path Path to the ChIP-Seq file. The ChIP-Seq file contains genomic coordinates of regulatory active regions and can be obtained from experimental data or online less specific online resources. This file is optional, although an empty dummy file has to be used, if no ChIP-Seq data is supplied. See the documentation for the exact format of this file.
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#'
#' @details
#' Step 1 of the revana workflow ...
#' * validates the paths file
#' * validates all reference files
#' * validates all sample files that are referenced in the paths file
#' * creates the results folder structure in the output directory
#'
#' @export
run_step1_per_cohort <- function(paths_file_path,
                                 output_dir,
                                 gene_annotation_ref_file_path,
                                 gene_annotation_exons_ref_file_path,
                                 TAD_file_path,
                                 fimo_motif_ref_path = NULL,
                                 motif_id_tf_gene_name_table_path = NULL,
                                 genehancer_ref_file_path = NULL,
                                 chipseq_file_path = NULL,
                                 run_tf_binding_site_analysis = FALSE,
                                 reference_genome= NULL) {
  # check directories ---------------------------------------------
  check_output_dir_existence(output_dir)

  # import paths file  --------------------------------------------
  paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

  # validate reference files  --------------------------------
  validate_reference_files(
    gene_annotation_ref_file_path,
    gene_annotation_exons_ref_file_path,
    fimo_motif_ref_path,
    motif_id_tf_gene_name_table_path,
    TAD_file_path,
    genehancer_ref_file_path,
    chipseq_file_path,
    run_tf_binding_site_analysis
  )

  # validate all referenced files  --------------------------------
  validate_files(paths)

  # validate dependencies  --------------------------------
  if(is.null(reference_genome)){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      stop(
        "Package \"BSgenome.Hsapiens.UCSC.hg19\" must be installed if you plan to run the TF Binding Site Analysis Feature.",
        call. = FALSE
      )
    }
  }
  

  # create all sample directories ---------------------------------
  sample_ids <- paths$sample_id
  sample_dir_paths <- file.path(output_dir, sample_ids)
  sample_dir_paths %>% purrr::walk(~ dir.create(.x))
}