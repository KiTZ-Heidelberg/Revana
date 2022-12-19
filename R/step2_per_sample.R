#' Run Step 2 of the Revana workflow on a cohort/subgroup of tumor sample data - Parallelized implementation
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param sample_id ID of the sample
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param gene_annotation_exons_ref_file_path Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param fimo_motif_ref_path Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param motif_id_tf_gene_name_table_path Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#'
#' @details
#'
#' This function runs step 2 on a single sample.
#'
#' For each included sample step 2 of the revana workflow ...
#' * processes SNP marker input and calculates allele specific expression (ASE)
#' * processes somatic SNV/InDel input
#' * runs TF binding site analysis (if run_tf_binding_site_analysis == TRUE)
#' * processes copy number input and calculates gene copy numbers and coverage ratios
#' * processes expression input
#' * processes copy number aberration input (CNAs)
#' * processes structural variants (SVs)
#'
run_step2_per_sample <- function(paths_file_path,
                                 output_dir,
                                 sample_id,
                                 gene_annotation_ref_file_path,
                                 gene_annotation_exons_ref_file_path,
                                 fimo_motif_ref_path = NULL,
                                 motif_id_tf_gene_name_table_path = NULL,
                                 run_tf_binding_site_analysis = TRUE,
                                 reference_genome = NULL,
                                 verbose = FALSE) {

  log_msg("started step 2", verbose = verbose, sample_id = sample_id)

  # check directories --------------------------------------------
  check_output_dir_existence(output_dir)
  sample_dir <- file.path(output_dir, sample_id)
  check_sample_dir_existence(sample_dir, sample_id)

  # import paths file  --------------------------------------------
  paths <- import_paths_file(paths_file_path, check_file_table_headers = T)
  this_sample_id <- sample_id
  paths_this_sample <- paths %>% dplyr::filter(sample_id == this_sample_id)

  # validate all referenced files of this sample  ---------------
  validate_gene_annotation_ref_file_path(gene_annotation_ref_file_path)
  validate_files(paths_this_sample)
  if (nrow(paths_this_sample) > 1) {
    stop("More than one row exists for this sample_id")
  }


  # CALCULATE ASE ################################################

  # import required data  ---------------------------------------
  markers <- import_marker_file(paths_this_sample$marker_file)
  CNAs <- import_CNA_file(paths_this_sample$CNA_file)
  expression <- import_expression_file(paths_this_sample$expression_file)

  # import required references  ----------------------------------
  gene_annotation <- import_gene_annotation_ref_file(gene_annotation_ref_file_path)

  log_msg("ASE imports completed", verbose = verbose, sample_id = sample_id)

  # process data  ------------------------------------------------
  processed_markers <- process_markers(markers, CNAs)
  log_msg("markers processed", verbose = verbose, sample_id = sample_id)
  rm(markers)
  gc() # free memory

  runs <- find_marker_runs(processed_markers)
  log_msg("marker runs calculated", verbose = verbose, sample_id = sample_id)

  geneanno_marker_summary <- calculate_ASE(processed_markers, gene_annotation, runs)
  log_msg("ASE calculated", verbose = verbose, sample_id = sample_id)

  # store results ------------------------------------------------
  saveRDS(
    processed_markers,
    file.path(sample_dir, paste0(sample_id, ".markers.Rds"))
  )

  readr::write_tsv(
    processed_markers,
    file.path(sample_dir, paste0(sample_id, ".markers.txt"))
  )

  saveRDS(
    runs,
    file.path(sample_dir, paste0(sample_id, ".runs.Rds"))
  )

  readr::write_tsv(
    runs,
    file.path(sample_dir, paste0(sample_id, ".runs.txt"))
  )

  saveRDS(
    geneanno_marker_summary,
    file.path(sample_dir, paste0(sample_id, ".geneanno_marker_summary.Rds"))
  )

  readr::write_tsv(
    geneanno_marker_summary,
    file.path(sample_dir, paste0(sample_id, ".geneanno_marker_summary.txt"))
  )

  rm(processed_markers, runs, geneanno_marker_summary)
  gc() # free memory

  log_msg("ASE files written", verbose = verbose, sample_id = sample_id)

  # PROCESS SNVs ################################################

  # import required data  ---------------------------------------
  somatic_SNVs <- import_somatic_SNV_file(paths_this_sample$somatic_SNV_file)

  # store data --------------------------------------------------
  saveRDS(
    somatic_SNVs,
    file.path(sample_dir, paste0(sample_id, ".somatic_SNV.Rds"))
  )

  readr::write_tsv(
    somatic_SNVs,
    file.path(sample_dir, paste0(sample_id, ".somatic_SNV.txt"))
  )

  if (run_tf_binding_site_analysis) {

    log_msg("TF binding site analysis started", verbose = verbose, sample_id = sample_id)

    # import required references  ----------------------------------
    motif_id_tf_gene_name_table <- import_motif_id_tf_gene_name_ref_file(motif_id_tf_gene_name_table_path)

    # process data  ------------------------------------------------
    somatic_SNV_tf_binding_data <- process_somatic_SNVs_chunkwise(
      somatic_SNVs,
      expression,
      motif_id_tf_gene_name_table,
      fimo_in_file_path_prefix = file.path(sample_dir, paste0(sample_id, ".FIMO.in.")),
      fimo_out_file_path_prefix = file.path(sample_dir, paste0(sample_id, ".FIMO.out.")),
      fimo_log_file_path_prefix = file.path(sample_dir, paste0(sample_id, ".FIMO.log.")),
      fimo_in_file_path_suffix = ".txt",
      fimo_out_file_path_suffix = ".txt",
      fimo_log_file_path_suffix = ".txt",
      fimo_motif_ref_path = fimo_motif_ref_path,
      SNV_chunk_size = 3000,
      delete_FIMO_out_right_away = TRUE,
      reference_genome = reference_genome
    )

    # store results ------------------------------------------------
    saveRDS(
      somatic_SNV_tf_binding_data,
      file.path(sample_dir, paste0(sample_id, ".somatic_SNV_tf_binding_data.Rds"))
    )

    readr::write_tsv(
      somatic_SNV_tf_binding_data,
      file.path(sample_dir, paste0(sample_id, ".somatic_SNV_tf_binding_data.txt"))
    )

    rm(motif_id_tf_gene_name_table, somatic_SNV_tf_binding_data)
    gc() # free memory

    log_msg("TF binding site analysis completed", verbose = verbose, sample_id = sample_id)
  }

  rm(somatic_SNVs)
  gc() # free memory

  log_msg("somatic SNVs processed", verbose = verbose, sample_id = sample_id)


  # PROCESS COPY NUMBERS ########################################
  # import required data  ---------------------------------------
  copy_numbers <- import_copy_number_file(paths_this_sample$copy_number_file)


  # import required references  ----------------------------------
  exons <- import_gene_annotation_exon_ref_file(gene_annotation_exons_ref_file_path)
  exon_units <- get_non_overlapping_exon_units(exons)

  # process data  ------------------------------------------------
  processed_copy_numbers <- process_copy_numbers(copy_numbers, sample_id)
  copy_number_by_gene <- calculate_cov_ratio_average_per_exon_unit(processed_copy_numbers, exon_units) %>%
    calculate_cov_ratio_average_per_gene()

  # store results ------------------------------------------------
  saveRDS(
    processed_copy_numbers,
    file.path(sample_dir, paste0(sample_id, ".copy_number.Rds"))
  )

  readr::write_tsv(
    processed_copy_numbers,
    file.path(sample_dir, paste0(sample_id, ".copy_number.txt"))
  )

  saveRDS(
    copy_number_by_gene,
    file.path(sample_dir, paste0(sample_id, ".copy_number_by_gene.Rds"))
  )

  readr::write_tsv(
    copy_number_by_gene,
    file.path(sample_dir, paste0(sample_id, ".copy_number_by_gene.txt"))
  )

  log_msg("copy numbers processed", verbose = verbose, sample_id = sample_id)



  # PROCESS EXPRESSION ##########################################

  # process data  ------------------------------------------------
  expression_processed <- process_expression(
    expression,
    gene_annotation,
    copy_number_by_gene
  )

  # store results ------------------------------------------------
  saveRDS(
    expression_processed,
    file.path(sample_dir, paste0(sample_id, ".expression.Rds"))
  )

  readr::write_tsv(
    expression_processed,
    file.path(sample_dir, paste0(sample_id, ".expression.txt"))
  )

  rm(expression, expression_processed, gene_annotation, copy_numbers, processed_copy_numbers, copy_number_by_gene)
  gc() # free memory

  log_msg("expression processed", verbose = verbose, sample_id = sample_id)


  # PROCESS CNAs ##########################################

  # process data  ------------------------------------------------
  CNAs_processed <- process_CNAs(CNAs)

  # store results ------------------------------------------------
  saveRDS(
    CNAs_processed,
    file.path(sample_dir, paste0(sample_id, ".CNA.Rds"))
  )

  readr::write_tsv(
    CNAs_processed,
    file.path(sample_dir, paste0(sample_id, ".CNA.txt"))
  )

  rm(CNAs, CNAs_processed)
  gc() # free memory

  log_msg("CNAs processed", verbose = verbose, sample_id = sample_id)



  # PROCESS SVs ##########################################

  # import required data  ---------------------------------------
  SVs <- import_SV_file(paths_this_sample$SV_file)

  # process data  ------------------------------------------------
  SVs_processed <- process_SVs(SVs)

  # store results ------------------------------------------------
  saveRDS(
    SVs_processed,
    file.path(sample_dir, paste0(sample_id, ".SV.Rds"))
  )

  readr::write_tsv(
    SVs_processed,
    file.path(sample_dir, paste0(sample_id, ".SV.txt"))
  )

  rm(SVs_processed, SVs)
  gc() # free memory

  log_msg("SVs processed", verbose = verbose, sample_id = sample_id)

  log_msg("completed step 2", verbose = verbose, sample_id = sample_id)
}


#' Run Step 2 of the Revana workflow on a cohort/subgroup of tumor sample data - Parallelized implementation
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param gene_annotation_exons_ref_file_path Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param fimo_motif_ref_path Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param motif_id_tf_gene_name_table_path Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#'
#'
#' @details
#'
#' This implementation of step 2 runs the samples in parallel. This is the optimized and default implementation.
#'
#' For each included sample step 2 of the revana workflow ...
#' * processes SNP marker input and calculates allele specific expression (ASE)
#' * processes somatic SNV/InDel input
#' * runs TF binding site analysis (if run_tf_binding_site_analysis == TRUE)
#' * processes copy number input and calculates gene copy numbers and coverage ratios
#' * processes expression input
#' * processes copy number aberration input (CNAs)
#' * processes structural variants (SVs)
#'
#' @export

run_step2_per_cohort <- function(paths_file_path,
                                 output_dir,
                                 gene_annotation_ref_file_path,
                                 gene_annotation_exons_ref_file_path,
                                 fimo_motif_ref_path = NULL,
                                 motif_id_tf_gene_name_table_path = NULL,
                                 run_tf_binding_site_analysis = TRUE,
                                 reference_genome = NULL,
                                 verbose = FALSE) {
  # import paths file  --------------------------------------------
  paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

  n_samples <- nrow(paths)

  if(.Platform$OS.type == "unix") {
    future::plan(future::multicore)
  }else{
    future::plan(future::multisession)
  }

  cat("Step 2: running in parallel on ")
  cat(future::nbrOfWorkers())
  cat(" workers\n")

  progressr::with_progress({
    p <- progressr::progressor(steps = n_samples)
    unique(paths$sample_id) %>% furrr::future_walk(function(x) {
      # uncomment for easier debugging:
      # withCallingHandlers(message=handle_message, warning = handle_warning, {
      run_step2_per_sample(
        paths_file_path,
        output_dir,
        x,
        gene_annotation_ref_file_path,
        gene_annotation_exons_ref_file_path,
        fimo_motif_ref_path,
        motif_id_tf_gene_name_table_path,
        run_tf_binding_site_analysis,
        reference_genome,
        verbose
      )
      # })

      p() # progress ++
    })
  })
}

#' Run Step 2 of the Revana workflow on a cohort/subgroup of tumor sample data - Serial implementation
#'
#' @param paths_file_path Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.
#' @param output_dir Where are the analysis results to be stored
#' @param gene_annotation_ref_file_path Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param gene_annotation_exons_ref_file_path Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.
#' @param fimo_motif_ref_path Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param motif_id_tf_gene_name_table_path Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.
#' @param run_tf_binding_site_analysis Should TF binding site analysis be conducted?
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#'
#'
#' @details
#'
#' This implementation of step 2 runs the samples serially. This is not the default implementation.
#'
#' For each included sample step 2 of the revana workflow ...
#' * processes SNP marker input and calculates allele specific expression (ASE)
#' * processes somatic SNV/InDel input
#' * runs TF binding site analysis (if run_tf_binding_site_analysis == TRUE)
#' * processes copy number input and calculates gene copy numbers and coverage ratios
#' * processes expression input
#' * processes copy number aberration input (CNAs)
#' * processes structural variants (SVs)
#'
#' @export

run_step2_per_cohort_serial <- function(paths_file_path,
                                        output_dir,
                                        gene_annotation_ref_file_path,
                                        gene_annotation_exons_ref_file_path,
                                        fimo_motif_ref_path = NULL,
                                        motif_id_tf_gene_name_table_path = NULL,
                                        run_tf_binding_site_analysis = TRUE,
                                        reference_genome = NULL,
                                        verbose = FALSE) {
  # import paths file  --------------------------------------------
  paths <- import_paths_file(paths_file_path, check_file_table_headers = T)

  n_samples <- nrow(paths)

  cat("Step 2: running serially\n")

  progressr::with_progress({
    p <- progressr::progressor(steps = n_samples)
    for (s_id in paths$sample_id) {
      run_step2_per_sample(
        paths_file_path,
        output_dir,
        s_id,
        gene_annotation_ref_file_path,
        gene_annotation_exons_ref_file_path,
        fimo_motif_ref_path,
        motif_id_tf_gene_name_table_path,
        run_tf_binding_site_analysis,
        reference_genome,
        verbose
      )
      p()
    }
  })
}