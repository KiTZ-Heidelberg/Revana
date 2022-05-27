calculate_ref_length <- function(somatic_SNVs) {
  somatic_SNVs_with_ref_length <- somatic_SNVs %>% dplyr::mutate(
    ref_length = stringr::str_length(ref)
  )
  return(somatic_SNVs_with_ref_length)
}

mark_SNVs_close_to_adjacent <- function(somatic_SNVs, sequence_padding = 20) {
  somatic_SNVs_marked <- somatic_SNVs %>%
    dplyr::mutate(
      pos_diff_to_prev = dplyr::if_else(chrom == dplyr::lag(chrom, n = 1), pos - dplyr::lag(pos, n = 1), NA_integer_),
      pos_diff_to_next = dplyr::if_else(chrom == dplyr::lead(chrom, n = 1), dplyr::lead(pos, n = 1) - pos, NA_integer_)
    ) %>%
    dplyr::mutate(
      close_to_prev = dplyr::if_else(is.na(pos_diff_to_prev), F, pos_diff_to_prev <= (sequence_padding - 1)),
      close_to_next = dplyr::if_else(is.na(pos_diff_to_next), F, pos_diff_to_next <= (sequence_padding - 1))
    ) %>%
    dplyr::select(!pos_diff_to_prev) %>%
    dplyr::select(!pos_diff_to_next)
  return(somatic_SNVs_marked)
}

add_sequence_to_somatic_SNVs <- function(somatic_SNVs, reference_genome = NULL, sequence_padding = 20) {
  somatic_snv_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
    somatic_SNVs,
    seqnames.field = "chrom",
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = T
  )
  somatic_SNVs_GRanges_before <- somatic_snv_GRanges
  GenomicRanges::start(somatic_SNVs_GRanges_before) <- GenomicRanges::start(somatic_SNVs_GRanges_before) - sequence_padding
  GenomicRanges::end(somatic_SNVs_GRanges_before) <- GenomicRanges::end(somatic_SNVs_GRanges_before) - 1


  somatic_SNVs_GRanges_after <- somatic_snv_GRanges
  GenomicRanges::end(somatic_SNVs_GRanges_after) <- GenomicRanges::end(somatic_SNVs_GRanges_after) + somatic_snv_GRanges$ref_length + sequence_padding
  GenomicRanges::start(somatic_SNVs_GRanges_after) <- GenomicRanges::start(somatic_SNVs_GRanges_after) + somatic_snv_GRanges$ref_length

  # set common seqlevelStyle: UCSC => "chr1"
  GenomeInfoDb::seqlevelsStyle(somatic_SNVs_GRanges_before) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(somatic_SNVs_GRanges_after) <- "UCSC"

  # define reference genome
  if (is.null(reference_genome)) {
    actual_reference_genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else {
    actual_reference_genome <- reference_genome
  }


  # get sequences

  seq_before_snv <- as.character(Biostrings::getSeq(actual_reference_genome, somatic_SNVs_GRanges_before))
  seq_after_snv <- as.character(Biostrings::getSeq(actual_reference_genome, somatic_SNVs_GRanges_after))

  somatic_SNVs_with_seq <- somatic_SNVs %>%
    dplyr::mutate(
      seq_ref = paste0(seq_before_snv, ref, seq_after_snv),
      seq_alt = paste0(seq_before_snv, alt, seq_after_snv)
    ) %>%
    dplyr::mutate(
      is_seq_potentially_incorrect = (close_to_prev | close_to_next)
    )

  # NOTE:   An implementation of sequence generation, that includes adjacent SNVs
  #         would provide further value and could be implemented in the future.
  #         Here it would also be necessary, to consider the WGS allele frequency in
  #         order to prevent matching SNVs that lie on different alleles.

  return(somatic_SNVs_with_seq)
}

run_fimo <- function(somatic_SNVs_with_seq,
                     fimo_in_file_path,
                     fimo_out_file_path,
                     fimo_log_file_path,
                     fimo_motif_ref_path,
                     selected_motifs) {
  # open file in write mode
  fimo_in_file_con <- file(fimo_in_file_path, open = "w")

  # reformat somatic SNVs for FIMO and write to FIMO IN FILE
  somatic_SNVs_with_seq %>%
    purrr::pwalk(
      function(chrom, pos, ref, alt, ref_length, close_to_prev, close_to_next, seq_ref, seq_alt, is_seq_potentially_incorrect) {
        ref_id <- paste0(chrom, ".", format(pos, scientific = F), ".", ref, ".", alt, ".ref")
        mut_id <- paste0(chrom, ".", format(pos, scientific = F), ".", ref, ".", alt, ".mut")

        writeLines(paste0(">", ref_id, "\n"), con = fimo_in_file_con)
        writeLines(paste0(seq_ref, "\n"), con = fimo_in_file_con)
        writeLines(paste0(">", mut_id, "\n"), con = fimo_in_file_con)
        writeLines(paste0(seq_alt, "\n"), con = fimo_in_file_con)
      }
    )
  close(fimo_in_file_con)

  # FIMO MOTIFS COMMAND OPTION
  fimo_motifs_command_option <- selected_motifs %>%
    purrr::map(~ paste0("--motif \"", .x, "\"")) %>%
    unlist() %>%
    paste0(collapse = " ")

  # RUN FIMO
  run_fimo_command <- paste0(
    "fimo --verbosity 1 --thresh 1 --no-qvalue --text ",
    fimo_motifs_command_option,
    " ",
    fimo_motif_ref_path,
    " ",
    fimo_in_file_path,
    " > ",
    fimo_out_file_path,
    " 2>> ",
    fimo_log_file_path
  )

  system(run_fimo_command)
  system(paste0("echo >> ", fimo_out_file_path))
}


add_tf_gene_name <- function(fimo_data_tidy, motif_id_tf_gene_name_table) {
  motif_id_tf_gene_name_df <- motif_id_tf_gene_name_table %>%
    dplyr::select(motif_id, tf_gene_name)
  fimo_data_tidy_with_tf_gene_name <- fimo_data_tidy %>%
    dplyr::left_join(motif_id_tf_gene_name_df, by = c("motif_id"))
  return(fimo_data_tidy_with_tf_gene_name)
}

add_tf_expression <- function(fimo_data_tidy_with_tf_gene_name, expression_data) {
  tf_expression <- expression_data %>%
    dplyr::select(tf_gene_name = gene_name, tf_FPKM = FPKM)
  fimo_data_tidy_with_tf_expression <- fimo_data_tidy_with_tf_gene_name %>%
    dplyr::left_join(tf_expression, by = c("tf_gene_name"))
  return(fimo_data_tidy_with_tf_expression)
}

evaluate_tf_binding_sites <- function(fimo_data_tidy_with_tf_expression, thres_tf_FPKM = 10) {
  data_tf_binding_sites_called <- fimo_data_tidy_with_tf_expression %>%
    dplyr::mutate(relevant_tf_binding_site = (
      new_significant_binding_site & (tf_FPKM >= thres_tf_FPKM)
    ))
  return(data_tf_binding_sites_called)
}

post_fimo_analysis_data_table <- function(fimo_out_file_path, sequence_padding = 20, p_value_significance_threshold = 0.05) {
  fimo_data <- data.table::fread(fimo_out_file_path)
  # delete empty q-value column
  fimo_data[, `q-value` := NULL]
  # separate SNV4
  fimo_data[, c("chrom", "pos", "ref", "alt", "mut.or.ref") := data.table::tstrsplit(sequence_name, ".", fixed = TRUE)]
  # get mut.or.ref.length
  fimo_data[, mut.or.ref.length := ifelse(mut.or.ref == "mut", stringr::str_length(alt), stringr::str_length(ref))]

  fimo_data_filtered <- fimo_data[(start <= (sequence_padding + mut.or.ref.length)) & (stop >= (sequence_padding + 1))]
  # only min p-values
  fimo_data_filtered_wide <-
    fimo_data_filtered[fimo_data_filtered[, .I[which.min(`p-value`)], by = c("sequence_name", "motif_id")]$V1] %>%
    # # only min p-values
    # fimo_data_filtered_wide <-
    #   fimo_data[fimo_data[, .I[which.min(`p-value`)], by=c("sequence_name","motif_id")]$V1
    #             ][(start <= (sequence_padding + mut.or.ref.length))&(stop >= (sequence_padding+1))] %>%

    # to wide format
    data.table::dcast(motif_id + chrom + pos + ref + alt ~ mut.or.ref, value.var = c("mut.or.ref.length", "score", "p-value", "matched_sequence"))

  # release memory
  rm(fimo_data)
  gc()


  # calculate transcription factor binding differences
  fimo_data_filtered_wide[, score_diff := score_mut - score_ref]
  fimo_data_filtered_wide[, log.diff.p.values := log10(`p-value_ref`) - log10(`p-value_mut`)]
  fimo_data_filtered_wide[, `p-value_mut_significant` := `p-value_mut` <= p_value_significance_threshold]
  fimo_data_filtered_wide[, `p-value_ref_significant` := `p-value_ref` <= p_value_significance_threshold]
  fimo_data_filtered_wide[, new_significant_binding_site := (`p-value_mut_significant` & (!`p-value_ref_significant`))]

  return(data.table::setDT(fimo_data_filtered_wide))
}


select_TF_motifs <- function(motif_id_tf_gene_name_table, expression, FPKM_threshold) {
  motif_id_tf_gene_name_table %>%
    add_tf_expression(expression) %>%
    dplyr::filter(tf_FPKM >= FPKM_threshold) %>%
    dplyr::pull(motif_id)
}


#' Processessing of somatic SNVs - chunkwise implementation
#'
#' @param somatic_SNVs tibble of somatic SNVs
#' @param expression tibble of FPKM expression data
#' @param motif_id_tf_gene_name_table translation table from FIMO motifs to transcription factor gene names
#' @param fimo_in_file_path_prefix path to where should the FIMO input file be saved. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_out_file_path_prefix path to where  FIMO should save its output file. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_log_file_path_prefix path to where FIMO should save its log file. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_in_file_path_suffix path to where should the FIMO input file be saved. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_out_file_path_suffix path to where  FIMO should save its output file. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_log_file_path_suffix path to where FIMO should save its log file. Path is assembled like this PREFIX_CHUNK-NUMBER_SUFFIX
#' @param fimo_motif_ref_path path to the FIMO motif reference file
#' @param SNV_chunk_size input size of the SNV chunk (number of SNV rows to be used for each chunk)
#' @param delete_FIMO_out_right_away should the FIMO data be deleted right after using it? TRUE: Eliminates lots of unnecessary output data 
#' @param reference_genome Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)
#' @importFrom magrittr %>%
#' @importFrom data.table :=
process_somatic_SNVs_chunkwise <- function(somatic_SNVs,
                                           expression,
                                           motif_id_tf_gene_name_table,
                                           fimo_in_file_path_prefix,
                                           fimo_out_file_path_prefix,
                                           fimo_log_file_path_prefix,
                                           fimo_in_file_path_suffix,
                                           fimo_out_file_path_suffix,
                                           fimo_log_file_path_suffix,
                                           fimo_motif_ref_path,
                                           SNV_chunk_size = 10000,
                                           delete_FIMO_out_right_away = TRUE,
                                           reference_genome = NULL) {

  # Select motifs for tf that are sufficiently expressed in sample
  selected_motifs <- select_TF_motifs(motif_id_tf_gene_name_table, expression, FPKM_threshold = 5)


  n_somatic_SNV_rows <- nrow(somatic_SNVs)
  somatic_SNVs_chunked_list <- split(somatic_SNVs, ((seq(nrow(somatic_SNVs)) - 1) %/% SNV_chunk_size) + 1)
  max_chunk_index <- ((n_somatic_SNV_rows - 1) %/% SNV_chunk_size) + 1

  tf_binding_analysis_data_list <- vector("list", max_chunk_index)

  for (i in seq_len(max_chunk_index)) {
    gc()

    somatic_SNVs_of_chunk <- somatic_SNVs_chunked_list[[i]]
    fimo_in_file_path <- paste0(fimo_in_file_path_prefix, "_", i, "_", fimo_in_file_path_suffix)
    fimo_out_file_path <- paste0(fimo_out_file_path_prefix, "_", i, "_", fimo_out_file_path_suffix)
    fimo_log_file_path <- paste0(fimo_log_file_path_prefix, "_", i, "_", fimo_log_file_path_suffix)


    somatic_SNVs_of_chunk %>%
      calculate_ref_length() %>%
      mark_SNVs_close_to_adjacent() %>%
      add_sequence_to_somatic_SNVs(reference_genome) %>%
      run_fimo(
        fimo_in_file_path,
        fimo_out_file_path,
        fimo_log_file_path,
        fimo_motif_ref_path,
        selected_motifs
      )
    tf_binding_analysis_data_of_chunk <-
      post_fimo_analysis_data_table(fimo_out_file_path) %>%
      add_tf_gene_name(motif_id_tf_gene_name_table) %>%
      add_tf_expression(expression) %>%
      evaluate_tf_binding_sites()

    gc()
    if (delete_FIMO_out_right_away == TRUE) {
      file.remove(fimo_out_file_path)
    }
    tf_binding_analysis_data_list[[i]] <- tf_binding_analysis_data_of_chunk
  }

  tf_binding_analysis_data_table <- data.table::rbindlist(tf_binding_analysis_data_list)

  return(tf_binding_analysis_data_table)
}