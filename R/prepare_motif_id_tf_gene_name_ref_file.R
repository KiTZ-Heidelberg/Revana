prepare_motif_id_tf_gene_name_ref_file <- function (motif_id_tf_gene_name_file_path, motif_id_tf_gene_name_file_output_path) {
  # File import --------------------------
  motif_id_tf_gene_name_df <- readr::read_tsv(motif_id_tf_gene_name_file_path)

  # Filtering --------------------------
  motif_id_tf_gene_name_df_reformatted <- motif_id_tf_gene_name_df %>%
    # select columns
    dplyr::select(motif_id = Model, tf_gene_name = `Transcription factor`)

  # Write file to output -------------------
  readr::write_tsv(motif_id_tf_gene_name_df_reformatted, file = motif_id_tf_gene_name_file_output_path, col_names = T)
}
