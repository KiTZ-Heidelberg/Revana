#' FIMO Motif Id - TF Gene Name Conversion Table File
#'
#' @param annotation_file_path path to the annotation file from HOCOMOCO
#' @param output_path output path for th created reference file
#' @export
#' @importFrom magrittr %>%
prepare_motif_id_tf_gene_name_table_from_HOCOMOCO_annotation_tsv <- function(annotation_file_path, output_path) {
    readr::read_tsv(annotation_file_path) %>%
        dplyr::select(
            motif_id = Model,
            tf_gene_name = `Transcription Factor`
        ) %>%
        readr::write_tsv(output_path)
}