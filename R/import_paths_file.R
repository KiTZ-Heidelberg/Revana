#' Import provided paths file
#'
#' @param paths_file_path path to the paths_file
#' @param check_file_table_headers logical value that whether the path file should be checked for containing correct headers. Defaults to TRUE
#' @param col_types overwrite the column types of the imported paths file.
#' 
#'
#' @return tibble with the provided paths of all samples
#' @export
#'
#'
import_paths_file <- function(paths_file_path, check_file_table_headers = T, col_types = "ccccccc") {
  # check if paths file exists
  check_file_existence(paths_file_path, name_of_file_type = "paths file")

  # check paths file header
  required_cols <- c(
    "sample_id",
    "marker_file",
    "CNA_file",
    "somatic_SNV_file",
    "expression_file",
    "SV_file",
    "copy_number_file"
  )

  if (check_file_table_headers == T) {
    check_file_table_header(paths_file_path, required_cols = required_cols, name_of_file_type = "paths file")
  }

  # import paths file
  paths <- readr::read_tsv(
    file = paths_file_path,
    col_types = col_types
  )
  return(paths)
}