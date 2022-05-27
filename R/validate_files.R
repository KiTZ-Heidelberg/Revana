validate_files <- function(paths_tibble) {
  # check for sample duplicates --------------
  duplicate_samples <- paths_tibble %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(sample_id)

  if (length(duplicate_samples) > 0) {
    stop(
      paste0(
        "\n\tDuplicate sample_ids were provided:\n\t",
        paste("\t", duplicate_samples, sep = "", collapse = "\n"),
        "\n"
      )
    )
  }


  # check for existence ----------------------
  # marker file
  check_referenced_marker_path <- function(path) check_file_existence(path, name_of_file_type = "marker file")
  purrr::map(paths_tibble$marker_file, check_referenced_marker_path)

  # CNA file
  check_referenced_CNA_path <- function(path) check_file_existence(path, name_of_file_type = "CNA file")
  purrr::map(paths_tibble$CNA_file, check_referenced_CNA_path)

  # copy number file
  check_referenced_copy_number_path <- function(path) check_file_existence(path, name_of_file_type = "copy_number file")
  purrr::map(paths_tibble$copy_number_file, check_referenced_copy_number_path)

  # somatic SNVs file
  check_referenced_somatic_SNV_path <- function(path) check_file_existence(path, name_of_file_type = "somatic SNV file")
  purrr::map(paths_tibble$somatic_SNV_file, check_referenced_somatic_SNV_path)

  # expression file
  check_referenced_expression_path <- function(path) check_file_existence(path, name_of_file_type = "expression file")
  purrr::map(paths_tibble$expression_file, check_referenced_expression_path)

  # SV file
  check_referenced_SV_path <- function(path) check_file_existence(path, name_of_file_type = "SV file")
  purrr::map(paths_tibble$SV_file, check_referenced_SV_path)

  # check for right headers ------------------
  # marker file
  required_cols_marker_file <- c("chrom", "pos", "ref", "alt", "reads.WGS.ref", "reads.WGS.alt", "reads.RNA.ref", "reads.RNA.alt")
  check_referenced_marker_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_marker_file, name_of_file_type = "marker file")
  purrr::map(paths_tibble$marker_file, check_referenced_marker_file_header)

  # CNA file
  required_cols_CNA_file <- c("chrom", "start", "end", "copy_number", "CNA_type", "log2")
  check_referenced_CNA_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_CNA_file, name_of_file_type = "CNA file")
  purrr::map(paths_tibble$CNA_file, check_referenced_CNA_file_header)

  # copy number file
  required_cols_copy_number_file <- c("chrom", "start", "end", "cov_ratio")
  check_referenced_copy_number_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_copy_number_file, name_of_file_type = "copy number file")
  purrr::map(paths_tibble$copy_number_file, check_referenced_copy_number_file_header)

  # somatic SNVs file
  required_cols_somatic_SNV_file <- c("chrom", "pos", "ref", "alt")
  check_referenced_somatic_SNV_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_somatic_SNV_file, name_of_file_type = "somatic SNV file")
  purrr::map(paths_tibble$somatic_SNV_file, check_referenced_somatic_SNV_file_header)

  # expression file
  required_cols_expression_file <- c("gene_name", "FPKM")
  check_referenced_expression_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_expression_file, name_of_file_type = "expression file")
  purrr::map(paths_tibble$expression_file, check_referenced_expression_file_header)

  # SV file
  required_cols_SV_file <- c("chrom1", "pos1", "chrom2", "pos2", "sv_type", "eventInversion")
  check_referenced_SV_file_header <- function(path) check_file_table_header(path, required_cols = required_cols_SV_file, name_of_file_type = "SV file")
  purrr::map(paths_tibble$SV_file, check_referenced_SV_file_header)
}