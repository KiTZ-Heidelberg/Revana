#' Combine results paths file from several analysed subgroups into one
#'
#' @param list_of_results_paths_files A list of the results paths files for each analysed subcohort. Each file is created as output paths summary by revana::run()
#' @param output_path_merged_results_file Path, where the new combined results paths file is to be stored
#'
#'
#' @export


merge_results_paths_files <- function(list_of_results_paths_files, output_path_merged_results_file) {
    if (!is.list(list_of_results_paths_files)) {
        stop("merge_results_paths_files requires argument list_of_results_paths_files to be a list\n")
    }

    list_of_results_paths_files %>%
        purrr::map(readr::read_tsv) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        readr::write_tsv(file = output_path_merged_results_file)
}