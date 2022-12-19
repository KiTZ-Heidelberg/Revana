#' Merge RDS chunks containing data frames or data tables into one single RDS file
#'
#' @param pattern pattern to match the RDS chunk files
#' @param folder folder to search for RDS chunk files
#' @param output_file_path where to store the merged output file (full_path)
#' @param verbose should verbose logging be activated? (TRUE / FALSE)
#' @param max_size_all_chunks the maximum size of all the chunks combined that would still be merged
#' @export

merge_Rds_chunks <- function(pattern, folder, output_file_path, verbose = FALSE, max_size_all_chunks = 2000000000 ) {
    matched_files <- list.files(path = folder, pattern=pattern, recursive = FALSE, full.names = TRUE)

    number_of_files <- length(matched_files)
    if(number_of_files == 0){
        log_msg("Merging Rds chunks: 0 files matched", verbose = verbose, log_time = FALSE)
        return(NULL)
    }
    log_msg(paste0("Merging Rds chunks: ",number_of_files, " files matched"), verbose = verbose, log_time = FALSE)

    file_sizes <- matched_files %>% purrr::map(~file.info(.x)$size) %>% unlist()
    total_size <- sum(file_sizes)

    if(total_size > max_size_all_chunks) {
        log_msg(paste0("Merging Rds chunks: Merging aborted - File size of all chunks exceeded max size of ", max_size_all_chunks), verbose = verbose, log_time = FALSE)
    }else{
        imported_objects <- matched_files %>% purrr::map(~readRDS(.x))
        merged_object <- data.table::rbindlist(imported_objects)
        saveRDS(merged_object, file= output_file_path)
        log_msg(paste0("Merging Rds chunks: ",number_of_files, " files merged"), verbose = verbose, log_time = TRUE)
    }
}