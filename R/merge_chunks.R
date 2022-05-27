#' Merge RDS chunks containing data frames or data tables into one single RDS file
#'
#' @param pattern pattern to match the RDS chunk files
#' @param folder folder to search for RDS chunk files
#' @param output_file_path where to store the merged output file (full_path)
#' @param max_size_all_chunks the maximum size of all the chunks combined that would still be merged
#' @export

merge_Rds_chunks <- function(pattern, folder, output_file_path, max_size_all_chunks = 2000000000) {
    matched_files <- list.files(path = folder, pattern=pattern, recursive = FALSE, full.names = TRUE)

    number_of_files <- length(matched_files)
    if(number_of_files == 0){
        cat("0 files matched.\n")
        return(NULL)
    }
    cat(paste0("matched ", number_of_files, " files\n"))

    file_sizes <- matched_files %>% purrr::map(~file.info(.x)$size) %>% unlist()
    total_size <- sum(file_sizes)

    if(total_size > max_size_all_chunks) {
        cat(paste0('File size of all chunks exceeded max size of ', 2000000000, "\n"))
    }else{
        imported_objects <- matched_files %>% purrr::map(~readRDS(.x))
        merged_object <- data.table::rbindlist(imported_objects)
        saveRDS(merged_object, file= output_file_path)
        cat("MERGING COMPLETED\n")
    }
}