#' Processing of Copy Number Data
#'
#' @param copy_numbers tibble of somatic copy number data
#' @param sample_id ID of the sample. Only used for adequate error messages
#' @importFrom magrittr %>%
process_copy_numbers <- function(copy_numbers, sample_id) {
    if(nrow(copy_numbers)> 1){
        # remove overlapping ends
        LEN <- length(copy_numbers$start)
        overlapping_ends <- (copy_numbers$end >= c(copy_numbers$start[2:LEN], Inf)) &
            (copy_numbers$chrom == c(copy_numbers$chrom[2:LEN], ""))
        cat(paste0("Removed ", sum(overlapping_ends), " overlapping ends\n"))
        # cut one base from end overlapping
        copy_numbers$end[which(overlapping_ends)] <- copy_numbers$end[which(overlapping_ends)] - 1

        # if there are_still_overlapping_ends THROW ERROR
        remaining_overlapping_ends <- (copy_numbers$end >= c(copy_numbers$start[2:LEN], Inf)) &
            (copy_numbers$chrom == c(copy_numbers$chrom[2:LEN], ""))

        if (sum(remaining_overlapping_ends) > 0) {
            cat(paste0(sum(remaining_overlapping_ends), " overlapping ends remain\n"))
            stop(paste0("copy number file of ", sample_id, " has overlapping ranges and cannot be used!"))
        }
    }
    return(copy_numbers)
}