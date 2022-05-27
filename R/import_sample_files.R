#' Imports provided marker file
#' @description The provided file should contain genomewide SNPs of both germline and somatic origin
#' @param path
#'
#' @return tibble
#' @noRd
import_marker_file <- function(path) {
    col_types <- "cicciiii"
    markers <- readr::read_tsv(file = path, col_types = col_types)
    return(markers)
}


#' Imports provided somatic SNV file
#' @description The provided file should contain genomewide SNVs and InDels of somatic origin
#' @param path
#'
#' @return tibble
#' @noRd
import_somatic_SNV_file <- function(path) {
    col_types <- "cicc"
    somatic_SNVs <- readr::read_tsv(file = path, col_types = col_types)
    return(somatic_SNVs)
}

#' Imports provided SV file
#' @description The provided file should contain genomewide SVs of somatic origin
#' @param path
#'
#' @return tibble
#' @noRd
import_SV_file <- function(path) {
    col_types <- "cicicc"
    SVs <- readr::read_tsv(file = path, col_types = col_types)
    return(SVs)
}

#' Imports provided CNA file
#' @description The provided file should contain genomewide Copy Number Aberrations
#' @param path
#'
#' @return tibble
#' @noRd
import_CNA_file <- function(path) {
    col_types <- "ciincn"
    CNAs <- readr::read_tsv(file = path, col_types = col_types)
    return(CNAs)
}

#' Imports provided copy_number file
#' @description The provided file should contain genomewide Copy Number Aberrations
#' @param path
#'
#' @return tibble
#' @noRd
import_copy_number_file <- function(path) {
    col_types <- "ciin"
    copy_numbers <- readr::read_tsv(file = path, col_types = col_types)
    return(copy_numbers)
}

#' Imports provided RNA expression file
#' @description The provided file should contain FPKM values per gene
#' @param path
#'
#' @return tibble
#' @noRd
import_expression_file <- function(path) {
    col_types <- "cd"
    expression <- readr::read_tsv(file = path, col_types = col_types)
    return(expression)
}