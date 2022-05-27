#' Imports provided ChIP-Seq file
#' @description The provided file should contain ChIP-Seq data
#' @param path
#'
#' @return tibble
#' @noRd
import_chipseq_file <- function(path) {
    if (is.null(path)) {
        chipseq <- data.frame(
            chrom = character(0),
            start = integer(0),
            end = integer(0),
            cluster = character(0)
        )
    } else {
        col_types <- "ciic"
        chipseq <- readr::read_tsv(file = path, col_types = col_types)
    }
    return(chipseq)
}

#' Imports reference gene annotation
#' @description The provided file should contain a subset of columns of the GENCODE gene annotation
#' @param path
#'
#' @return tibble
#' @noRd
import_gene_annotation_ref_file <- function(path) {
    col_types <- "ciiiccccccl"
    gene_annotation_ref <- readr::read_tsv(file = path, col_types = col_types)
    return(gene_annotation_ref)
}

#' Imports reference exon annotation
#' @description The provided file should contain a subset of columns of the GENCODE gene annotation
#' @param path
#'
#' @return tibble
#' @noRd
import_gene_annotation_exon_ref_file <- function(path) {
    col_types <- "ciic"
    gene_annotation_ref_exon <- readr::read_tsv(file = path, col_types = col_types)
    return(gene_annotation_ref_exon)
}

#' Imports reference gene annotation
#' @description The provided file should contain a subset of columns of the GENCODE gene annotation
#' @param path
#'
#' @return tibble
#' @noRd
import_genehancer_file <- function(path) {
    if (is.null(path)) {
        genehancer <- data.frame(
            chrom = character(0),
            feature_name = character(0),
            start = integer(0),
            end = integer(0),
            score = numeric(0),
            genehancer_id = character(0),
            connected_gene = character(0),
            connected_gene_score = numeric(0),
            is_elite = logical(0),
            is_association_elite = logical(0)
        )
    } else {
        col_types <- "cciinccnll"
        genehancer <- readr::read_tsv(file = path, col_types = col_types)
    }

    return(genehancer)
}

#' Imports FIMO motif id to gene name translation table
#' @description The provided file should contain motif ids and matching gene names
#' @param path
#'
#' @return tibble
#' @noRd
import_motif_id_tf_gene_name_ref_file <- function(path) {
    col_types <- "cc"
    motif_id_tf_gene_name <- readr::read_tsv(file = path, col_types = col_types)
    return(motif_id_tf_gene_name)
}

#' Imports provided TAD file
#' @description The provided file should contain TAD data (Topologically Associated Domains)
#' @param path
#'
#' @return tibble
#' @noRd
import_TAD_file <- function(path) {
    col_types <- "cii"
    TADs <- readr::read_tsv(file = path, col_types = col_types)
    return(TADs)
}