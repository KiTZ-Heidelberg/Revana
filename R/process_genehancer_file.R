#' Processing of GeneHancer Reference File
#'
#' @param genehancer tibble of GeneHancers
#' @param gene_annotation tibble of GeneAnnotation
#' @importFrom magrittr %>%

process_genehancer <- function(genehancer, gene_annotation) {
    genehancer_of_geneanno_genes <- genehancer %>%
        dplyr::semi_join(
            gene_annotation,
            by = c("connected_gene" = "gene_name")
        )
    return(genehancer_of_geneanno_genes)
}