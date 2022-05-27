process_expression <- function(expression, gene_annotation, copy_number_by_gene) {
    expression_unique_by_gene_name <- expression %>%
        # eliminate duplicates
        # keep max FPKM
        dplyr::group_by(gene_name) %>%
        dplyr::slice_max(order_by = FPKM, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    gene_annotation_gene_names <- gene_annotation %>% dplyr::select(gene_name)
    expression_of_geneanno_genes <- gene_annotation_gene_names %>%
        dplyr::left_join(expression_unique_by_gene_name, by = c("gene_name")) %>%
        dplyr::left_join(copy_number_by_gene, by = c("gene_name")) %>%
        dplyr::mutate(FPKM_copy_number_normalized = FPKM / avg_cov_ratio)
    return(expression_of_geneanno_genes)
}