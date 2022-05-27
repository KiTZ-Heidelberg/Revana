get_non_overlapping_exon_units <- function(exons) {
    exons_granges_list <- GenomicRanges::makeGRangesListFromDataFrame(exons, split.field = "gene_name")

    # reduce to unique non-overlapping units
    non_overlapping_exons_granges_list <- GenomicRanges::reduce(exons_granges_list)

    # convert back to df
    non_overlapping_exons_grange <- IRanges::stack(non_overlapping_exons_granges_list, index.var = "gene_name")
    non_overlapping_exons <- data.frame(
        chrom = GenomicRanges::seqnames(non_overlapping_exons_grange),
        start = GenomicRanges::start(non_overlapping_exons_grange),
        end = GenomicRanges::end(non_overlapping_exons_grange),
        gene_name = non_overlapping_exons_grange$gene_name
    )
    return(non_overlapping_exons)
}