calculate_cov_ratio_average_per_exon_unit <- function(copy_numbers, exon_units) {
    if(nrow(copy_numbers) > 0){
        exon_units_granges <- GenomicRanges::makeGRangesFromDataFrame(df = exon_units, seqnames.field = "chrom", start.field = "start", end.field = "end", ignore.strand = T)
        copy_numbers_granges <- GenomicRanges::makeGRangesFromDataFrame(df = copy_numbers, seqnames.field = "chrom", start.field = "start", end.field = "end", ignore.strand = T, keep.extra.columns = T)

        GenomeInfoDb::seqlevelsStyle(exon_units_granges) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(copy_numbers_granges) <- "UCSC"

        # ranges without copy number available -> copy number = 1
        GenomicRanges::setdiff(exon_units_granges, copy_numbers_granges) -> missing_copy_number_ranges
        if(length(missing_copy_number_ranges$seqnames) > 0 ){
            missing_copy_number_ranges$cov_ratio <- 1
            copy_numbers_granges_complete <- c(copy_numbers_granges, missing_copy_number_ranges)
        }else{
                copy_numbers_granges_complete <- copy_numbers_granges
        }

        # calculate binned average
        copy_numbers_RleList <- GenomicRanges::mcolAsRleList(copy_numbers_granges_complete, "cov_ratio")
        copy_numbers_per_exon_unit_granges <- GenomicRanges::binnedAverage(bins = exon_units_granges, numvar = copy_numbers_RleList, varname = "cov_ratio")
        exons_units_with_copy_number <- cbind(exon_units, data.frame(cov_ratio = copy_numbers_per_exon_unit_granges$cov_ratio))
    }else{
        exons_units_with_copy_number <- exon_units
        exons_units_with_copy_number$cov_ratio <- 1
    }
    return(exons_units_with_copy_number)
}

calculate_cov_ratio_average_per_gene <- function(exons_units_with_copy_number) {
    exons_units_with_copy_number %>%
        dplyr::mutate(length = (end - start)) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(length_all_exons = sum(length)) %>%
        dplyr::mutate(gene_fraction = (length / length_all_exons)) %>%
        dplyr::summarise(avg_cov_ratio = sum(cov_ratio * gene_fraction))
}