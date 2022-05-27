add_chipseq_info_to_SNVs <- function (SNVs, chipseq) {
    data.table::setDT(SNVs)
    data.table::setDT(chipseq)

    if(nrow(SNVs) == 0) {
        # returns empty table if SNVs are empty
        SNVs[, overlaps_with_chipseq := character(0)]

        return (SNVs)
    }

    if(nrow(chipseq) == 0) {
        # returns FALSE column if chipseq is empty
        SNVs[, overlaps_with_chipseq := FALSE]

        return (SNVs)
    }

    # convert data tables to GRanges
    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq,
        seqnames.field="chipseq_chrom",
        start.field="chipseq_start",
        end.field="chipseq_end",
        ignore.strand=TRUE
    )

    SNV_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        SNVs,
        seqnames.field="chrom",
        start.field="pos",
        end.field="pos",
        ignore.strand=TRUE
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(chipseq_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SNV_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    n_overlaps <- GenomicRanges::countOverlaps(query = SNV_GRanges, subject = chipseq_GRanges)
    any_overlaps <- n_overlaps > 0

    SNVs[, overlaps_with_chipseq := any_overlaps]
}

add_chipseq_info_to_SNV_tf_binding_data <- function (SNV_tf_binding_data, chipseq) {
    data.table::setDT(SNV_tf_binding_data)
    data.table::setDT(chipseq)

    if(nrow(SNV_tf_binding_data) == 0) {
        # returns empty table if SNV_tf_binding_data are empty
        SNV_tf_binding_data[, overlaps_with_chipseq := character(0)]

        return (SNV_tf_binding_data)
    }

    if(nrow(chipseq) == 0) {
        # returns FALSE column if chipseq is empty
        SNV_tf_binding_data[, overlaps_with_chipseq := FALSE]

        return (SNV_tf_binding_data)
    }

    # convert data tables to GRanges
    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq,
        seqnames.field="chipseq_chrom",
        start.field="chipseq_start",
        end.field="chipseq_end",
        ignore.strand=TRUE
    )

    SNV_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        SNV_tf_binding_data,
        seqnames.field="chrom",
        start.field="pos",
        end.field="pos",
        ignore.strand=TRUE
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(chipseq_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SNV_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    n_overlaps <- GenomicRanges::countOverlaps(query = SNV_GRanges, subject = chipseq_GRanges)
    any_overlaps <- n_overlaps > 0

    SNV_tf_binding_data[, overlaps_with_chipseq := any_overlaps]
}