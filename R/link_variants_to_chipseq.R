link_CNA_gene_combinations_to_chipseq <- function(CNA_gene_combinations, chipseq_with_TADs){
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = CNA_gene_combinations,
        dt2 = chipseq_with_TADs,
        dt1.seqnames.field1 = "cna_chrom",
        dt1.start.field1 = "cna_start",
        dt1.end.field1 = "cna_start",
        dt1.seqnames.field2 = "cna_chrom",
        dt1.start.field2 = "cna_end",
        dt1.end.field2 = "cna_end",
        dt2.seqnames.field  = "chipseq_chrom",
        dt2.start.field  = "chipseq_variant_overlap_window_start",
        dt2.end.field = "chipseq_variant_overlap_window_end"
    ))

}

link_SV_gene_combinations_to_chipseq <- function(CNA_gene_combinations, chipseq_with_TADs){
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = CNA_gene_combinations,
        dt2 = chipseq_with_TADs,
        dt1.seqnames.field1 = "sv_break1_chrom",
        dt1.start.field1 = "sv_break1_pos",
        dt1.end.field1 = "sv_break1_pos",
        dt1.seqnames.field2 = "sv_break2_chrom",
        dt1.start.field2 = "sv_break2_pos",
        dt1.end.field2 = "sv_break2_pos",
        dt2.seqnames.field  = "chipseq_chrom",
        dt2.start.field  = "chipseq_variant_overlap_window_start",
        dt2.end.field = "chipseq_variant_overlap_window_end"
    ))

}
link_CNAs_to_chipseq <- function (CNAs, chipseq_with_cis_activation_summary_TADs) {
    data.table::setDT(CNAs)

    if((nrow(CNAs) == 0) | (nrow(chipseq_with_cis_activation_summary_TADs) == 0)) {
        # returns empty table if CNAs are empty
        data.table::setnames(
            CNAs,
            c("chrom", "start", "end"),
        c("cna_chrom", "cna_start", "cna_end")
        )

        return (cbind(CNAs[FALSE, ], chipseq_with_cis_activation_summary_TADs[FALSE, ]))
    }
    
    # convert data tables to GRanges
    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq_with_cis_activation_summary_TADs,
        seqnames.field="chipseq_chrom",
        start.field="variant_overlap_window_start_chipseq",
        end.field="variant_overlap_window_end_chipseq",
        ignore.strand=TRUE
    )

    CNA_start_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        CNAs,
        seqnames.field="chrom",
        start.field="start",
        end.field="start",
        ignore.strand=TRUE
    )

    CNA_end_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        CNAs,
        seqnames.field="chrom",
        start.field="end",
        end.field="end",
        ignore.strand=TRUE
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(chipseq_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(CNA_start_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(CNA_end_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    overlaps1 <- GenomicRanges::findOverlaps(query = CNA_start_GRanges, subject = chipseq_GRanges)
    overlaps2 <- GenomicRanges::findOverlaps(query = CNA_end_GRanges, subject = chipseq_GRanges)

    overlaps1_df <- data.frame(query = S4Vectors::queryHits(overlaps1), subject = S4Vectors::subjectHits(overlaps1))
    overlaps2_df <- data.frame(query = S4Vectors::queryHits(overlaps2), subject = S4Vectors::subjectHits(overlaps2))

    # combine overlaps:
    # overlap between gene and one of the breakpoints is sufficient
    all_overlaps <- rbind(overlaps1_df, overlaps2_df)
    all_overlaps_unique <- all_overlaps[!duplicated(all_overlaps), ]

    # rename CNA cols for combination table
    data.table::setnames(
        CNAs,
        c("chrom", "start", "end"),
        c("cna_chrom", "cna_start", "cna_end")
    )

    CNA_chipseq_combinations <- cbind(
        CNAs[all_overlaps_unique$query],
        chipseq_with_cis_activation_summary_TADs[all_overlaps_unique$subject]
    )
    return(CNA_chipseq_combinations[order(gene_name)])
}

link_SVs_to_chipseq <- function (SVs, chipseq_with_cis_activation_summary_TADs) {
    data.table::setDT(SVs)
    data.table::setDT(chipseq_with_cis_activation_summary_TADs)

    if((nrow(SVs) == 0) | (nrow(chipseq_with_cis_activation_summary_TADs) == 0)) {
        # returns empty table if SVs are empty
        data.table::setnames(
            SVs,
            c("chrom1", "pos1", "chrom2", "pos2"),
            c("sv_break1_chrom", "sv_break1_pos", "sv_break2_chrom", "sv_break2_pos")
        )

        return (cbind(SVs[FALSE, ], chipseq_with_cis_activation_summary_TADs[FALSE, ]))
    }
    
    # convert data tables to GRanges
    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq_with_cis_activation_summary_TADs,
        seqnames.field="chipseq_chrom",
        start.field="variant_overlap_window_start_chipseq",
        end.field="variant_overlap_window_end_chipseq",
        ignore.strand=TRUE
    )

    SV_break1_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        SVs,
        seqnames.field="chrom1",
        start.field="pos1",
        end.field="pos1",
        ignore.strand=TRUE
    )

    SV_break2_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        SVs,
        seqnames.field="chrom2",
        start.field="pos2",
        end.field="pos2",
        ignore.strand=TRUE
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(chipseq_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SV_break1_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SV_break2_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and SVs
    overlaps1 <- GenomicRanges::findOverlaps(query = SV_break1_GRanges, subject = chipseq_GRanges)
    overlaps2 <- GenomicRanges::findOverlaps(query = SV_break2_GRanges, subject = chipseq_GRanges)

    overlaps1_df <- data.frame(query = S4Vectors::queryHits(overlaps1), subject = S4Vectors::subjectHits(overlaps1))
    overlaps2_df <- data.frame(query = S4Vectors::queryHits(overlaps2), subject = S4Vectors::subjectHits(overlaps2))

    # combine overlaps:
    # overlap between gene and one of the breakpoints is sufficient
    all_overlaps <- rbind(overlaps1_df, overlaps2_df)
    all_overlaps_unique <- all_overlaps[!duplicated(all_overlaps), ]

    # rename SV cols for combination table
    data.table::setnames(
        SVs,
        c("chrom1", "pos1", "chrom2", "pos2"),
        c("sv_break1_chrom", "sv_break1_pos", "sv_break2_chrom", "sv_break2_pos")
    )

    SV_chipseq_combinations <- cbind(
        SVs[all_overlaps_unique$query],
        chipseq_with_cis_activation_summary_TADs[all_overlaps_unique$subject]
    )
    return(SV_chipseq_combinations[order(gene_name)])
}

link_SNVs_to_chipseq <- function (SNVs, chipseq_with_cis_activation_summary) {
    return(merge_data_tables_by_overlap(
        dt1 = SNVs,
        dt2 = chipseq_with_cis_activation_summary,
        seqnames.field1="snv_chrom",
        start.field1="snv_pos",
        end.field1="snv_pos",
        seqnames.field2 = "chipseq_chrom",
        start.field2="chipseq_start",
        end.field2="chipseq_end"
    ))
}