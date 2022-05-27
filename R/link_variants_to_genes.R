link_SVs_to_genes <- function(SVs, cis_activation_summary_TADs) {
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = SVs,
        dt2 = cis_activation_summary_TADs,
        dt1.seqnames.field1 = "sv_break1_chrom",
        dt1.start.field1 = "sv_break1_pos",
        dt1.end.field1 = "sv_break1_pos",
        dt1.seqnames.field2 = "sv_break2_chrom",
        dt1.start.field2 = "sv_break2_pos",
        dt1.end.field2 = "sv_break2_pos",
        dt2.seqnames.field  = "chrom",
        dt2.start.field  = "variant_overlap_window_start",
        dt2.end.field = "variant_overlap_window_end"
    ))
}

link_CNAs_to_genes <- function (CNAs, cis_activation_summary_TADs) {
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = CNAs,
        dt2 = cis_activation_summary_TADs,
        dt1.seqnames.field1 = "cna_chrom",
        dt1.start.field1 = "cna_start",
        dt1.end.field1 = "cna_start",
        dt1.seqnames.field2 = "cna_chrom",
        dt1.start.field2 = "cna_end",
        dt1.end.field2 = "cna_end",
        dt2.seqnames.field  = "chrom",
        dt2.start.field  = "variant_overlap_window_start",
        dt2.end.field = "variant_overlap_window_end"
    ))
}

link_SNVs_to_genes <- function (SNVs, cis_activation_summary_TADs) {
    return(merge_data_tables_by_overlap(
        dt1 = SNVs,
        dt2 = cis_activation_summary_TADs,
        seqnames.field1="snv_chrom",
        start.field1="snv_pos",
        end.field1="snv_pos",
        seqnames.field2 = "chrom",
        start.field2="variant_overlap_window_start",
        end.field2="variant_overlap_window_end"
    ))
}

link_SNV_tf_binding_data_to_genes <- function (SNV_tf_binding_data, cis_activation_summary_TADs) {
    data.table::setDT(SNV_tf_binding_data)

    if(nrow(SNV_tf_binding_data) == 0) {
        # returns empty table if SNV_tf_binding_data are empty
        data.table::setnames(
            SNV_tf_binding_data,
            c("chrom", "pos"),
            c("snv_chrom", "snv_pos")
        )

        return (cbind(SNV_tf_binding_data[FALSE, ], cis_activation_summary_TADs[FALSE, ]))
    }

    # convert data tables to GRanges
    genes_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        cis_activation_summary_TADs,
        seqnames.field="chrom",
        start.field="variant_overlap_window_start",
        end.field="variant_overlap_window_end",
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
    GenomeInfoDb::seqlevelsStyle(genes_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SNV_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    overlaps <- GenomicRanges::findOverlaps(query = SNV_GRanges, subject = genes_GRanges)


    # rename SNV cols for combination table
    data.table::setnames(
        SNV_tf_binding_data,
        c("chrom", "pos"),
        c("snv_chrom", "snv_pos")
    )


    SNV_tf_binding_data_gene_combinations <- cbind(
        SNV_tf_binding_data[S4Vectors::queryHits(overlaps)],
        cis_activation_summary_TADs[S4Vectors::subjectHits(overlaps)]
    )


    return(SNV_tf_binding_data_gene_combinations[order(gene_name)])
}

