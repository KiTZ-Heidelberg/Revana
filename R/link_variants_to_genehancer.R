merge_genehancer_with_cis_activation_summary <- function(genehancer, cis_activation_summary) {
    genehancer_with_cis_activation_summary <- genehancer %>%
        dplyr::inner_join(cis_activation_summary, by = c("connected_gene" = "gene_name"), suffix = c("_genehancer", "_gene"))

    return(genehancer_with_cis_activation_summary)
}

link_SNVs_to_genehancer <- function(SNVs, genehancer_with_cis_activation_summary) {
    return(merge_data_tables_by_overlap(
        dt1 = SNVs,
        dt2 = genehancer_with_cis_activation_summary,
        seqnames.field1 = "snv_chrom",
        start.field1 = "snv_pos",
        end.field1 = "snv_pos",
        seqnames.field2 = "chrom_genehancer",
        start.field2 = "start_genehancer",
        end.field2 = "end_genehancer"
    ))
}

link_SNV_tf_binding_data_to_genehancer <- function(SNV_tf_binding_data, genehancer_with_cis_activation_summary) {
    data.table::setDT(SNV_tf_binding_data)
    data.table::setDT(genehancer_with_cis_activation_summary)

    if ((nrow(SNV_tf_binding_data) == 0) | (nrow(genehancer_with_cis_activation_summary) == 0)) {
        # returns empty table if SNV_tf_binding_data are empty
        data.table::setnames(
            SNV_tf_binding_data,
            c("chrom", "pos"),
            c("snv_chrom", "snv_pos")
        )

        return(cbind(SNV_tf_binding_data[FALSE, ], genehancer_with_cis_activation_summary[FALSE, ]))
    }

    # convert data tables to GRanges
    genehancer_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        genehancer_with_cis_activation_summary,
        seqnames.field = "chrom_genehancer",
        start.field = "start_genehancer",
        end.field = "end_genehancer",
        ignore.strand = TRUE
    )

    SNV_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        SNV_tf_binding_data,
        seqnames.field = "chrom",
        start.field = "pos",
        end.field = "pos",
        ignore.strand = TRUE
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(genehancer_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(SNV_GRanges) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    overlaps <- GenomicRanges::findOverlaps(query = SNV_GRanges, subject = genehancer_GRanges)

    # rename SNV cols for combination table
    data.table::setnames(
        SNV_tf_binding_data,
        c("chrom", "pos"),
        c("snv_chrom", "snv_pos")
    )


    SNV_tf_binding_data_genehancer_combinations <- cbind(
        SNV_tf_binding_data[S4Vectors::queryHits(overlaps)],
        genehancer_with_cis_activation_summary[S4Vectors::subjectHits(overlaps)]
    )


    return(SNV_tf_binding_data_genehancer_combinations)
}

link_CNAs_to_genehancer <- function(CNAs, genehancer_with_cis_activation_summary_TADs) {
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = CNAs,
        dt2 = genehancer_with_cis_activation_summary_TADs,
        dt1.seqnames.field1 = "cna_chrom",
        dt1.start.field1 = "cna_start",
        dt1.end.field1 = "cna_start",
        dt1.seqnames.field2 = "cna_chrom",
        dt1.start.field2 = "cna_end",
        dt1.end.field2 = "cna_end",
        dt2.seqnames.field = "chrom_genehancer",
        dt2.start.field = "variant_overlap_window_start_genehancer",
        dt2.end.field = "variant_overlap_window_end_genehancer"
    ))
}

link_SVs_to_genehancer <- function(SVs, genehancer_with_cis_activation_summary_TADs) {
    return(merge_data_tables_by_overlap_dt1_has_two_Ranges(
        dt1 = SVs,
        dt2 = genehancer_with_cis_activation_summary_TADs,
        dt1.seqnames.field1 = "sv_break1_chrom",
        dt1.start.field1 = "sv_break1_pos",
        dt1.end.field1 = "sv_break1_pos",
        dt1.seqnames.field2 = "sv_break2_chrom",
        dt1.start.field2 = "sv_break2_pos",
        dt1.end.field2 = "sv_break2_pos",
        dt2.seqnames.field = "chrom_genehancer",
        dt2.start.field = "variant_overlap_window_start_genehancer",
        dt2.end.field = "variant_overlap_window_end_genehancer"
    ))
}