add_TAD_boundaries_to_cis_activation_summary <- function(cis_activation_summary, TADs) {
    # maybe implement gene padding to take into account distance to promoter region
    # also maybe make padding dependent on strand orientation
    data.table::setDT(TADs)

    genes_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        cis_activation_summary,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )
    TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        TADs,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )

    data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))
    overlaps <- GenomicRanges::findOverlaps(query = TAD_GRanges, subject = genes_GRanges)
    # union_data = all TAD gene combinations
    union_data <- dplyr::bind_cols(
        cis_activation_summary[S4Vectors::subjectHits(overlaps), .(gene_name)],
        TADs[S4Vectors::queryHits(overlaps)]
    )

    TAD_info_per_gene <- union_data[
        , .(
            n_overlapping_TADs = .N,
            min_TAD_start = min(TAD_start),
            max_TAD_end = max(TAD_end)
        ),
        by = "gene_name"
    ]

    cis_activation_summary_TADs <- merge(x = cis_activation_summary, y = TAD_info_per_gene, by = "gene_name", all.x = TRUE) # left join
    # cis_activation_summary %>%
    #     dplyr::left_join(TAD_info_per_gene, by = c("gene_name"))

    cis_activation_summary_TADs[is.na(n_overlapping_TADs), n_overlapping_TADs := 0]

    return(cis_activation_summary_TADs)
}

add_variant_overlap_window <- function(cis_activation_summary_TADs) {
    cis_activation_summary_TADs[, variant_overlap_window_start := pmin(start, min_TAD_start, na.rm = TRUE)]
    cis_activation_summary_TADs[, variant_overlap_window_end := pmax(end, max_TAD_end, na.rm = TRUE)]
}

add_TAD_boundaries_to_gene_annotation <- function(gene_annotation, TADs) {
    data.table::setDT(TADs)
    data.table::setDT(gene_annotation)

    genes_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        gene_annotation,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )
    TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        TADs,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )

    data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))
    overlaps <- GenomicRanges::findOverlaps(query = TAD_GRanges, subject = genes_GRanges)
    # union_data = all TAD gene combinations
    union_data <- dplyr::bind_cols(
        gene_annotation[S4Vectors::subjectHits(overlaps), .(gene_name)],
        TADs[S4Vectors::queryHits(overlaps)]
    )


    data.table::setDT(union_data)
    print(class(union_data))


    TAD_info_per_gene <- union_data[
        , .(
            n_overlapping_TADs = .N,
            min_TAD_start = min(TAD_start),
            max_TAD_end = max(TAD_end)
        ),
        by = "gene_name"
    ]

    gene_annotation_TADs <- merge(x = gene_annotation, y = TAD_info_per_gene, by = "gene_name", all.x = TRUE) # left join

    gene_annotation_TADs[is.na(n_overlapping_TADs), n_overlapping_TADs := 0]

    return(gene_annotation_TADs)
}

add_variant_overlap_window_to_gene_annotation <- function(gene_annotation) {
    gene_annotation[, variant_overlap_window_start := pmin(start, min_TAD_start, na.rm = TRUE)]
    gene_annotation[, variant_overlap_window_end := pmax(end, max_TAD_end, na.rm = TRUE)]
}

add_TAD_boundaries_to_genehancer_with_cis_activation_summary <- function(genehancer_with_cis_activation_summary, TADs) {
    # maybe implement gene padding to take into account distance to promoter region
    # also maybe make padding dependent on strand orientation
    data.table::setDT(TADs)
    data.table::setDT(genehancer_with_cis_activation_summary)

    if ((nrow(TADs) == 0) | (nrow(genehancer_with_cis_activation_summary) == 0)) {
        # returns empty table if TADs or genehancers are empty
        data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))

        empty_TAD_info_per_gene <- data.table::setDT(data.frame(
            n_overlapping_TADs_genehancer = integer(),
            min_TAD_start_genehancer = integer(),
            max_TAD_end_genehancer = integer()
        ))

        return(cbind(genehancer_with_cis_activation_summary[FALSE, ], empty_TAD_info_per_gene[FALSE, ]))
    }

    genehancer_with_cis_activation_summary_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        genehancer_with_cis_activation_summary,
        seqnames.field = "chrom_genehancer",
        start.field = "start_genehancer",
        end.field = "end_genehancer",
        ignore.strand = TRUE
    )
    TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        TADs,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )

    data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))
    overlaps <- GenomicRanges::findOverlaps(query = TAD_GRanges, subject = genehancer_with_cis_activation_summary_GRanges)
    # union_data = all TAD gene combinations
    union_data <- dplyr::bind_cols(
        genehancer_with_cis_activation_summary[S4Vectors::subjectHits(overlaps), .(chrom_genehancer, start_genehancer, end_genehancer, connected_gene)],
        TADs[S4Vectors::queryHits(overlaps)]
    )

    TAD_info_per_gene <- union_data[
        , .(
            n_overlapping_TADs_genehancer = .N,
            min_TAD_start_genehancer = min(TAD_start),
            max_TAD_end_genehancer = max(TAD_end)
        ),
        by = c("chrom_genehancer", "start_genehancer", "end_genehancer", "connected_gene")
    ]

    genehancer_with_cis_activation_summary_TADs <- merge(x = genehancer_with_cis_activation_summary, y = TAD_info_per_gene, by = c("chrom_genehancer", "start_genehancer", "end_genehancer", "connected_gene"), all.x = TRUE) # left join
    # cis_activation_summary %>%
    #     dplyr::left_join(TAD_info_per_gene, by = c("gene_name"))

    genehancer_with_cis_activation_summary_TADs[is.na(n_overlapping_TADs_genehancer), n_overlapping_TADs_genehancer := 0]

    return(genehancer_with_cis_activation_summary_TADs)
}

add_variant_overlap_window_to_genehancer_with_cis_activation_summary <- function(genehancer_with_cis_activation_summary) {
    genehancer_with_cis_activation_summary[, variant_overlap_window_start_genehancer := pmin(start_genehancer, min_TAD_start_genehancer, na.rm = TRUE)]
    genehancer_with_cis_activation_summary[, variant_overlap_window_end_genehancer := pmax(end_genehancer, max_TAD_end_genehancer, na.rm = TRUE)]
}

add_TAD_boundaries_to_chipseq_with_cis_activation_summary <- function(chipseq_with_cis_activation_summary, TADs = NULL) {
    chipseq_with_cis_activation_summary_TADs <- chipseq_with_cis_activation_summary

    chipseq_with_cis_activation_summary_TADs[, n_overlapping_TADs_chipseq := n_overlapping_TADs]
    chipseq_with_cis_activation_summary_TADs[, min_TAD_start_chipseq := min_TAD_start]
    chipseq_with_cis_activation_summary_TADs[, max_TAD_end_chipseq := max_TAD_end]

    return(chipseq_with_cis_activation_summary_TADs)
}

add_variant_overlap_window_to_chipseq_with_cis_activation_summary <- function(chipseq_with_cis_activation_summary_TADs) {
    chipseq_with_cis_activation_summary_TADs[, variant_overlap_window_start_chipseq := variant_overlap_window_start]
    chipseq_with_cis_activation_summary_TADs[, variant_overlap_window_end_chipseq := variant_overlap_window_end]
}

add_TAD_boundaries_to_chipseq <- function(chipseq, TADs) {
    data.table::setDT(TADs)
    data.table::setDT(chipseq)

    if ((nrow(TADs) == 0) | (nrow(chipseq) == 0)) {
        # returns empty table if TADs or genehancers are empty
        data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))

        empty_TAD_info_per_chipseq <- data.table::setDT(data.frame(
            n_overlapping_TADs_chipseq = integer(),
            min_TAD_start_chipseq = integer(),
            max_TAD_end_chipseq = integer()
        ))

        return(cbind(chipseq[FALSE, ], empty_TAD_info_per_chipseq[FALSE, ]))
    }


    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq,
        seqnames.field = "chipseq_chrom",
        start.field = "chipseq_start",
        end.field = "chipseq_end",
        ignore.strand = TRUE
    )
    TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        TADs,
        seqnames.field = "chrom",
        start.field = "start",
        end.field = "end",
        ignore.strand = TRUE
    )

    data.table::setnames(TADs, c("TAD_chrom", "TAD_start", "TAD_end"))

    overlaps <- GenomicRanges::findOverlaps(query = TAD_GRanges, subject = chipseq_GRanges)
    # union_data = all TAD gene combinations
    union_data <- dplyr::bind_cols(
        chipseq[S4Vectors::subjectHits(overlaps), .(chipseq_chrom, chipseq_start, chipseq_end)],
        TADs[S4Vectors::queryHits(overlaps)]
    )

    data.table::setDT(union_data)


    TAD_info_per_chipseq <- union_data[
        , .(
            n_overlapping_TADs_chipseq = .N,
            min_TAD_start_chipseq = min(TAD_start),
            max_TAD_end_chipseq = max(TAD_end)
        ),
        by = c("chipseq_chrom", "chipseq_start", "chipseq_end")
    ]

    chipseq_TADs <- merge(x = chipseq, y = TAD_info_per_chipseq, by = c("chipseq_chrom", "chipseq_start", "chipseq_end"), all.x = TRUE) # left join

    chipseq_TADs[is.na(n_overlapping_TADs_chipseq), n_overlapping_TADs_chipseq := 0]

    return(chipseq_TADs)
}

add_variant_overlap_window_to_chipseq <- function(chipseq) {
    chipseq[, chipseq_variant_overlap_window_start := pmin(chipseq_start, min_TAD_start_chipseq, na.rm = TRUE)]
    chipseq[, chipseq_variant_overlap_window_end := pmax(chipseq_end, max_TAD_end_chipseq, na.rm = TRUE)]
}

add_TAD_combination_to_SVs <- function(SVs, TADs) {
    a <- merge_data_tables_by_overlap(
        dt1 = SVs,
        dt2 = TADs,
        seqnames.field1 = "sv_break1_chrom",
        start.field1 = "sv_break1_pos",
        end.field1 = "sv_break1_pos",
        seqnames.field2 = "chrom",
        start.field2 = "start",
        end.field2 = "end"
    )

    data.table::setnames(a, old = c("chrom", "start", "end"), new = c("TAD_sv_break1_chrom", "TAD_sv_break1_start", "TAD_sv_break1_end"))

    b <- merge_data_tables_by_overlap(
        dt1 = a,
        dt2 = TADs,
        seqnames.field1 = "sv_break2_chrom",
        start.field1 = "sv_break2_pos",
        end.field1 = "sv_break2_pos",
        seqnames.field2 = "chrom",
        start.field2 = "start",
        end.field2 = "end"
    )

    data.table::setnames(b, old = c("chrom", "start", "end"), new = c("TAD_sv_break2_chrom", "TAD_sv_break2_start", "TAD_sv_break2_end"))

    c <- b %>%
        dplyr::mutate(
            TAD_sv_break1_name = paste0(TAD_sv_break1_chrom, ":", TAD_sv_break1_start, "-", TAD_sv_break1_end),
            TAD_sv_break2_name = paste0(TAD_sv_break2_chrom, ":", TAD_sv_break2_start, "-", TAD_sv_break2_end)
        ) %>%
        dplyr::mutate(TAD_combination = paste0(pmin(TAD_sv_break1_name, TAD_sv_break2_name), " <-> ", pmax(TAD_sv_break1_name, TAD_sv_break2_name)))

    d <- SVs %>%
        dplyr::left_join(
            c,
            by = c(
                "sv_break1_chrom",
                "sv_break1_pos",
                "sv_break2_chrom",
                "sv_break2_pos",
                "sv_type",
                "eventInversion"
            )
        )

    return(d)
}