add_TAD_boundaries_to_gene_annotation <- function(gene_annotation, TADs) {
    data.table::setDT(TADs)
    data.table::setDT(gene_annotation)

    genes_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        gene_annotation,
        seqnames.field="chrom",
        start.field="start",
        end.field="end",
        ignore.strand=TRUE
    )
    TAD_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        TADs,
        seqnames.field="chrom",
        start.field="start",
        end.field="end",
        ignore.strand=TRUE
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
        )
        , by = "gene_name"
    ]

    gene_annotation_TADs <- merge(x = gene_annotation, y = TAD_info_per_gene, by = "gene_name", all.x = TRUE) # left join

    gene_annotation_TADs[is.na(n_overlapping_TADs), n_overlapping_TADs:= 0]

    return(gene_annotation_TADs)
}

add_variant_overlap_window_to_gene_annotation <- function (gene_annotation) {
    gene_annotation[ , variant_overlap_window_start := pmin(start, min_TAD_start, na.rm = TRUE)]
    gene_annotation[ , variant_overlap_window_end := pmax(end, max_TAD_end, na.rm = TRUE)]
    return(gene_annotation)
}

link_chipseq_to_genes <- function (chipseq, gene_annotation) {
    data.table::setDT(chipseq)
    data.table::setDT(gene_annotation)

    if(nrow(chipseq) == 0) {
        # returns empty table if chipseq are empty
        data.table::setnames(
            gene_annotation,
            old = c("chrom","start","end"),
            new = c("gene_annotation_chrom", "gene_annotation_start", "gene_annotation_end"))

        return (cbind(gene_annotation[FALSE, .(gene_name)], chipseq[FALSE, ]))
    }

    genes_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        gene_annotation,
        seqnames.field="chrom",
        start.field="variant_overlap_window_start",
        end.field="variant_overlap_window_end",
        ignore.strand=TRUE
    )
    chipseq_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
        chipseq,
        seqnames.field="chrom",
        start.field="start",
        end.field="end",
        ignore.strand=TRUE
    )

    data.table::setnames(gene_annotation, old = c("chrom","start","end"), new = c("gene_annotation_chrom", "gene_annotation_start", "gene_annotation_end"))
    overlaps <- GenomicRanges::findOverlaps(query = chipseq_GRanges, subject = genes_GRanges)

    chipseq_gene_combinations <- dplyr::bind_cols(
        gene_annotation[S4Vectors::subjectHits(overlaps), .(gene_name)],
        chipseq[S4Vectors::queryHits(overlaps)]
    )

    return(chipseq_gene_combinations)
}

merge_chipseq_gene_combination_with_cis_activation_summary <- function (chipseq_gene_combinations, cis_activation_summary) {
    data.table::setDT(chipseq_gene_combinations)
    data.table::setDT(cis_activation_summary)

    data.table::setnames(chipseq_gene_combinations, old = c("chrom","start","end"), new = c("chipseq_chrom", "chipseq_start", "chipseq_end"))

    chipseq_with_cis_activation_summary <- merge(x = chipseq_gene_combinations, y = cis_activation_summary, by = "gene_name")
    
    return(chipseq_with_cis_activation_summary)
}