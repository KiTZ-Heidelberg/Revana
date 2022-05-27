add_RNA_read_count_to_markers <- function(RNA_bam_file, marker_file_without_RNA_read_count, new_marker_file_path) {
    # validate columns
    col_types <- "ciccii"
    markers <- readr::read_tsv(marker_file_without_RNA_read_count, col_types = col_types)


    # process by chromosome

    chromosomes <- markers %>%
        dplyr::group_by(chrom) %>%
        dplyr::summarize(min_pos = min(pos), max_pos = max(pos))


    markers_with_pileup_list <- list()


    for (i in seq_len(nrow(chromosomes))) {
        chrom_i <- chromosomes$chrom[i]
        start_i <- chromosomes$min_pos[i]
        end_i <- chromosomes$max_pos[i]

        markers_of_chrom_i <- markers %>% dplyr::filter(chrom == chrom_i)

        data.table::setDT(markers_of_chrom_i)

        grange_i <- GenomicRanges::GRanges(seqnames = chrom_i, ranges = IRanges::IRanges(start = start_i, end = end_i))

        pileup_output <- Rsamtools::pileup(
            RNA_bam_file,
            scanBamParam = Rsamtools::ScanBamParam(
                which = grange_i
            ),
            pileupParam = Rsamtools::PileupParam(
                max_depth = 10000,
                min_mapq = 30,
                distinguish_strands = FALSE,
                include_deletions = FALSE,
                include_insertions = FALSE
            )
        )
        data.table::setDT(pileup_output)
        data.table::setnames(pileup_output, c("chrom", "pos", "nucleotide", "count", "which_label"))
        markers_with_pileup_long <- data.table::merge.data.table(
            x = markers_of_chrom_i,
            y = pileup_output,
            by = c("chrom", "pos"),
            all.x = TRUE
        )

        # markers_with_pileup_long <- dplyr::left_join(markers_of_chrom_i, pileup_output, by = c("chrom" = "seqnames", "pos" = "pos"))
        markers_with_pileup <- markers_with_pileup_long %>%
            dplyr::group_by(chrom, pos, ref, alt, reads.WGS.ref, reads.WGS.alt) %>%
            dplyr::summarise(
                reads.RNA.ref = sum(count[nucleotide == ref]),
                reads.RNA.alt = sum(count[nucleotide == alt])
            ) %>%
            dplyr::ungroup() %>%
            tidyr::replace_na(list(reads.RNA.ref = 0, reads.RNA.alt = 0))

        markers_with_pileup_list[[i]] <- markers_with_pileup
    }

    markers_with_pileup_all_chroms <- dplyr::bind_rows(markers_with_pileup_list)

    readr::write_tsv(markers_with_pileup_all_chroms, file = new_marker_file_path)
}