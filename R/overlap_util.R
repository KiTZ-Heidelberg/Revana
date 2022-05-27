merge_data_tables_by_overlap <- function(
    dt1,
    dt2,
    seqnames.field1,
    start.field1,
    end.field1,
    seqnames.field2,
    start.field2,
    end.field2
    ) {
        data.table::setDT(dt1)
        data.table::setDT(dt2)

        if((nrow(dt1) == 0)|(nrow(dt2) == 0)) {
            # returns empty table if dt1 or dt2 are empty

            return (cbind(dt1[FALSE, ], dt2[FALSE, ]))
        }

        # convert data tables to GRanges
        dt1_GRange <- GenomicRanges::makeGRangesFromDataFrame(
            dt1,
            seqnames.field=seqnames.field1,
            start.field=start.field1,
            end.field=end.field1,
            ignore.strand=TRUE
        )

        dt2_GRange <- GenomicRanges::makeGRangesFromDataFrame(
            dt2,
            seqnames.field=seqnames.field2,
            start.field=start.field2,
            end.field=end.field2,
            ignore.strand=TRUE
        )

        # set common seqlevelStyle: UCSC => "chr1"
        GenomeInfoDb::seqlevelsStyle(dt1_GRange) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(dt2_GRange) <- "UCSC"

        # find overlaps between gene variant overlap window and CNAs
        overlaps <- GenomicRanges::findOverlaps(query = dt1_GRange, subject = dt2_GRange)

        #create merged table
        combinations <- cbind(
            dt1[S4Vectors::queryHits(overlaps)],
            dt2[S4Vectors::subjectHits(overlaps)]
        )

        return(combinations)
    }

merge_data_tables_by_overlap_dt1_has_two_Ranges <- function(
    dt1,
    dt2,
    dt1.seqnames.field1,
    dt1.start.field1,
    dt1.end.field1,
    dt1.seqnames.field2,
    dt1.start.field2,
    dt1.end.field2,
    dt2.seqnames.field,
    dt2.start.field,
    dt2.end.field
    ) {
        data.table::setDT(dt1)
        data.table::setDT(dt2)

    if((nrow(dt1) == 0)|(nrow(dt2) == 0)) {
        # returns empty table if dt1 or dt2 are empty

        return (cbind(dt1[FALSE, ], dt2[FALSE, ]))
    }
    
    # convert data tables to GRanges
    # convert data tables to GRanges
    dt1_GRange1 <- GenomicRanges::makeGRangesFromDataFrame(
        dt1,
        seqnames.field=dt1.seqnames.field1,
        start.field=dt1.start.field1,
        end.field=dt1.end.field1,
        ignore.strand=TRUE
    )

    dt1_GRange2 <- GenomicRanges::makeGRangesFromDataFrame(
        dt1,
        seqnames.field=dt1.seqnames.field2,
        start.field=dt1.start.field2,
        end.field=dt1.end.field2,
        ignore.strand=TRUE
    )

    dt2_GRange <- GenomicRanges::makeGRangesFromDataFrame(
        dt2,
        seqnames.field=dt2.seqnames.field,
        start.field=dt2.start.field,
        end.field=dt2.end.field,
        ignore.strand=TRUE
    )


    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(dt1_GRange1) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(dt1_GRange2) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(dt2_GRange) <- "UCSC"


    # find overlaps between gene varinant overlap window and CNAs
    overlaps1 <- GenomicRanges::findOverlaps(query = dt1_GRange1, subject = dt2_GRange)
    overlaps2 <- GenomicRanges::findOverlaps(query = dt1_GRange2, subject = dt2_GRange)

    overlaps1_df <- data.frame(query = S4Vectors::queryHits(overlaps1), subject = S4Vectors::subjectHits(overlaps1))
    overlaps2_df <- data.frame(query = S4Vectors::queryHits(overlaps2), subject = S4Vectors::subjectHits(overlaps2))

    # combine overlaps:
    # overlap between gene and one of the breakpoints is sufficient
    all_overlaps <- rbind(overlaps1_df, overlaps2_df)
    all_overlaps_unique <- all_overlaps[!duplicated(all_overlaps), ]

    combinations <- cbind(
        dt1[all_overlaps_unique$query],
        dt2[all_overlaps_unique$subject]
    )
    return(combinations)
}