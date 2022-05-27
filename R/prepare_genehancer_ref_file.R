#' Prepare the GeneHancer reference file
#'
#' @param genehancer_data Data frame containing the GeneHancer data (can be imported from Excel or GFF files)
#' @param output_path Output path for the Genehancer reference file
#' @param genehancer_element_elite_status_file_path If available provide elite status data for GeneHancer elements. This argument takes in the file path, not the dataframe
#' @param genehancer_gene_associations_scores_file_path If available provide elite association data. This argument takes in the file path, not the dataframe
#' @param keep_only_double_elites Should Revana only consider double elite GeneHancer regions
#' @param skip_lifting_over By default data is lifted over from GRCh38 (hg38) to GRCh37 (hg19) coordinates. If this argument is set to true, original coordinates will be kept. Set this to TRUE if you plan on running Revana on GRCh38 data
#'
#'
#' @export
#'
prepare_genehancer_ref_file <- function(genehancer_data,
                                        output_path,
                                        genehancer_element_elite_status_file_path = NULL,
                                        genehancer_gene_associations_scores_file_path = NULL,
                                        keep_only_double_elites = FALSE,
                                        skip_lifting_over = FALSE) {
    if (is.null(genehancer_element_elite_status_file_path)) {
        genehancer_element_elite_status_data <- data.frame(
            genehancer_id = character(0),
            is_elite = logical(0)
        )
    } else {
        genehancer_element_elite_status_data <- readr::read_tsv(genehancer_element_elite_status_file_path) %>%
            dplyr::mutate(is_elite = as.logical(is_elite)) %>%
            dplyr::select(
                genehancer_id = GHid,
                is_elite = is_elite
            )
    }

    if (is.null(genehancer_gene_associations_scores_file_path)) {
        genehancer_gene_associations_scores_data <- data.frame(
            genehancer_id = character(0),
            connected_gene = character(0),
            is_association_elite = logical(0)
        )
    } else {
        genehancer_gene_associations_scores_data <- readr::read_tsv(genehancer_gene_associations_scores_file_path) %>%
            dplyr::mutate(is_elite = as.logical(is_elite)) %>%
            dplyr::select(
                genehancer_id = GHid,
                connected_gene = symbol,
                is_association_elite = is_elite
            )
    }

    if ("#chrom" %in% names(genehancer_data)) {
        genehancer_data <- genehancer_data %>% dplyr::rename(chrom = `#chrom`)
    }

    genehancer_data_long <- genehancer_data %>%
        dplyr::select(
            chrom,
            feature_name = `feature name`,
            start,
            end,
            score,
            attributes
        ) %>%
        dplyr::mutate(
            genehancer_id = stringr::str_remove(
                stringr::str_remove(
                    stringr::str_extract(attributes, pattern = "^genehancer_id=.+?;"),
                    pattern = "genehancer_id="
                ),
                pattern = ";$"
            )
        ) %>%
        dplyr::mutate(
            attributes = stringr::str_remove(attributes, pattern = "^genehancer_id=.+?;")
        ) %>%
        dplyr::mutate(
            attributes = stringr::str_remove(attributes, pattern = "^connected_gene=")
        ) %>%
        tidyr::separate_rows(attributes, sep = ";connected_gene=") %>%
        tidyr::separate(col = attributes, into = c("connected_gene", "connected_gene_score"), sep = ";score=") %>%
        # add genehancer_element_elite_status
        dplyr::left_join(genehancer_element_elite_status_data, by = c("genehancer_id")) %>%
        tidyr::replace_na(list("is_elite" = FALSE)) %>%
        dplyr::left_join(genehancer_gene_associations_scores_data, by = c("genehancer_id", "connected_gene")) %>%
        tidyr::replace_na(list("is_association_elite" = FALSE))

    if (keep_only_double_elites == TRUE) {
        genehancer_data_long <- genehancer_data_long %>%
            dplyr::filter((is_elite & is_association_elite) == TRUE)
    }

    if (skip_lifting_over == FALSE) {
        chain_path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
        chain <- rtracklayer::import.chain(chain_path)

        genehancer_data_long <- tibble::rowid_to_column(genehancer_data_long, "ID")

        genehancer_data_long_granges <- GenomicRanges::makeGRangesFromDataFrame(
            genehancer_data_long,
            keep.extra.columns = TRUE,
            ignore.strand = TRUE
        )

        GenomeInfoDb::seqlevelsStyle(genehancer_data_long_granges) <- "UCSC"

        genehancer_data_long_granges_hg19 <- unlist(rtracklayer::liftOver(genehancer_data_long_granges, chain))

        genehancer_data_long_grange_as_df_hg19 <- data.frame(
            chrom_new = GenomicRanges::seqnames(genehancer_data_long_granges_hg19),
            start_new = GenomicRanges::start(genehancer_data_long_granges_hg19),
            end_new = GenomicRanges::end(genehancer_data_long_granges_hg19),
            ID = genehancer_data_long_granges_hg19$ID
        )

        genehancer_data_long <- genehancer_data_long %>%
            dplyr::left_join(genehancer_data_long_grange_as_df_hg19, by = "ID") %>%
            tidyr::drop_na(chrom_new) %>%
            tidyr::drop_na(start_new) %>%
            tidyr::drop_na(end_new) %>%
            dplyr::mutate(chrom = chrom_new, start = start_new, end = end_new) %>%
            dplyr::select(!chrom_new) %>%
            dplyr::select(!start_new) %>%
            dplyr::select(!end_new) %>%
            dplyr::select(!ID)
    }

    readr::write_tsv(genehancer_data_long, file = output_path)
}