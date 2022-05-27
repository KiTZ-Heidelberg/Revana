calculate_marker_allele_frequency <- function (markers) {
  markers_with_AF <- markers %>%
    dplyr::mutate(
      WGS_alt_AF = reads.WGS.alt/(reads.WGS.ref + reads.WGS.alt),
      RNA_alt_AF = reads.RNA.alt/(reads.RNA.ref + reads.RNA.alt)
      )
  return(markers_with_AF)
}

filter_heterozygous_markers <- function (
  markers,
  min_total_coverage_WGS = 10,
  min_ref_or_mut_coverage_RNA = 5,
  # min_ref_or_mut_coverage_RNA = 5, => for classic ASE analysis
  min_ref_count_WGS = 3,
  min_alt_count_WGS = 3,
  min_alt_AF_WGS = 0.3,
  max_alt_AF_WGS = 0.7
  ) {
  filtered_markers <- markers %>%
    dplyr::filter(
      (reads.WGS.ref + reads.WGS.alt) >= min_total_coverage_WGS,
      (reads.RNA.ref + reads.RNA.alt) >= min_ref_or_mut_coverage_RNA,
      reads.WGS.ref >= min_ref_count_WGS,
      reads.WGS.alt >= min_alt_count_WGS,
      WGS_alt_AF >= min_alt_AF_WGS,
      WGS_alt_AF <= max_alt_AF_WGS
    )
  return(filtered_markers)
}

add_copy_number_tags_to_markers <- function (markers, CNAs) {
  markers_GRanges <- GenomicRanges::GRanges(seqnames = markers$chrom, ranges = IRanges::IRanges(start = markers$pos, width = 1))
  CNA_GRanges <- GenomicRanges::GRanges(seqnames = CNAs$chrom, ranges = IRanges::IRanges(start = CNAs$start, end = CNAs$end))

  # set common seqlevelStyle: UCSC => "chr1"
  GenomeInfoDb::seqlevelsStyle(markers_GRanges) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(CNA_GRanges) <- "UCSC"

  is_marker_copy_number_alterated <- IRanges::overlapsAny(markers_GRanges, CNA_GRanges)
  markers_tagged <- markers %>% dplyr::mutate(copynumber_tag = dplyr::if_else(is_marker_copy_number_alterated, "cnv", "diploid"))
  return(markers_tagged)
}



add_allelic_imbalance_measures_to_markers <- function (markers) {
  # still super slow: maybe use precalculated lookup table in future!
  calculate_pvalue <-
    function (reads.RNA.alt,
              reads.RNA.ref,
              balanced.expression) {
        covg <- reads.RNA.alt + reads.RNA.ref
        sigma <- 10.8 * (1 - exp(-1 * covg / 105))

        dist_binom <- dbinom(seq(0, covg), covg, balanced.expression)
        dist_norm  <- dnorm(seq(-1000, 1000), mean = 0, sd = sigma)
        dist_conv <- convolve(dist_binom, dist_norm, type = "open")

        if (reads.RNA.alt > covg * balanced.expression) {
          p_value <- sum(dist_conv[(1001 + reads.RNA.alt):length(dist_conv)])
        } else {
          p_value <- sum(dist_conv[1:(1001 + reads.RNA.alt)])
        }
        if (p_value < 0) {
          return_p_value <- 0
        } else if (p_value > 1) {
          return_p_value <- 1
        } else{
          return_p_value <- p_value
        }
        return(return_p_value)
    }

  balanced.expression <- dplyr::if_else(markers$copynumber_tag == "cnv", markers$WGS_alt_AF, 0.5)
  p_values <- list(markers$reads.RNA.alt, markers$reads.RNA.ref, balanced.expression) %>% purrr::pmap(~calculate_pvalue(..1,..2,..3))

  returned_markers <- markers %>%
    dplyr::mutate(
      delta_abs = dplyr::if_else(copynumber_tag == "cnv", abs(RNA_alt_AF-WGS_alt_AF), abs(RNA_alt_AF-0.5)),
      p_value = unlist(p_values)
      )

  return(returned_markers)
}

classify_ASE_of_single_markers <- function (
  markers,
  delta_abs_threshold_single_marker_diploid = 0.3,
  delta_abs_threshold_single_marker_cnv = 0.2,
  p_value_threshold = 0.05
  ) {
  returned.markers <-markers %>%
    dplyr::mutate(
      is_ASE_marker = (p_value <= p_value_threshold) & (dplyr::if_else(
        copynumber_tag == "cnvloh",
        delta_abs >= delta_abs_threshold_single_marker_cnv,
        delta_abs >= delta_abs_threshold_single_marker_diploid
        )))
  return(returned.markers)
}

#' Processessing of Single Nucleotide Markers
#'
#' @param markers tibble of markers
#' @param CNAs tibble of CNAs
#' @importFrom magrittr %>%
process_markers <- function(
  markers,
  CNAs
) {
    processed_markers <- markers %>%
        calculate_marker_allele_frequency() %>%
        filter_heterozygous_markers() %>%
        add_copy_number_tags_to_markers(CNAs) %>%
        add_allelic_imbalance_measures_to_markers() %>%
        classify_ASE_of_single_markers ()
    return(processed_markers)
}