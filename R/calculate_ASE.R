merge_markers_with_geneanno <- function (markers, geneanno) {
  markers_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
    markers,
    seqnames.field = "chrom",
    start.field = "POS",
    end.field = "POS",
    keep.extra.columns = F
  )
  geneanno_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
    geneanno,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = F
  )

  # set common seqlevelStyle: UCSC => "chr1"
  GenomeInfoDb::seqlevelsStyle(markers_GRanges) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(geneanno_GRanges) <- "UCSC"

  overlaps <- GenomicRanges::findOverlaps(markers_GRanges, geneanno_GRanges)

  geneanno_gene_names <- geneanno %>% dplyr::select(gene_name)
  markers_by_gene <- dplyr::bind_cols(
    dplyr::slice(markers, S4Vectors::queryHits(overlaps)),
    dplyr::slice(geneanno_gene_names, S4Vectors::subjectHits(overlaps))
  )
  return(markers_by_gene)
}

summarize_markers_by_gene <- function (markers_by_gene, geneanno) {
  marker_summary_by_gene <- markers_by_gene %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      n_markers = dplyr::n(),
      n_ASE_markers = sum(is_ASE_marker == T),
      n_cnv_markers = sum(copynumber_tag == "cnv"),
      p_all = paste(p_value, collapse = ","),
      delta_abs_all = paste(delta_abs, collapse = ","),
      copynumber_tag_all = paste(copynumber_tag, collapse = ","),
      combined_p_value = exp(sum(log(p_value))/dplyr::n()),
      mean_delta_abs = mean(delta_abs)
    )

  marker_summary_by_gene_extended <- dplyr::left_join(
    geneanno,
    marker_summary_by_gene, by = c("gene_name")) %>%
    tidyr::replace_na(replace = list(n_markers = 0, n_ASE_markers = 0, n_cnv_markers = 0, copynumber_tag_all = "", delta_abs_all = "", p_all = ""))

  return(marker_summary_by_gene_extended)
}

calculate_FDR_via_ABH <- function (geneanno_marker_summary) {
  if (nrow(geneanno_marker_summary) == 1) {
    returned_gene_summary <- geneanno_marker_summary %>% dplyr::mutate(Bonferroni = combined_p_value, ABH = combined_p_value)
  }else {
    raw.p <- geneanno_marker_summary$combined_p_value
    multtest.out <- multtest::mt.rawp2adjp(raw.p,c("Bonferroni","ABH"), na.rm = T)

    adj.p.matrix <- multtest.out$adj
    indexes <- order(multtest.out$index)

    adj.p.matrix <- adj.p.matrix[indexes,]
    ABH <- adj.p.matrix[,"ABH"]
    Bonferroni <- adj.p.matrix[,"Bonferroni"]
    returned_gene_summary <- geneanno_marker_summary %>% dplyr::mutate(Bonferroni = Bonferroni, ABH = ABH)
  }
  return(returned_gene_summary)
}

classify_biallelic_expression_OHE_matrix_filter <- function (geneanno_marker_summary) {
  returned_geneanno_marker_summary <- geneanno_marker_summary %>%
    dplyr::mutate(is_biallelic = ((combined_p_value >= 0.05) & (n_ASE_markers == 0))) %>%
    # for combined_p_value = NA due to lack of markers => assume biallelic expression
    tidyr::replace_na(list(is_biallelic = T))
  return(returned_geneanno_marker_summary)
}

check_genes_for_ASE_marker_run <- function(geneanno_marker_summary, runs){
  if(nrow(runs)>0){
    geneanno_marker_summary_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
      geneanno_marker_summary,
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = F
    )

    runs_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
      runs,
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = F
    )

    # set common seqlevelStyle: UCSC => "chr1"
    GenomeInfoDb::seqlevelsStyle(geneanno_marker_summary_GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(runs_GRanges) <- "UCSC"

    # cis-X approach
    count_overlap_results <- GenomicRanges::countOverlaps(geneanno_marker_summary_GRanges, runs_GRanges, type ="any")

    # more conservative approach
    # count_overlap_results <- GenomicRanges::countOverlaps(geneanno_marker_summary_GRanges, runs_GRanges, type ="within")
    does_overlap <- count_overlap_results >= 1
  
    returned_geneanno_marker_summary <- geneanno_marker_summary %>%
      dplyr::mutate(ASE_marker_run_overlap = does_overlap)
  }else{
    returned_geneanno_marker_summary <- geneanno_marker_summary %>%
      dplyr::mutate(ASE_marker_run_overlap = FALSE)
  }

  return(returned_geneanno_marker_summary)
}

calculate_ASE <- function (processed_markers, geneanno, runs) {
  markers_by_gene <- processed_markers %>%
    filter_heterozygous_markers(min_ref_or_mut_coverage_RNA = 10) %>%
    merge_markers_with_geneanno(geneanno)

  geneanno_marker_summary <- summarize_markers_by_gene(markers_by_gene, geneanno) %>%
    calculate_FDR_via_ABH() %>%
    classify_biallelic_expression_OHE_matrix_filter() %>%
    check_genes_for_ASE_marker_run(runs)

  return(geneanno_marker_summary)
}
