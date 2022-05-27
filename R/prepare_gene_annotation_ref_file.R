#' Prepare the gene annotation reference files from default online sources
#'
#' @param gencode_gtf_file_path Path to the Gencode GTF file. Please download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz
#' @param imprint_genes_file_path Path to the Gencode GTF file. Please download from https://www.geneimprint.com/site/genes-by-species or SOURCE
#' @param cancer_gene_file_path Path to the cancer gene census file. Please download from https://cancer.sanger.ac.uk/census or SOURCE
#' @param ref_file_output_path Output path for the gene annotion reference file as required for Revana
#' @param ref_file_output_path_exons Output path for the exon annotion reference file as required for Revana
#'
#' @export
#'
prepare_gene_annotation_ref_file <- function(gencode_gtf_file_path, imprint_genes_file_path, cancer_gene_file_path, ref_file_output_path, ref_file_output_path_exons) {
  # File import --------------------------
  gencode_gtf_file_con <- file(gencode_gtf_file_path, open = "r")
  gtfimport_df <- data.frame(rtracklayer::import(gencode_gtf_file_con))
  close(gencode_gtf_file_con)




  # GENES ##############################
  # Filtering --------------------------
  gtfimport_df_filtered <- gtfimport_df %>%
    # keep only genes
    dplyr::filter(type == "gene") %>%
    # keep only longest transcript
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(order_by = width, n = 1, with_ties = F) %>%
    dplyr::ungroup() %>%
    # select columns
    dplyr::select(chrom = seqnames, start, end, width, strand, gene_name, gene_type)

  # add imprinting information ---------
  if (!(imprint_genes_file_path == "")) {
    imprint_genes <- readr::read_tsv(imprint_genes_file_path)
    imprint_genes_formatted <- imprint_genes %>%
      dplyr::mutate(gene_name = dplyr::if_else(is.na(Aliases) | (Aliases == ""), Gene, paste0(Gene, ", ", Aliases))) %>%
      tidyr::separate_rows(gene_name, sep = ",\\s*") %>%
      dplyr::select(gene_name, imprinting_status = Status, imprinting_expressed_allele = `Expressed Allele`) %>%
      dplyr::distinct(gene_name, .keep_all = T)
  } else {
    imprint_genes_formatted <- data.frame(
      gene_name = character(0),
      imprinting_status = character(0),
      imprinting_expressed_allele = character(0)
    )
  }

  gtfimport_df_filtered_with_imprint <- gtfimport_df_filtered %>%
    dplyr::left_join(imprint_genes_formatted, by = c("gene_name")) %>%
    tidyr::replace_na(list(imprinting_status = "no_imprinting"))


  # add cancer_gene information ---------
  if (!(cancer_gene_file_path == "")) {
    cancer_genes <- readr::read_tsv(cancer_gene_file_path)

    cancer_genes_formatted <- cancer_genes %>%
      dplyr::filter(!(`Gene Symbol` %in% c("IGH", "IGK", "IGL", "HLA-A"))) %>%
      dplyr::mutate(is_cancer_gene = TRUE) %>%
      dplyr::select(gene_name = Synonyms, cancer_gene_role_in_cancer = `Role in Cancer`, is_cancer_gene) %>%
      tidyr::separate_rows(gene_name, sep = ",") %>%
      dplyr::distinct(gene_name, .keep_all = T)
  } else {
    cancer_genes_formatted <- data.frame(
      gene_name = character(0),
      cancer_gene_role_in_cancer = character(0),
      is_cancer_gene = logical(0)
    )
  }

  gtfimport_df_filtered_with_imprint_with_cancer_gene <- gtfimport_df_filtered_with_imprint %>%
    dplyr::left_join(cancer_genes_formatted, by = c("gene_name")) %>%
    tidyr::replace_na(list(is_cancer_gene = FALSE))

  # Write file to output -------------------
  readr::write_tsv(
    gtfimport_df_filtered_with_imprint_with_cancer_gene,
    file = ref_file_output_path,
    col_names = T
  )

  # EXONS ###############################
  gtfimport_df_filtered_exons <- gtfimport_df %>%
    # keep only genes
    dplyr::filter(type == "exon") %>%
    # select columns
    dplyr::select(chrom = seqnames, start, end, gene_name)


  # Write file to output -------------------
  readr::write_tsv(
    gtfimport_df_filtered_exons,
    file = ref_file_output_path_exons,
    col_names = T
  )
}