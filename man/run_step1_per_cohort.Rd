% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step1_per_cohort.R
\name{run_step1_per_cohort}
\alias{run_step1_per_cohort}
\title{Run Step 1 of the Revana workflow on a cohort/subgroup of tumor sample data}
\usage{
run_step1_per_cohort(
  paths_file_path,
  output_dir,
  gene_annotation_ref_file_path,
  gene_annotation_exons_ref_file_path,
  TAD_file_path,
  fimo_motif_ref_path = NULL,
  motif_id_tf_gene_name_table_path = NULL,
  genehancer_ref_file_path = NULL,
  chipseq_file_path = NULL,
  run_tf_binding_site_analysis = FALSE,
  reference_genome = NULL
)
}
\arguments{
\item{paths_file_path}{Path to the paths file. The paths file contains all relevant paths for all the samples to be included in the analysis. See the documentation for the exact format of this file.}

\item{output_dir}{Where are the analysis results to be stored}

\item{gene_annotation_ref_file_path}{Path to the reference file for the gene annotation. See the documentation for the exact format of this file and how to obtain it.}

\item{gene_annotation_exons_ref_file_path}{Path to the reference file for the exon annotation. See the documentation for the exact format of this file and how to obtain it.}

\item{TAD_file_path}{Path to the TAD file. The TAD file contains genomic coordinates of Topologically Associated Domains (TADs). See the documentation for the exact format of this file.}

\item{fimo_motif_ref_path}{Path to the reference file containing the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.}

\item{motif_id_tf_gene_name_table_path}{Path to the reference for the TF gene name annotation of the FIMO motifs. See the documentation for the exact format of this file and how to obtain it.}

\item{genehancer_ref_file_path}{Path to the reference for GeneHancer. See the documentation for the exact format of this file and how to obtain it.}

\item{chipseq_file_path}{Path to the ChIP-Seq file. The ChIP-Seq file contains genomic coordinates of regulatory active regions and can be obtained from experimental data or online less specific online resources. This file is optional, although an empty dummy file has to be used, if no ChIP-Seq data is supplied. See the documentation for the exact format of this file.}

\item{run_tf_binding_site_analysis}{Should TF binding site analysis be conducted?}

\item{reference_genome}{Which genome should be used. Default (NULL) uses Human Genome build GRCh37 (hg19)}
}
\description{
Run Step 1 of the Revana workflow on a cohort/subgroup of tumor sample data
}
\details{
Step 1 of the revana workflow ...
\itemize{
\item validates the paths file
\item validates all reference files
\item validates all sample files that are referenced in the paths file
\item creates the results folder structure in the output directory
}
}
