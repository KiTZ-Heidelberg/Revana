% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_report.R
\name{create_results_HTML_report}
\alias{create_results_HTML_report}
\title{Create Interactive HTML report for Revana results}
\usage{
create_results_HTML_report(
  HTML_report_output_dir_path,
  output_paths_file_path,
  use_cache = F,
  create_cache = F,
  max_n_genes_to_analyze_further = 500,
  additional_genes_to_analyze_further = character(0),
  skip_creating_toggle_images = F,
  has_run_tf_binding_site_analysis = F,
  add_custom_HTML_to_gene_subdocument = NULL,
  add_custom_HTML_to_gene_subdocument_sample_data = NULL,
  genome_name = "hg19"
)
}
\arguments{
\item{HTML_report_output_dir_path}{path to the directory where the HTML report and its figures should be stored}

\item{output_paths_file_path}{path to the results path file - created by revana::run() . If several subgroups should be summarized the respective results paths files need to be merged in advance}

\item{use_cache}{debugging option to use cache for loaded data}

\item{create_cache}{debugging option to create cache for loaded data, if not already in place}

\item{max_n_genes_to_analyze_further}{max number of genes that are analyzed in 'by gene' analysis}

\item{additional_genes_to_analyze_further}{additional genes that are to be included in 'by gene' analysis}

\item{skip_creating_toggle_images}{should the cis-activated samples per gene - alias toggle plots - be created? defaults to TRUE}

\item{has_run_tf_binding_site_analysis}{Was TF binding site analysis conducted?}

\item{add_custom_HTML_to_gene_subdocument}{plugin function for custom HTML in by gene subdocument ()}

\item{add_custom_HTML_to_gene_subdocument_sample_data}{plugin function for custom HTML in by gene subdocument - sample data section}

\item{genome_name}{name of the reference genome - tested with "hg19" and "hg38"}
}
\description{
Create Interactive HTML report for Revana results
}
