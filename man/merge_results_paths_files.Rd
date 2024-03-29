% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_results_paths_files.R
\name{merge_results_paths_files}
\alias{merge_results_paths_files}
\title{Combine results paths file from several analysed subgroups into one}
\usage{
merge_results_paths_files(
  list_of_results_paths_files,
  output_path_merged_results_file
)
}
\arguments{
\item{list_of_results_paths_files}{A list of the results paths files for each analysed subcohort. Each file is created as output paths summary by revana::run()}

\item{output_path_merged_results_file}{Path, where the new combined results paths file is to be stored}
}
\description{
Combine results paths file from several analysed subgroups into one
}
