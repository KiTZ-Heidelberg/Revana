
#' Asserting file existence
#'
#' @description throws an error if a necessary file doesn't exist under the provided path
#'
#' @param path provided path for the file that gets checked for existence
#' @param name_of_file_type name of the file type
#'
#' @return
#'
#' @examples
#' check_file_existence("/invalidpath/invalidfile.txt", name_of_file_type = "important text file")
#' check_file_existence("/invalidpath/invalidfile.txt")
#' @noRd
check_file_existence <- function(path, name_of_file_type = "") {
  if (name_of_file_type == "") {
    error_message_part_1 <- paste0("A required file does not exist.")
  } else {
    error_message_part_1 <- paste0("The provided file for ", name_of_file_type, " does not exist and is required.")
  }
  if (path == "") {
    stop(paste0("\n\t", error_message_part_1, "\n\tYou provided an empty string\n"))
  }
  if (!file.exists(path)) {
    stop(paste0("\n\t", error_message_part_1, "\n\tYou provided the following path: \n\t", path, "\n"))
  }
}

check_file_existence_auto_created <- function(path, name_of_file_type = "") {
  if (name_of_file_type == "") {
    error_message_part_1 <- paste0("A required output file of one of the preceding steps does not exist.")
  } else {
    error_message_part_1 <- paste0("The file ", name_of_file_type, "a required output of  one of the preceding steps and does not exist.")
  }
  if (path == "") {
    stop(paste0("\n\t", error_message_part_1, "\n\tThe expected path is an empty string\n"))
  }
  if (!file.exists(path)) {
    stop(paste0("\n\t", error_message_part_1, "\n\tExpected path: \n\t", path, "\n"))
  }
}

check_file_table_header <- function(path, required_cols, name_of_file_type = "") {
  con <- file(path, open = "r")
  line1 <- readLines(con, 1)
  close.connection(con)
  colnames_from_header <- stringr::str_split(line1, pattern = "\\t")[[1]]
  n_cols <- length(colnames_from_header)
  if (name_of_file_type == "") {
    error_message <- paste0(
      "The provided header of an imported file did not match the required columns: \n",
      "File: \n",
      path,
      "\n",
      "Required columns: \n",
      paste("\t", required_cols, sep = "", collapse = "\n"),
      "\n",
      "Provided columns: \n",
      paste("\t", colnames_from_header, sep = "", collapse = "\n"),
      "\n"
    )
  } else {
    error_message <- paste0(
      "The provided header of the file for ", name_of_file_type, " did not match the required columns: \n",
      "File: \n",
      path,
      "\n",
      "Required columns: \n",
      paste("\t", required_cols, sep = "", collapse = "\n"),
      "\n",
      "Provided columns: \n",
      paste("\t", colnames_from_header, sep = "", collapse = "\n"),
      "\n"
    )
  }
  if (!identical(colnames_from_header, required_cols)) {
    stop(error_message)
  }
}

check_data_structure <- function(tibble, required_cols, name_of_file_type = "") {
  if (name_of_file_type == "") {
    error_message <- paste0(
      "The provided columns of an imported file did not match the required columns: \n",
      paste("\t", required_cols, sep = "", collapse = "\n"),
      "\n"
    )
  } else {
    error_message <- paste0(
      "The provided columns for ", name_of_file_type, " did not match the required columns: \n",
      paste("\t", required_cols, sep = "", collapse = "\n"),
      "\n"
    )
  }
  if (!identical(colnames(tibble), required_cols)) {
    stop(error_message)
  }
}

check_output_dir_existence <- function(output_dir) {
  if (output_dir == "") {
    stop(paste0("\n\tThe output directory does not exist.\n\tYou provided an empty string\n"))
  }
  if (!dir.exists(output_dir)) {
    stop(paste0("\n\tThe output directory does not exist.\n\tPlease create the directory or provide a path to an existing directory.\n\tYou provided the following path: \n\t", output_dir, "\n"))
  }
}

check_sample_dir_existence <- function(sample_dir, sample_id) {
  if (sample_dir == "") {
    stop(paste0("\n\tThe sample directory for sample ", sample_id, " does not exist.\n\tAn empty string was provided\n"))
  }
  if (!dir.exists(sample_dir)) {
    stop(paste0("\n\tThe sample directory for sample ", sample_id, " does not exist.\n"))
  }
}