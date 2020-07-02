sequence_input.character <- function(input){
  # character input should be file path unless input == "shuffle"
  if (input == "shuffle") return(input)

  cmdr::cmd_files_exist(input)
  return(input)
}

sequence_input.DNAStringSet <- function(input){
  # stringset is written to temporary fasta file
  write_fasta(input)
}

sequence_input.BStringSet <- function(input){
  write_fasta(input)
}

sequence_input.AAStringSet <- function(input){
  write_fasta(input)
}

motif_input.character <- function(input, path = NULL){
  # character input must be file path,
  # subsequent type checking will be carried out by commandline utility
  # path != NULL, copy & rename to this name
  if (length(input) > 1){
    stop("input must be character of length == 1")
  }
  if (!is.null(path)){
    if (path != "" | !is.na(path)){
      file.copy(input, path)
      input <- path
    }
  }

  cmdr::cmd_files_exist(input)

  out <- list(metadata = NULL,
              path = input)
  return(out)
}

motif_input.data.frame <- function(input, path = tempfile(fileext = ".meme")){

  if (!("motif" %in% names(input))) {
    stop("input data.frame must contain \"motif\" column")
  }

  if (!is_universalmotif_list(input$motif)) {
    stop("motif column is not in universalmotif format")
  }

  path <- input$motif %>%
    write_meme_input_path(path = path)

  out <- list(metadata = input,
              path = path)
  return(out)
}

motif_input.list <- function(input, path = tempfile(fileext = ".meme")){
  # check list is universalmotif list
  if (!is_universalmotif_list(input)) error_universalmotif_list(list)

  df <- as_universalmotif_dataframe(input)

  path <- input %>%
    write_meme_input_path(path = path)

  out <- list(metadata = df,
              path = path)

  return(out)
}

motif_input.universalmotif <- function(input, path = tempfile(fileext = ".meme")){

  df <- as_universalmotif_dataframe(input) %>%
    data.frame

  path <- input %>%
    write_meme_input_path(path = path)

  out <- list(metadata = df,
              path = path)

  return(out)
}

#' Helper for exporting universalmotif data to temp files
#'
#' @param input universalmotif list
#' @param path if path = "" (or NULL or NA) use tempfile, default: tempfile
#'
#' @return
#'
#' @noRd
write_meme_input_path <- function(input, path){
  if (is.null(path)){
      path <- tempfile(fileext = ".meme")
  }
  if (path == "" | is.null(path) | is.na(path)) {
    path <- tempfile(fileext = ".meme")
  }

  input %>%
    write_meme_list(path = as.character(path))
}

#' filter out control entries from list input
#'
#' @param input list of Biostrings::XStringSet objects
#' @param control character vector matching names in input
#'
#' @return list with `input` and `control` values
#'
#' @noRd
split_input_control <- function(input, control){
  # input = list
  # control = character vector
  input_minus_control <- input %>%
    cmdr::cmd_list_drop_named(control)

  control_seq <- input %>%
    cmdr::cmd_list_keep_named(control) %>%
    Biostrings::BStringSetList(., use.names = FALSE) %>%
    unlist

  input <- input_minus_control
  control <- control_seq

  return(
    list(input = input,
         control = control)
        )

}

#' Correctly handle input/control input logic when input is a list
#'
#' This is the backend that allows operations like using a name in input as control
#' or when passing a list of sequences to control to use the pool of them. used
#' in runDreme and runAme, or any other command that takes sequences as input and control
#'
#' @param input sequences
#' @param control control
#'
#' @return list w/ correct input & controls
#'
#' @noRd
sequence_input_control_list <- function(input, control){

  if (is.character(control)){
    if (all(control %in% names(input))){
      x <- split_input_control(input, control)

      input <- x$input
      control <- x$control
    }
  }

  if (is.list(control)){
    ctrl <- Biostrings::BStringSetList(control)
    control <- unlist(ctrl)
  }

  if (is(control, "BStringSetList")){
    control <- unlist(control)
  }

  return(list(input = input,
              control = control))
}
