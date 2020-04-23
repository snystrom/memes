sequence_input.character <- function(input){
  # character input should be file path unless input == "shuffle"
  if (input == "shuffle") return(input)

  dotargs::check_files_exist(input)
  return(input)
}

sequence_input.DNAStringSet <- function(input){
  # stringset is written to temporary fasta file
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

  dotargs::check_files_exist(input)

  out <- list(metadata = NULL,
              path = input)
  return(out)
}

motif_input.data.frame <- function(input, path = tempfile(fileext = ".meme")){
  # check data.frame is dreme_out object
  if (!is_dreme_results(input)) warn_dreme_results(input)

  path <- input$motifs %>%
    write_meme_input_path(path = path)

  out <- list(metadata = input,
              path = path)
  return(out)
}

motif_input.list <- function(input, path = tempfile(fileext = ".meme")){
  # check list is universalmotif list
  if (!is_universalmotif_list(input)) error_universalmotif_list(list)

  df <- universalmotif_to_meme_df(input)

  path <- input %>%
    write_meme_input_path(path = path)

  out <- list(metadata = df,
              path = path)

  return(out)
}

motif_input.universalmotif <- function(input, path = tempfile(fileext = ".meme")){

  df <- universalmotif_to_meme_df(input) %>%
    data.frame

  path <- input %>%
    write_meme_input_path(path = path)

  out <- list(metadata = df,
              path = path)

  return(out)
}

#' Helper for exporting universalmotif data to temp files
#'
#' @param input
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
