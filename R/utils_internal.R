#' Returns valid path to meme bin directory or supported meme utility
#'
#' @param path user override path to meme bin (optional)
#' @param util one of c(NULL,"dreme","ame","fimo","tomtom")
#'
#' @return valid path to meme/bin or meme utility
#'
#' @examples
#'
#' @noRd
handle_meme_path <- function(path = NULL, util = NULL){
  f <- dotargs::build_path_handler(environment_var = "MEME_PATH",
                                       option_name = "meme_bin",
                                       default_path = "~/meme/bin",
                                       utils = c("dreme", "ame", "fimo", "tomtom")
                                       )
  f(path, util)
}

#' Helper for writing unique directory names
#'
#' @param input input file path
#' @param control control file path
#'
#' @return directory name in the style of: "input_vs_output".
#'
#' @examples
#' outdir_name("condition1.fa", "backgroundSequence.fa")
#' outdir_name("path/to/condition1.fa", "backgroundSequence.fa")
#'
#' @noRd
outdir_name <- function(input, control){

  paste0(dirname(input), "/",
         basename(tools::file_path_sans_ext(input)),
         "_vs_",
         basename(tools::file_path_sans_ext(control)))
}


#' Converts xml attrs to data-frame
#'
#' @param xml xml object
#'
#' @return
#'
#' @examples
#'
#' @noRd
attrs_to_df <- function(xml, ...) {
  # parsing DREME output
  # converts xml attributes to dataframe
  # where each column is an attribute
  xml2::xml_attrs(xml) %>%
    data.frame() %>%
    t() %>%
    data.frame(row.names = NULL, ...)
}


#' Writes universalmotif list to tempfile by default
#'
#' light wrapper on universalmotif::write_meme which returns path to file
#' written. Defaults to writing temporary file.
#'
#' @param list universalmotif list
#' @param path path to write
#' @param version meme version to append to header (default: 5)
#'
#' @return valid path
#'
#' @noRd
write_meme_list <- function(list, path = tempfile(fileext = ".meme"), version = 5){
  list %>%
    universalmotif::write_meme(path, append = F, overwrite = T, version = version)

  dotargs::check_files_exist(path)

  return(path)
}

#' Copy a file to temp location for testing
#'
#' I want to test using some preexisting files but don't want to update the git
#' history for them, so this copies to temp location.
#'
#' @param path
#'
#' @return tempfile path
#'
#' @noRd
duplicate_file <- function(path){
  dupFile <- tempfile()
  file.copy(path, dupFile)
  return(dupFile)
}

#' Check if processx process completed successfully
#'
#' @param processx_out output of processx::run(error_on_status = F)
#'
#' @return NULL if exit status 0, otherwise print all stdout + stderr
#'
#' @examples
#'
#' @noRd
process_check_error <- function(processx_out){
  if (processx_out$status != 0){
    cat(processx_out$stdout)
    stop(processx_out$stderr)
  }
  return(NULL)
}
