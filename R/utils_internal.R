utils::globalVariables(".")

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
  f <- cmdlr::build_path_handler(environment_var = "MEME_BIN",
                                       option_name = "meme_bin",
                                       default_path = "~/meme/bin",
                                       utils = c("dreme", "ame", "fimo", "tomtom", "meme")
                                       )
  f(path, util)
}

#' Searches for valid MEME database file
#'
#' Default search heirarchy: Sys.getEnv("MEME_DB) > getOption("meme_db") > user-defined path
#'
#' @param path optional path to tomtom database (either `character(1)` or
#'   `list()` if named list, use names as database name)
#'
#' @return valid path to tomtom database
#'
#' @examples
#'
#' @noRd
handle_meme_database_path <- function(path = NULL){
  # database can be path, or universalmotif list, (or vector: c(motifList, path))
  # names will be used as file names for non file-path entries

  if (!is.null(path)){
    if (path == ""){
      stop("path cannot be an empty string")
    }
  }

  if (any(is.data.frame(path))) {
    stop("data.frame is not a supported input type, if this is a dreme results object, try passing it inside a list like: database = list(results)")
  }

  if (!is.character(path) & !is.list(path) & !is.null(path)){
    stop("path must be character or list")
  }

  if (length(path) > 1 | is.list(path)){
    paths <- purrr::imap(path, ~{
      # Resolve how to name database entries:
      if(.y != "" & !is.character(.x)) {
        # rename non-path inputs to index# or name (if defined by user)
        out <- file.path(tempdir(), .y)
      } else if (.y != "" & !is.numeric(.y)) {
        # use current file name & path if user does not define a new name for path entries
        # (allows path inputs when all entries unnamed to not get renamed to their index position)
        out <- file.path(tempdir(), .y)
      } else{
        # Otherwise, use type-specific path default
        out <- NULL
      }

      motif_input(.x, out)
    }) %>%
      purrr::map_chr("path") %>%
      purrr::set_names(NULL)
    return(paths)
  }


  # Allows setting option to a universalmotif object
  # and return path
  if (is.null(path) & !is.null(getOption("meme_db"))) {
    if (!is.character(getOption("meme_db"))){
      db <- getOption("meme_db")
      x <- motif_input(db)
      return(x$path)
    }
  }

  # Otherwise search envrionment variable / option definition
  f <- cmdlr::build_path_handler(environment_var = "MEME_DB",
                                   option_name = "meme_db")
  f(path = path)
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

  cmdlr::check_files_exist(path)

  return(path)
}

#' Copy a file to temp location for testing
#'
#' I want to test using some preexisting files but don't want to update the git
#' history for them, so this copies to temp location.
#'
#' @param path path to file to duplicate to tempfile
#'
#' @return tempfile path
#'
#' @noRd
duplicate_file <- function(path){
  dupFile <- tempfile()
  file.copy(path, dupFile)
  return(dupFile)
}


#' Normalize rank
#'
#' For groups of different size, it is inappropriate to compare rank position in a heatmap, for example
#'
#' This function converts rank as a fraction of total ranks for better between-group comparisons.
#'
#' @param rank
#'
#' @return
#'
#' @examples
#' rank_normalize(c(1,3,5))
#'
#' @noRd
rank_normalize <- function(rank){
  if (length(rank) == 1) {
    # Rank 1 is highest rank
    return(0)
  }
  (rank - 1) / (max(rank) - 1)
}
