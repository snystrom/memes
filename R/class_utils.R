# Functions for validating/creating/etc. custom internal classes


#' Unfinished constructor for a dreme_results data.frame
#'
#' NOTE: `motif` is only a list() not universalmotif list when it's empty.
#'
#' @return an empty dreme_results valid data.frame
#'
#' @noRd
new_dreme_results <- function(){
  data.frame(rank = integer(),
             name = character(),
             altname = character(),
             seq = character(),
             length = integer(),
             nsites = integer(),
             positive_hits = integer(),
             negative_hits = integer(),
             pvalue = numeric(),
             evalue = numeric(),
             unerased_evalue = numeric(),
             positive_total = integer(),
             negative_total = integer(),
             pos_frac = numeric(),
             neg_frac = numeric(),
             motif = list(),
             stringsAsFactors = F
             )
}

#' Validate & throw errors for dreme_results data.frames
#'
#' @param res data.frame to check if contains valid dreme_results
#'
#' @return NULL if successful, otherwise will return one of several errors
#'   describing why this object is not valid. Tries to be as helpful as possible
#'   to the user when describing properties that are incompatible with the data type.
#' @export
#'
#' @noRd
error_dreme_results <- function(res){
  spec_dreme_res <- new_dreme_results()
  if (!(all(names(spec_dreme_res) %in% names(res)))) {

    missingNames <- names(spec_dreme_res)[!names(spec_dreme_res) %in% names(res)]
    nameString <- paste(missingNames, collapse = ", ")
    stop(paste0("object is not a valid dreme results data.frame. Missing columns: ", nameString))
  }

  is_universalmotif <- purrr::map_lgl(res$motif, ~{class(.x) == "universalmotif"}) %>%
    purrr::set_names(NULL)

  if (length(is_universalmotif) == 0) stop("motif column is empty")

  if (!all(is_universalmotif)) stop("not all objects in motif column are of type universalmotif")
  return(NULL)
}

#' Bool is/isnot valid dreme_results
#'
#' @param res data.frame
#'
#' @return TRUE or FALSE
#'
#' @noRd
is_dreme_results <- function(res){
  spec_dreme_res <- new_dreme_results()

  # all names exist
  if (!(all(names(spec_dreme_res) %in% names(res)))) {
    return(FALSE)
  }

  # motif column is universalmotif type
  is_universalmotif_list(res$motif)
}

#' Check whether input is a universalmotif data.frame
#'
#' @param res data.frame
#'
#' @return TRUE or FALSE
#'
is_universalmotif_dataframe <- function(res){
  spec_df <- universalmotif::create_motif() %>%
    as_universalmotif_dataframe()

  if (!all(names(spec_df) %in% names(res))) {
    return(FALSE)
  }

  is_universalmotif_list(res$motif)
}

#' Bool is/isnot list of universalmotif objects
#'
#' @param list a list
#'
#' @return TRUE or FALSE if all members are universalmotif objects.
#'
#' @noRd
is_universalmotif_list <- function(list){
  if (length(list) == 0) {return(FALSE)}

  purrr::map_lgl(list, ~{class(.x) == "universalmotif"}) %>%
    purrr::set_names(NULL) %>%
    all
}

#' Validate & throw errors for universalmotif list
#'
#' @param list a list
#'
#' @return NULL or informative error describing why list is invalid
#'
#' @noRd
error_universalmotif_list <- function(list){

  if (is_universalmotif_list(list)) return(NULL)

  check_universalmotif <- purrr::map_lgl(list, ~{class(.x) == "universalmotif"}) %>%
    purrr::set_names(NULL)

  # warn no objects are motif
  if (sum(check_universalmotif) == 0) stop("no entries in list are of type universalmotif")

  # warn some objects not motif
  if (!sum(check_universalmotif) == length(check_universalmotif)) {
    bad_index <- which(!check_universalmotif)
    bad_index <- paste0("c(" , paste(bad_index, collapse = ", "), ")")
    stop(paste0("some entries in list are not of type universalmotif.\nIndices of bad entries: ", bad_index))
  }
}


#' Convert universalmotif to data.frame with motif column
#'
#' @param motif universalmotif object or list of universalmotifs
#'
#' @return data.frame with all motif metadata + `motif` column containing the universalmotif object
#'
#' @export
#'
#' @examples
#' universalmotif::create_motif() %>%
#'   as_universalmotif_dataframe()
as_universalmotif_dataframe <- function(motif, na.rm = FALSE){
  data <- universalmotif::summarise_motifs(motif, na.rm = na.rm)

  if (class(motif) == "universalmotif"){
    data$motif <- list(motif)
  } else if (class(motif) == "list"){
    data$motif <- motif
  }
  return(data)
}

#' Convert universalmotif data.frames back into universalmotifs
#'
#' @param data a universalmotif_dataframe (output from
#'   `as_universalmotif_dataframe()`, or `runDreme()`)
#'
#' @return universalmotif list from motifs, updated to reflect the data.frame
#'   column values.
#' @export
#'
#' @examples
#' df <- universalmotif::create_motif() %>%
#'   as_universalmotif_dataframe() %>%
#'   dplyr::mutate(altname = "new_alt_name")
#'
#' motifs <- as_universalmotif(df)
as_universalmotif <- function(data){
  data %<>%
    update_motifs()

  return(data$motif)
}
