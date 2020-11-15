# Functions for validating/creating/etc. custom internal classes

#' Unfinished constructor for a tomtom_results data.frame
#'
#' NOTE: `motif` is only a list() not universalmotif list when it's empty.
#'
#' @return a minimal empty tomtom_results valid data.frame
#'
#' @noRd
new_tomtom_results <- function(){
  data.frame(name = NA_character_,
             altname = NA_character_,
             motif = NA,
             best_match_name = NA_character_,
             best_match_altname = NA_character_,
             best_match_offset = NA_integer_,
             best_match_pvalue = NA_real_,
             best_match_evalue = NA_real_,
             best_match_qvalue = NA_real_,
             best_match_strand = NA_character_,
             best_match_motif = NA,
             tomtom = NA,
             stringsAsFactors = FALSE
             )

}

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
             stringsAsFactors = FALSE
             )
}

#' Validate & throw errors for dreme_results data.frames
#'
#' @param res data.frame to check if contains valid dreme_results
#'
#' @return NULL if successful, otherwise will return one of several errors
#'   describing why this object is not valid. Tries to be as helpful as possible
#'   to the user when describing properties that are incompatible with the data type.
#'
#' @noRd
error_dreme_results <- function(res){
  spec_dreme_res <- new_dreme_results()
  if (!(all(names(spec_dreme_res) %in% names(res)))) {

    missingNames <- names(spec_dreme_res)[!names(spec_dreme_res) %in% names(res)]
    nameString <- paste(missingNames, collapse = ", ")
    stop(paste0("object is not a valid dreme results data.frame. Missing columns: ", nameString))
  }

  is_universalmotif <- purrr::map_lgl(res$motif, ~{is(.x, "universalmotif")}) %>%
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

#' Bool is/isnot valid tomtom_results
#'
#' @param res data.frame
#'
#' @return TRUE or FALSE
#'
#' @noRd
is_tomtom_results <- function(res){
  spec_tomtom_res <- new_tomtom_results()

  # all names exist
  if (!(all(names(spec_tomtom_res) %in% names(res)))) {
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

  purrr::map_lgl(list, ~{is(.x, "universalmotif")}) %>%
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

  check_universalmotif <- purrr::map_lgl(list, ~{is(.x, "universalmotif")}) %>%
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
#' The universalmotif data.frame structure is an R data frame where
#' `universalmotif` metadata slots are represented as columns, and individual
#' motifs are stored along rows. Additionally, the `universalmotif`
#' representation of the motif is stored in a nested column named `motif`. Users
#' can perform arbitrary manipulations to these data.frames such as adding,
#' removing, or renaming columns. When manipulating columns that correspond to
#' slot names of a `universalmotif` object, this stages changes which can be
#' propagated to the `universalmotif` object in the `motif` column by calling
#' [update_motifs()]. Users can convert the data.frame back to pure
#' `universalmotif` format using [as_universalmotif()].
#'
#' Columns which are linked to `universalmotif` format are:
#' `name`, `altname`, `family`, `organism`, `consensus`, `alphabet`, `strand`,
#' `icscore`, `nsites`, `bkgsites`, `pval`, `qval`, `eval`
#' 
#' Note that changing the above columns will result in changes to the
#' `universalmotif` representation when calling [update_motifs()] or
#' [as_universalmotif()]
#'
#' @param motif universalmotif object or list of universalmotifs
#' @param na.rm whether to include undefined columns for empty `universalmotif` slots (default: FALSE).
#'
#' @return data.frame with all motif metadata as columns and a special `motif`
#'   column containing the universalmotif object representation of each motif.
#'
#' @seealso [update_motifs()] for synchronizing the data.frame values with the
#'   `motif` column, and [as_universalmotif()] to convert back to `universalmotif` format.
#' 
#' @export
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' motif_df <- as_universalmotif_dataframe(motif)
as_universalmotif_dataframe <- function(motif, na.rm = FALSE){
  data <- universalmotif::summarise_motifs(motif, na.rm = na.rm)

  if (is(motif, "universalmotif")){
    data$motif <- list(motif)
  } else if (is(motif, "list")){
    data$motif <- motif
  }
  return(data)
}

#' Convert universalmotif data.frames back into universalmotifs
#' 
#' This function converts universalmotif data.frames into `universalmotif`
#' format, first by updating the `motif` metadata to reflect the current values
#' of the linked columns, then extracting the updated `universalmotif` objects.
#' Columns which do not correspond to `universalmotif` slot names are dropped
#' and not propagated to the `universalmotif` output.
#' 
#' Columns which are propagated to `universalmotif` format:
#' `name`, `altname`, `family`, `organism`, `consensus`, `alphabet`, `strand`,
#' `icscore`, `nsites`, `bkgsites`, `pval`, `qval`, `eval`
#'
#' @param data a universalmotif_dataframe (output from
#'   `as_universalmotif_dataframe()`, or `runDreme()`)
#'
#' @return universalmotif list from motifs, updated to reflect the data.frame
#'   column values.
#'   
#' @seealso [as_universalmotif_dataframe()]
#' @export
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' df <- as_universalmotif_dataframe(motif)
#' df <- dplyr::mutate(df, altname = "new_alt_name")
#'
#' motifs <- as_universalmotif(df)
as_universalmotif <- function(data){
  data %<>%
    update_motifs()

  return(data$motif)
}
