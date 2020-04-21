# Functions for validating/creating/etc. custom internal classes


#' Unfinished constructor for a dreme_results data.frame
#'
#' NOTE: `motifs` is only a list() not universalmotif list when it's empty.
#'
#' @return an empty dreme_results valid data.frame
#'
#' @examples
#' @noRd
new_dreme_results <- function(){
  data.frame(rank = integer(),
             id = character(),
             alt = character(),
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
             motifs = list(),
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
#' @examples
#' @noRd
error_dreme_results <- function(res){
  spec_dreme_res <- new_dreme_results()
  if (!(all(names(spec_dreme_res) %in% names(res)))) {

    missingNames <- names(spec_dreme_res)[!names(spec_dreme_res) %in% names(res)]
    nameString <- paste(missingNames, collapse = ", ")
    stop(paste0("object is not a valid dreme results data.frame. Missing columns: ", nameString))
  }

  is_universalmotif <- purrr::map_lgl(res$motifs, ~{class(.x) == "universalmotif"}) %>%
    purrr::set_names(NULL)

  if (length(is_universalmotif) == 0) stop("motifs column is empty")

  if (!all(is_universalmotif)) stop("not all objects in motif column are of type universalmotif")
  return(NULL)
}

#' Bool is/isnot valid dreme_results
#'
#' @param res data.frame
#'
#' @return TRUE or FALSE
#'
#' @examples
#' @noRd
is_dreme_results <- function(res){
  spec_dreme_res <- new_dreme_results()

  # all names exist
  if (!(all(names(spec_dreme_res) %in% names(res)))) {
    return(FALSE)
  }

  # motif column is universalmotif type
  is_universalmotif_list(res$motifs)
}

#' Bool is/isnot list of universalmotif objects
#'
#' @param list a list
#'
#' @return TRUE or FALSE if all members are universalmotif objects.
#'
#' @examples
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
#' @examples
#' @noRd
error_universalmotif_list <- function(list){

  if (is_universalmotif(list)) return(NULL)

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
#' @return data.frame with all motif metadata + `motifs` column containing universalmotif object
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' universalmotif::create_motif() %>%
#'   universalmotif_to_df()
#' @noRd
universalmotif_to_meme_df <- function(motif){
  # this is needed to overcome limitation of bind_rows causing error with list columns
  data <- switch(class(motif),
                 list = purrr::map(motif, universalmotif::as.data.frame) %>%
                          dplyr::bind_rows(),
                 universalmotif = universalmotif::as.data.frame(motif))

  df <- data %>%
    dplyr::rename("id" = "name",
                  "alt" = "altname")

  df$motifs <- motif
  return(df)
}
