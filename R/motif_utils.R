#' Update motif columns by data values
#'
#' **NOTE** ignores NA values
#'
#' **NOTE** This feature is experimental and subject to change based on
#' user-feedback. Please provide feedback at
#' \link{https://github.com/snystrom/memes/issues/31}
#'
#' @param .data data.frame of results
#' @param ... quoted name-value pairs of columns to edit. format: "motif_column_name" = "data_column_name"
#' @param .motif column of universalmotif objects to edit (default: "motif")
#' @param .override named list or character vector where names are
#'   motif_column_name, values are data_column_name, used instead of ... (useful for programatically passing values)
#'
#' @return `.data` where `.motif` entries have been updated to values passed to `...` or `.override`
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' df <- as_universalmotif_dataframe(motif)
#' df <- dplyr::mutate(df, id = "newName")
#' df <- mutate_motif(df, "name" = "id")
#' # renamed motif
#' df$motif
#' @noRd
mutate_motif <- function(.data, ..., .motif = "motif", .override = NULL){
  #dots <- enquos(...)
  #return(dots)

  # motif col must exist & be list
  stopifnot(.motif %in% names(.data))
  stopifnot(is_universalmotif_list(.data[[.motif]]))

  dots <- cmdr::cmd_args_dots()

  if (!is.null(.override)){
    dots <- as.list(.override)
  }

  # check all values exist
  stopifnot(unlist(dots) %in% names(.data))
  stopifnot(names(dots) %in% names(universalmotif::summarise_motifs(.data[[.motif]], na.rm = FALSE)))

  # foreach entry, replace for each motif
  purrr::imap(dots, ~{
    data_col <- .x
    motif_slot <- .y

    # i tracks rows in dataframe
    i <- 1
    .data[[.motif]] <<- purrr::map(.data[[.motif]], ~{
      value <- .data[i,data_col]

      # Skip NA value replacement
      if (is.na(value)) {
        i <<- i+1
        return(.x)
      }

      .x[motif_slot] <- value
      i <<- i+1
      return(.x)
    })

  })
  return(.data)

}

#' Update the `motif` column to data.frame values
#'
#' @param .data data.frame with `motif` column
#'
#' @return .data where `motif` column has been updated to reflect the values
#'   from columns sharing names with unprotected universalmotif slots. Names of
#'   `motif` list are updated to reflect name.
#' @export
#'
#' @details
#'
#' **NOTE** that `consensus`, `alphabet`, `multifreq`, and `icscore` are protected columns and
#' cannot be updated.
#'
#' ## Table of values updated
#'
#' | `motif`  | `data.frame` |
#' |:--------:|:------------:|
#' | name     | name         |
#' | altname  | altname      |
#' | family   | family       |
#' | organism | organism     |
#' | strand   | strand       |
#' | nsites   | nsites       |
#' | bkgsites | bkgsites     |
#' | pval     | pval         |
#' | qval     | qval         |
#' | eval     | eval         |
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' df <- as_universalmotif_dataframe(motif)
#' df <- dplyr::mutate(df, id = "newName")
#' df <- update_motifs(df)
#' # renamed motif
#' df$motif
update_motifs <- function(.data){
  names_lookup <- c("name" = "name",
                    "altname" = "altname",
                    "family" = "family",
                    "organism" = "organism",
                    "strand" = "strand",
                    "nsites" = "nsites",
                    "bkgsites" = "bkgsites",
                    "pval" = "pval",
                    "qval" = "qval",
                    "eval" = "eval"
                    )

  to_mutate <- names_lookup[names_lookup %in% names(.data)]

  .data %<>%
    mutate_motif(.override = to_mutate)

  names(.data$motif) <- .data$name

  return(.data)
}
