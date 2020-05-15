#' Force best tomtom match by id
#'
#' @param res results from runTomTom
#' @param matches named vector where name is the input motif id, and value is the match_name to use as the new best match
#'
#' @return `res` with new best_* columns and re-ranked tomtom data in the `tomtom` list column for the updated entries.
#' @export
#'
#' @examples
#' # TODO: UPDATE THIS TO A WORKING EXAMPLE
#' motif <- universalmotif::create_motif()
#' db <- system.file("extdata/db/fly_factor_survey_id.meme")
#' res <- runTomTom(motif)
#' force_best_match(res, c("id" = "update"))
force_best_match <- function(res, matches){
  # matches is named vector, name = best_match
  purrr::iwalk(matches, ~{
    res[res$name == .y,]$tomtom[[1]] <<- res[res$name == .y,]$tomtom[[1]] %>%
      rank_tomtom_by_name(.x)
  })

  res %>%
    update_best_match
}

#' Update best match info by ranking of tomtom data
#'
#' This function updates the best_match columns based on the rankings on the tomtom list data.
#'
#' @param res results from runTomTom
#'
#' @return `res` with updated best_* columns
#' @export
#'
#' @examples
update_best_match <- function(res){
  res_nobest <- res %>%
    drop_best_match()

  new_tomtom <- res_nobest %>%
    tidyr::unnest(tomtom) %>%
    dplyr::select("name", "altname", dplyr::contains("match_"), "db_name") %>%
    nest_tomtom_results_best_top_row()

  res_nobest %>%
    dplyr::select(-"tomtom") %>%
    dplyr::left_join(new_tomtom, by = c("name", "altname")) %>%
    dplyr::select("name", "altname", dplyr::everything(),
                  dplyr::contains("best_"), "tomtom")
}

#' Drop best match info from tomtom results
#'
#' Convenience function for dropping all columns created by runTomTom prefixed
#' by "best_match_" and the "best_db_name" column. Keeps the "tomtom" data.frame
#' column. Can be useful if you want to unnest the data to resort the data
#'
#' @param res results of runTomTom
#'
#' @return
#' @export
#'
#' @examples
drop_best_match <- function(res){
  res %>%
    dplyr::select(-dplyr::contains("best_match_"),
                  -dplyr::any_of("best_db_name"))

}

#' Nest TomTom results columns into a data.frame column named "tomtom"
#'
#' This will also update the best_match information automatically to avoid any
#' ambiguities after manipulating unnested data. **NOTE:** that the resulting
#' columns may not be in the same order, so operations like `identical()` may
#' fail even though the column values are unchanged.
#'
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom magrittr %<>%
#' @examples
#' \dontrun{
#' #TODO: add better example
#' res <- runTomTom(motifs)
#' data <- tidyr::unnest(res, "tomtom")
#' identical(nest_tomtom(data), res)
#' }
nest_tomtom <- function(data){
  # Save motifs
  motif <- data$motif
  names(motif) <- data$name
  motif <- unique(motif)

  match_motif <- data$match_motif

  # tidyr::nest doesn't work with S4 because vctrs doesn't support it
  df <- data %>%
    dplyr::select(-dplyr::any_of(c("motif", "match_motif", "best_match_motif"))) %>%
    tidyr::nest(data = c("match_name",
                "match_altname",
                "match_pvalue",
                "match_evalue",
                "match_qvalue",
                "db_name")) %>%
    dplyr::rename("tomtom" = "data")

  # add back motifs
  df$motif <- motif
  # Add database motifs to tomtom
  i <- 1
  df$tomtom <- lapply(df$tomtom, function(x) {
    n <- nrow(x)
    j <- i + n - 1
    x$match_motif <- match_motif[i:j]
    i <<- j + 1
    return(x)
  })

  df %>%
    update_best_match()
}

#' Rank a tomtom results dataframe by match_name
#'
#' @param tomtom
#'
#' @param match_name
#'
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @importFrom dplyr desc
#'
#' @noRd
rank_tomtom_by_name <- function(tomtom, match_name){

  tomtom %>%
    dplyr::arrange(desc(match_name %in% !!match_name))

}
