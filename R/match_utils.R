#' Force best tomtom match by id
#'
#' @param res results from runTomTom
#' @param matches named vector where name is the input motif id, and value is the match_id to use as the new best match
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
  # matches is named vector, id = best_match
  purrr::iwalk(matches, ~{
    res[res$id == .y,]$tomtom[[1]] <<- res[res$id == .y,]$tomtom[[1]] %>%
      rank_tomtom_by_id(.x)
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

  res_notomtom <- res_nobest %>%

  new_tomtom <- res_nobest %>%
    tidyr::unnest(tomtom) %>%
    dplyr::select("id", "alt", dplyr::contains("match_"), "db_name") %>%
    nest_tomtom_results()

  res_nobest %>%
    dplyr::select(-"tomtom") %>%
    dplyr::left_join(new_tomtom, by = c("id", "alt"))

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
                  -"best_db_name")

}

#' Nest TomTom results columns into a data.frame column named "tomtom"
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #TODO: add better example
#' res <- runTomTom(motifs)
#' data <- tidyr::unnest(res, "tomtom")
#' identical(nest_tomtom(data), res)
#' }
nest_tomtom <- function(data){
  data %>%
    tidyr::nest(data = c("match_id",
                "match_alt",
                "match_pvalue",
                "match_evalue",
                "match_qvalue",
                "match_motif",
                "db_name")) %>%
    dplyr::rename("tomtom" = "data")
}

#' Rank a tomtom results dataframe by match_id
#'
#' @param tomtom
#'
#' @param match_id
#'
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @importFrom dplyr desc
#'
#' @noRd
rank_tomtom_by_id <- function(tomtom, match_id){

  tomtom %>%
    dplyr::arrange(desc(match_id %in% !!match_id))

}
