# import XML

#' Import tomtom data from previous run
#'
#' @param tomtom_xml_path path to tomtom.xml
#'
#' @return will return data.frame with input motifs & results for best match.
#'   `tomtom` list column contains full tomtom data for each input motif.
#'   NOTE: if tomtom detects no matches for any input motif, currently will
#'   throw warning & return NA values for `tomtom`, `best_match_id`, and
#'   `best_match_motifs`.
#' @export
#'
#' @examples
#' \dontrun{
#' importTomTomXML("path/to/tomtom.xml")
#' }
importTomTomXML <- function(tomtom_xml_path){
  tt_xml <- xml2::read_xml(tomtom_xml_path)

  query_data <- tt_xml %>%
    tomtom_query_motif_dfs()

  match_data <- tt_xml %>%
    get_tomtom_match_data()

  if (is.null(match_data)) {
    warning("TomTom detected no matches")
    #TODO: handle NULL w/ return data w/ NA for all tomtom columns,
    # this way hopefully it will break fewer pipelines
    query_data %>%
      dplyr::mutate(
        best_match_id = NA,
        best_motifs = NA,
        tomtom = NA) %>%
      return()
  }

  match_data %<>%
    dplyr::rename_at(c("offset", "pvalue", "evalue", "qvalue", "strand"), ~{paste0("match_", .x)}) %>%
    dplyr::select(-"rc")

  target_db_lookup <- tt_xml %>%
    get_tomtom_db_data %>%
    dplyr::select(db_idx, db_name)

  target_data <- get_tomtom_target_data(tt_xml) %>%
      dplyr::left_join(target_db_lookup, by = "db_idx")

  tomtom_results <- query_data %>%
    dplyr::left_join(match_data, by = "query_idx") %>%
    dplyr::left_join(target_data, by = "target_idx") %>%
    dplyr::select(-dplyr::contains("idx"))

  res <- dplyr::left_join(query_data, nest_tomtom_results(tomtom_results), by = c("id", "alt"))

  return(res)

}


#' Import Dreme output from previous run
#'
#' @param dreme_xml_path path to dreme.xml file
#'
#' @return data.frame with statistics for each discovered motif. The `motifs`
#'   column contains a universalmotif object representation in PCM format of
#'   each DREME motif. If no motifs are discovered, returns NULL.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' importDremeXML("dreme.xml")
#' }
importDremeXML <- function(dreme_xml_path){
  parseDreme(dreme_xml_path)
}
