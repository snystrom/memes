# import XML

#' Import tomtom data from previous run
#'
#' @param tomtom_xml_path path to tomtom.xml
#'
#' @return will return data.frame with input motifs & results for best match.
#'   `tomtom` list column contains full tomtom data for each input motif.
#'   NOTE: if tomtom detects no matches for any input motif, currently will
#'   throw warning & return NA values for `tomtom`, `best_match_name`, and
#'   `best_match_motifs`.
#'
#' @details tomtom list column format
#' the `tomtom` list column contains data.frames with the following format:
#'     - name: name of query PWM
#'     - altname: alternate name of query PWM
#'     - match_name: name of matched PWM
#'     - match_altname: alt name of matched PWM
#'     - match_pvalue: p-value of match
#'     - match_evalue: E-value of match
#'     - match_qvalue: q-value of match
#'     - match_offset: number of letters the query was offset from the target match
#'     - match_strand: whether the motif was found on input strand (+) or as reverse-complement (-)
#'     - db_name: database source of matched motif
#'     - match_motif: universalmotif object containing the PWM that was matched
#'
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
    message("TomTom detected no matches")
    #TODO: handle NULL w/ return data w/ NA for all tomtom columns,
    # this way hopefully it will break fewer pipelines
    null_match <- query_data %>%
      dplyr::mutate(
        best_match_name = NA,
        best_motifs = NA,
        tomtom = NA)

    return(null_match)
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
    dplyr::select(-dplyr::contains("idx")) %>%
    # Rename columns for max compatibility with universalmotif
    dplyr::rename("match_name" = "match_id",
                  "match_altname" = "match_alt")

  res <- dplyr::left_join(query_data, nest_tomtom_results(tomtom_results), by = c("name", "altname"))

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
