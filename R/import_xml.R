#' Import tomtom data from previous run
#'
#' @param tomtom_xml_path path to tomtom.xml
#'
#' @return will return data.frame with input motifs & results for best match.
#'   `tomtom` list column contains full tomtom data for each input motif.
#'   NOTE: if tomtom detects no matches for any input motif, currently will
#'   print a message & return NA values for `tomtom`, `best_match_name`, and
#'   `best_match_motif`.
#'
#' @seealso [runTomTom()]
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
#' tomtom_xml <- system.file("extdata/tomtom.xml", package = "memes")
#' importTomTomXML(tomtom_xml)
importTomTomXML <- function(tomtom_xml_path){
  parseTomTom(tomtom_xml_path, query_metadata = NULL)
}


#' Import Dreme output from previous run
#'
#' @param dreme_xml_path path to dreme.xml file
#'
#' @return data.frame with statistics for each discovered motif. The `motifs`
#'   column contains a universalmotif object representation in PCM format of
#'   each DREME motif. If no motifs are discovered, returns NULL.
#'
#' @seealso [runDreme()]
#' @export
#'
#' @examples
#' dreme_xml <- system.file("extdata/dreme.xml", package = "memes")
#' importDremeXML(dreme_xml)
importDremeXML <- function(dreme_xml_path){
  parseDreme(dreme_xml_path)
}
