#' Searches for valid TOMTOM database file
#'
#' Default search heirarchy: Sys.getEnv("TOMTOM_DB") > getOption("tomtom_db") > user-defined path
#'
#' @param path optional path to tomtom database
#'
#' @return valid path to tomtom database
#'
#' @examples
#'
#' @noRd
handle_tomtom_database_path <- function(path = NULL){
  f <- dotargs::build_path_handler(environment_var = "TOMTOM_DB",
                              option_name = "tomtom_db")
  f(path = path)
}

#' Run TomTom on target motifs
#'
#' @param input path to .meme format file of motifs or list of universalmotifs
#' @param database path to .meme format file to use as reference database (or list of universalmotifs)
#' @param outdir directory to store tomtom results (will be overwritten if
#'   exists). Default: location of input fasta file, or temporary location if using universalmotif list input.
#' @param thresh report matches less than or equal to this value. If evalue =
#'   TRUE (default), set an e-value threshold (default = 10). If evalue = FALSE,
#'   set a value between 0-1 (default = 0.5).
#' @param min_overlap only report matches that overlap by this value or more,
#'   unless input motif is shorter, in which case the shorter length is used as
#'   the minimum value
#' @param dist distance metric. Valid arguments: `allr | ed | kullback | pearson | sandelin | blic1 | blic5 | llr1 | llr5`.
#'   Default: `pearson`.
#' @param evalue whether to use E-value as significance threshold (default:
#'   `TRUE`). If evalue = FALSE, uses *q-value* instead.
#' @param meme_path path to "meme/bin/" (optional). If unset, will check R
#'   environment variable "TOMTOM_DB" (set in `.Renviron`), or option
#'   "tomtom_db" (set with `option(tomtom_db = "path/to/meme/bin")`)
#' @param ... additional flags passed to tomtom using {dotargs} formating (see
#'   [tomtom commandline reference](http://meme-suite.org/doc/tomtom.html?man_type=web) for details)
#'
#' @return data.frame of match results. Contains `match_motif` column of
#'   `universalmotif` objects with the matched PWM from the database.
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' runTomTom("searchmotifs.meme", "jasparMotifs.meme")
#' }
runTomTom <- function(input, database = NULL,
                      outdir = paste0(dirname(input), "/tomtom"),
                      thresh = 10,
                      min_overlap = 5,
                      dist = "pearson",
                      evalue = T,
                      meme_path = NULL, ...){
  # use TOMTOM server default values which differ from commandline defaults
  # TODO: email TOMTOM maintainers to ask if this is really a better default?

  command <- handle_meme_path(path = meme_path, util = "tomtom")

  database <- handle_tomtom_database_path(path = database)

  flags <- prepareTomTomFlags(outdir = outdir, thresh = thresh, min_overlap = min_overlap, dist = dist, evalue = evalue, ...)
  flags <- c(flags, input, database)

  processx::run(command, flags, spinner = T)

  tomtom_out <- dotargs::expected_outputs(c("tsv", "xml", "html"), "tomtom", outdir = outdir)

  tomtom_out %>%
    dotargs::check_files_exist()

  tomtom_results <- parseTomTom(tomtom_out$xml)

  return(tomtom_results)
}

#' Generate flags for TomTom commandline
#'
#' see runTomTom for explanation of inputs.
#'
#' Here I set the function arguments because the TOMTOM web server uses
#' different defaults from the commandline. I want to emulate the web server
#' function so need to always set these args to non-default values.
#'
#' @param outdir outdir
#' @param thresh threshold
#' @param min_overlap minimum overlap
#' @param dist distance function
#' @param evalue evalue
#' @param ... interpreted as {dotargs}
#'
#' @return
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
prepareTomTomFlags <- function(outdir, thresh, min_overlap, dist, evalue, ...){
  # lookup table converts arguments with - to _ so
  # user doesn't have to escape flags
  argsDict = c("outdir" = "oc",
               "min_overlap" = "min-overlap",
               "motif_pseudo" = "motif-pseudo",
               "no_ssc" = "no-ssc",
               "incomplete_scores" = "incomplete-scores")

  flags <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags(argsDict) %>%
    dotargs::crystallize_flags()

  return(flags)
}

#' Return query table from tomtom xml
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return data.frame of query data
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
get_tomtom_query_data <- function(tomtom_xml_data){
  xml2::xml_find_all(tomtom_xml_data, "//queries") %>%
    xml2::xml_children() %>%
    attrs_to_df() %>%
    dplyr::mutate(query_idx = (1:nrow(.) - 1)) %>%
    dplyr::rename("db_idx" = "db")
}

#' Get match table from tomtom xml
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return data.frame of match data
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
get_tomtom_match_data <- function(tomtom_xml_data){
  matches <- xml2::xml_find_all(tomtom_xml_data, "//matches") %>%
    xml2::xml_children()

  match_df <- purrr::map(matches, xml2::xml_children) %>%
    purrr::set_names(xml2::xml_attr(matches, "idx")) %>%
    purrr::map_dfr(attrs_to_df, stringsAsFactors = F, .id = "query_idx") %>%
    dplyr::mutate_at(c("pv", "ev", "qv"), as.double) %>%
    dplyr::mutate_at(c("query_idx", "idx", "off"), as.integer) %>%
    dplyr::mutate_at("rc", as.character()) %>%
    dplyr::rename("offset" = "off",
                  "pvalue" = "pv",
                  "evalue" = "ev",
                  "qvalue" = "qv",
                  "target_idx" = "idx") %>%
    dplyr::mutate(strand = ifelse(rc == "y", "-", "+"))

  return(match_df)
}

#' Get database info
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return information about each database file used in TomTom call
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
get_tomtom_db_data <- function(tomtom_xml_data){
  xml2::xml_find_all(tomtom_xml_data, "//target_dbs") %>%
    xml2::xml_children() %>%
    attrs_to_df(stringsAsFactors = F) %>%
    dplyr::mutate(db_idx = (1:nrow(.)) - 1) %>%
    dplyr::rename("db_name" = "name")
}

#' Get Target info & PWM
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return data.frame with target motif stats and a column "match_motif"
#'   containing a universalmotif format object for that motif.
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
get_tomtom_target_data <- function(tomtom_xml_data){
  targets <- xml2::xml_find_all(tomtom_xml_data, "//targets") %>%
    xml2::xml_children()

  target_df <- targets %>%
    attrs_to_df(stringsAsFactors = F) %>%
    dplyr::rename("db_idx" = "db") %>%
    dplyr::mutate_at(c("length", "nsites", "db_idx"), as.integer) %>%
    dplyr::mutate(target_idx = (1:nrow(.)) - 1)

  target_pfms <- targets %>%
    purrr::map(get_probability_matrix) %>%
    purrr::map(t)

  target_list <- target_df %>%
    split(.$id)

  pfmList <- purrr::map2(target_list, target_pfms, ~{
    universalmotif::create_motif(.y,
                                 type = "PCM",
                                 name = .x$id,
                                 altname = .x$alt,
                                 nsites = .x$nsites)

  })

  target_data <- target_df %>%
    dplyr::rename_at(c("id", "alt"), ~{paste0("match_", .x)}) %>%
    dplyr::select(dplyr::contains("idx"), dplyr::contains("match")) %>%
    dplyr::mutate(match_motif = pfmList)

  return(target_data)

}

#' Return tomtom stat information and PWMs for matched motifs
#'
#' @param tomtom_xml_path path to tomtom.xml output
#'
#' @return data.frame:
#'     - id: name of query PWM
#'     - alt: alternate name of query PWM
#'     - match_id: name of matched PWM
#'     - match_alt: alt name of matched PWM
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
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' parseTomTom("tomtom.xml")
#' }
parseTomTom <- function(tomtom_xml_path){
  tt_xml <- xml2::read_xml(tomtom_xml_path)

  query_data <- tt_xml %>%
    get_tomtom_query_data() %>%
    dplyr::select("query_idx", "id", "alt")

  match_data <- tt_xml %>%
    get_tomtom_match_data() %>%
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
    dplyr::select("id", "alt", "match_id", "match_alt", dplyr::contains("value"), "db_name", "match_motif")

  return(tomtom_results)

}
