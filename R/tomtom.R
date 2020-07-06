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
#'   environment variable "MEME_DB (set in `.Renviron`), or option
#'   "meme_db" (set with `option(meme_db = "path/to/meme/bin")`)
#' @param ... additional flags passed to tomtom using {cmdr} formating (see table below for details)
#'
#' @details Additional arguments
#'
#'  runTomTom() can accept all valid tomtom arguments passed to `...` as described in the
#'  [tomtom commandline reference](http://meme-suite.org/doc/tomtom.html?man_type=web). For
#'  convenience, below is a table of valid arguments, their default values, and
#'  their description.
#'
#'
#' | TomTom Flag       | allowed values | default | description                |
#' |:-----------------:|:--------------:|:-------:|:---------------------------|
#' | bfile             | file path      | `NULL`  | path to background model for converting frequency matrix to log-odds score (not used when `dist` is set to "ed", "kullback", "pearson", or "sandelin" |
#' | motif_pseudo      | `numeric`      | 0.1     | pseudocount to add to motifs |
#' | xalph             | `logical`      | FALSE   | convert alphabet of target database to alphabet of query database |
#' | norc              | `logical`      | FALSE   | Do not score reverse complements of motifs |
#' | incomplete_scores | `logical`      | FALSE   | Compute scores using only aligned columns |
#' | thresh            | `numeric`      | 0.5     | only report matches with significance values <= this value. Unless `evalue = TRUE`, this value must be < 1. |
#' | internal          | `logical`      | FALSE   | forces the shorter motif to be completely contained in the longer motif |
#' | min_overlap       | `integer`      | 1       | only report matches that overlap by this number of positions or more. If query motif is smaller than this value, its width is used as the min overlap for that query |
#' | time              | `integer`      | `NULL`  | Maximum runtime in CPU seconds (default: no limit) |
#'
#'
#' @return data.frame of match results. Contains `match_motif` column of
#'   `universalmotif` objects with the matched PWM from the database. If no
#'   matches are returned, `tomtom` and `best_match_motif` columns will be set
#'   to `NA` and a message indicating this will print.
#'
#' @details # Citation
#' If you use `runTomTom()` in your analysis, please cite:
#'
#' Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford
#' Noble, "Quantifying similarity between motifs", Genome Biology, 8(2):R24,
#' 2007. [full text](http://genomebiology.com/2007/8/2/R24)
#'
#' @details ## Licensing
#' The MEME Suite is free for non-profit use, but for-profit users should purchase a
#' license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
#'
#' @export
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' runTomTom("searchmotifs.meme", "jasparMotifs.meme")
#' }
runTomTom <- function(input, database = NULL,
                      outdir = "auto",
                      thresh = 10,
                      min_overlap = 5,
                      dist = "pearson",
                      evalue = T,
                      meme_path = NULL, ...){
  # use TOMTOM server default values which differ from commandline defaults
  # TODO: email TOMTOM maintainers to ask if this is really a better default?

  # save dreme results & join w/ tomtom results at end.
  # type validation happens below

  input <- motif_input(input)

  if (is.null(input$metadata)){
    # Allows .meme input files to
    # import query motif metadata
    # I use this solution instead of modifying motif_input
    # to allow motif_input on databases to not require importing the file
    # since these can be large
    input$metadata <- input$path %>%
      universalmotif::read_meme() %>%
      as_universalmotif_dataframe()
  }

  command <- handle_meme_path(path = meme_path, util = "tomtom")

  if (outdir == "auto") {outdir <- file.path(dirname(input$path), "tomtom")}

  database <- handle_meme_database_path(path = database)

  user_flags <- prepareTomTomFlags(outdir = outdir,
                                   thresh = thresh,
                                   min_overlap = min_overlap,
                                   dist = dist,
                                   evalue = evalue, ...)
  flags <- c(user_flags, input$path, database)

  ps_out <- processx::run(command, flags, spinner = T, error_on_status = F)
  ps_out %>%
    process_check_error(help_fun = ~{tomtom_help(command)},
                        user_flags = cmdr::cmd_help_parse_flags(user_flags),
                        flags_fun = ~{gsub("-", "_", .)},
                        default_help_fun = TRUE
                        )

  tomtom_out <- cmdr::cmd_output_expect(c("tsv", "xml", "html"), "tomtom", outdir = outdir)

  tomtom_out %>%
    cmdr::cmd_files_exist()

  tomtom_results <- parseTomTom(tomtom_out$xml)

  if (!is.null(input$metadata) & is.null(tomtom_results)) {
    warning("TomTom returned no matches")
    # TODO: add every tomtom results column w/ NA values?
    input$metadata$best_match_motif <- rep(NA, nrow(input$metadata))
    input$metadata$best_match_name <- rep(NA_character_, nrow(input$metadata))
    input$metadata$tomtom <- rep(NA, nrow(input$metadata))
    return(input$metadata)
  }

  if (!is.null(input$metadata) & !is.null(tomtom_results)){

    nest_tomtom <- nest_tomtom_results(tomtom_results)
    merge_res <- dplyr::left_join(input$metadata, nest_tomtom, by = c("name", "altname"))
    return(merge_res)
  }

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
#' @param ... interpreted as {cmdr}
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

  flags <- cmdr::cmd_args_all() %>%
    cmdr::cmd_args_to_flags(argsDict) %>%
    cmdr::cmd_list_crystallize()

  return(flags)
}

#' Returns tomtom help lines
#'
#' @param command path to tomtom. output of handle_meme_path(util = "tomtom")
#'
#' @return
#'
#' @noRd
tomtom_help <- function(command){
  processx::run(command, error_on_status = FALSE)$stderr
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
    attrs_to_df(stringsAsFactors = FALSE) %>%
    dplyr::mutate(query_idx = (1:nrow(.) - 1)) %>%
    dplyr::rename("db_idx" = "db")
}

#' Get match table from tomtom xml
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return data.frame of match data or NULL if no matches found
#'
#' @examples
#'
#' @importFrom magrittr %>%
#'
#' @noRd
get_tomtom_match_data <- function(tomtom_xml_data){
  matches <- xml2::xml_find_all(tomtom_xml_data, "//matches") %>%
    xml2::xml_children()

  if (length(matches) == 0){return(NULL)}

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

  target_df$pfm <- target_pfms

  target_data <- target_df %>%
    dplyr::mutate(match_motif = purrr::pmap(list(pfm, id, alt, nsites), ~{
      universalmotif::create_motif(..1,
                                   type = "PCM",
                                   name = ..2,
                                   altname = ..3,
                                   nsites = ..4)
    })) %>%
    dplyr::select(-"pfm") %>%
    dplyr::rename_at(c("id", "alt"), ~{paste0("match_", .x)}) %>%
    dplyr::select(dplyr::contains("idx"), dplyr::contains("match"))

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
#'     Returns NULL if not matches detected
#'
#'
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @examples
#' \dontrun{
#' parseTomTom("tomtom.xml")
#' }
#' @noRd
parseTomTom <- function(tomtom_xml_path){
  tt_xml <- xml2::read_xml(tomtom_xml_path)

  query_data <- tt_xml %>%
    get_tomtom_query_data() %>%
    dplyr::select(dplyr::any_of(c("query_idx", "id", "alt")))

  match_data <- tt_xml %>%
    get_tomtom_match_data()

  if (is.null(match_data)) {return(NULL)}

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
    # any_of() protects against alt not existing
    dplyr::select(dplyr::any_of(c("id", "alt")), "match_id", "match_alt", dplyr::contains("value"), "db_name", "match_motif") %>%
    # Rename columns for max compatibility with universalmotif
    dplyr::rename("name" = "id",
                  "match_name" = "match_id",
                  "match_altname" = "match_alt") %>%
    # allows renaming alt column only if exists
    dplyr::rename_all(dplyr::recode, alt = "altname")

  if (!("altname" %in% names(tomtom_results))) {
    # Instantiate altname column w/ NA values if not exists
    tomtom_results %<>%
      dplyr::mutate(altname = NA_character_)
  }

  return(tomtom_results)

}

#' Nest tomtom results & show only best match, store all others in `tomtom` list column
#'
#' @param tomtom_results
#'
#' @return data.frame with columns w/ all data for "best" match (defined by top
#'   hit, lowest pvalue). All other matches are nested into 'tomtom' column.
#'   Which is list of data.frames for each match too the given id.
#'
#' @noRd
nest_tomtom_results <- function(tomtom_results){

  tomtom_results %>%
    dplyr::group_by(name, altname) %>%
    tidyr::nest() %>%
    dplyr::mutate(best_match_info = purrr::map(data, ~{
      .x %>%
        dplyr::filter(match_evalue == min(match_evalue)) %>%
        head(1) %>%
        dplyr::rename_all(~{paste0("best_", .x)})
    })) %>%
    tidyr::unnest(best_match_info) %>%
    dplyr::rename("tomtom" = "data") %>%
    dplyr::mutate(tomtom = purrr::map(tomtom, data.frame)) %>%
    dplyr::select("name", "altname", dplyr::contains("best_"), dplyr::everything())
}


#' Return best_match data for top row of `tomtom`
#'
#' @param tomtom_results tomtom results object
#'
#' @return data.frame with columns w/ all data for "best" match (defined by first row of `tomtom` data).
#'   All other matches are nested into 'tomtom' column.
#'   Which is list of data.frames for each match too the given id.
#'
#' @noRd
nest_tomtom_results_best_top_row <- function(tomtom_results){
  tomtom_results %>%
    dplyr::group_by(name, altname) %>%
    tidyr::nest() %>%
    dplyr::mutate(best_match_info = purrr::map(data, ~{
      .x[1,] %>%
        dplyr::rename_all(~{paste0("best_", .x)})
    })) %>%
    tidyr::unnest(best_match_info) %>%
    dplyr::rename("tomtom" = "data") %>%
    dplyr::mutate(tomtom = purrr::map(tomtom, data.frame)) %>%
    dplyr::select("name", "altname", dplyr::contains("best_"), dplyr::everything())
}
