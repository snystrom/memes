#' Run TomTom on target motifs
#'
#' TomTom compares input motifs to a database of known, user-provided motifs to
#' identify matches.
#'
#' runTomTom will rank matches by significance and return a
#' best match motif for each input (whose properties are stored in the `best_match_*`
#' columns) as well as a ranked list of all possible matches stored in the
#' `tomtom` list column.
#'
#' @param input path to .meme format file of motifs, a list of universalmotifs,
#'   or a universalmotif data.frame object (such as the output of `runDreme()`)
#' @param database path to .meme format file to use as reference database (or list of universalmotifs)
#' @param outdir directory to store tomtom results (will be overwritten if
#'   exists). Default: location of input fasta file, or temporary location if using universalmotif input.
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
#' @param ... additional flags passed to tomtom using {cmdfun} formating (see table below for details)
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
#' @return data.frame of match results. Contains `best_match_motif` column of
#'   `universalmotif` objects with the matched PWM from the database, a series
#'   of `best_match_*` columns describing the TomTom results of the match, and a
#'   `tomtom` list column storing the ranked list of possible matches to each
#'   motif. If a universalmotif data.frame is used as input, these columns are
#'   appended to the data.frame. If no matches are returned, `tomtom` and
#'   `best_match_motif` columns will be set to `NA` and a message indicating
#'   this will print.
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
#' @rdname runTomTom
#'
#' @md
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' if (meme_is_installed()) {
#' motif <- universalmotif::create_motif("CCRAAAW")
#' database <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes")
#' 
#' runTomTom(motif, database)
#' }
runTomTom <- function(input, database = NULL,
                      outdir = "auto",
                      thresh = 10,
                      min_overlap = 5,
                      dist = "pearson",
                      evalue = TRUE,
                      meme_path = NULL, ...){
  UseMethod("runTomTom")
}

#' @export
#' @noRd
runTomTom.list <- function(input, database = NULL,
                    outdir = "auto",
                    thresh = 10,
                    min_overlap = 5,
                    dist = "pearson",
                    evalue = TRUE,
                    meme_path = NULL, ...){
  purrr::map(input, 
             runTomTom.default, 
              database = database,
              outdir = outdir,
              thresh = thresh,
              min_overlap = min_overlap,
              dist = dist,
              evalue = evalue,
              meme_path = meme_path, ...)
  
}

#' @export
#' @noRd
runTomTom.default <- function(input, database = NULL,
                      outdir = "auto",
                      thresh = 10,
                      min_overlap = 5,
                      dist = "pearson",
                      evalue = TRUE,
                      meme_path = NULL, ...){
  # use TOMTOM server default values which differ from commandline defaults?
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

  command <- search_meme_path(path = meme_path, util = "tomtom")

  if (outdir == "auto") {outdir <- file.path(dirname(input$path), "tomtom")}

  database <- search_meme_database_path(path = database)

  user_flags <- prepareTomTomFlags(outdir = outdir,
                                   thresh = thresh,
                                   min_overlap = min_overlap,
                                   dist = dist,
                                   evalue = evalue, ...)
  flags <- c(user_flags, input$path, database)

  ps_out <- processx::run(command, flags, spinner = TRUE, error_on_status = FALSE)
  ps_out %>%
    process_check_error(help_fun = ~{tomtom_help(command)},
                        user_flags = cmdfun::cmd_help_parse_flags(user_flags),
                        flags_fun = ~{gsub("-", "_", .)},
                        default_help_fun = TRUE
                        )

  tomtom_out <- cmdfun::cmd_file_expect("tomtom", c("tsv", "xml", "html"), outdir = outdir)

  tomtom_results <- parseTomTom(tomtom_out$xml, query_metadata = input$metadata)

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
#' @param ... interpreted as {cmdfun}
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
  argsDict <- c("outdir" = "oc",
               "min_overlap" = "min-overlap",
               "motif_pseudo" = "motif-pseudo",
               "no_ssc" = "no-ssc",
               "incomplete_scores" = "incomplete-scores")

  flags <- cmdfun::cmd_args_all() %>%
    cmdfun::cmd_list_interp(argsDict) %>%
    cmdfun::cmd_list_to_flags()

  return(flags)
}

#' Returns tomtom help lines
#'
#' @param command path to tomtom. output of search_meme_path(util = "tomtom")
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
    dplyr::mutate(query_idx = (seq_len(nrow(.)) - 1),
                  db = as.integer(.data$db)) %>%
    dplyr::rename("db_idx" = "db",
                  "name" = "id") %>%
    # allows renaming alt column only if exists
    dplyr::rename_all(dplyr::recode, alt = "altname")
}

#' Add tomtom query lookup table info to the optional query metadata data.frame
#'
#' @param query result from get_tomtom_query_data
#' @param metadata a universalmotif_dataframe of the query motifs
#'
#' @return
#' @noRd
#'
add_query_metadata <- function(query, metadata){

  if (any(c("query_idx", "db_idx") %in% names(metadata))) {
    # these colnames are privleged use in tomtom, so reserve any original values & convert back later
    # so they don't perturb join logic
    #metadata %<>%
    #  dplyr::rename_with(~{paste0(.x, ".original")}, dplyr::matches("query_idx$|db_idx$"))
    stop("input contains query_idx or db_idx colnames")
  }

  if (!("altname" %in% names(query))) {
    # Instantiate altname column w/ NA values if not exists
    query %<>%
      dplyr::mutate(altname = NA_character_)
  }

  query %<>%
    dplyr::select("name", "altname", "db_idx", "query_idx")

  query_with_metadata <- metadata %>%
    dplyr::left_join(query, by = c("name", "altname"))
    # return user-input idx cols if any
    #dplyr::rename_with(~{gsub("\\.original", "", .x)}, dplyr::matches("query_idx.original$|db_idx.original$"))

  return(query_with_metadata)
}


#' Get match table from tomtom xml
#'
#' @param tomtom_xml_data result from xml2::read_xml(tomtom_xml_path)
#'
#' @return data.frame of match data or NULL if no matches found
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @noRd
get_tomtom_match_data <- function(tomtom_xml_data){
  matches <- xml2::xml_find_all(tomtom_xml_data, "//matches") %>%
    xml2::xml_children()

  if (length(matches) == 0){return(NULL)}

  match_df <- purrr::map(matches, xml2::xml_children) %>%
    purrr::set_names(xml2::xml_attr(matches, "idx")) %>%
    purrr::map_dfr(attrs_to_df, stringsAsFactors = FALSE, .id = "query_idx") %>%
    dplyr::mutate_at(c("pv", "ev", "qv"), as.double) %>%
    dplyr::mutate_at(c("query_idx", "idx", "off"), as.integer) %>%
    dplyr::mutate_at("rc", as.character()) %>%
    dplyr::rename("offset" = "off",
                  "pvalue" = "pv",
                  "evalue" = "ev",
                  "qvalue" = "qv",
                  "target_idx" = "idx") %>%
    dplyr::mutate(strand = ifelse(.data$rc == "y", "-", "+")) %>%
    dplyr::rename_at(c("offset", "pvalue", "evalue", "qvalue", "strand"), ~{paste0("match_", .x)}) %>%
    dplyr::select(-"rc")

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
    attrs_to_df(stringsAsFactors = FALSE) %>%
    dplyr::mutate(db_idx = (seq_len(nrow(.)) - 1)) %>%
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
#' @importFrom rlang .data
#'
#' @noRd
get_tomtom_target_data <- function(tomtom_xml_data){
  targets <- xml2::xml_find_all(tomtom_xml_data, "//targets") %>%
    xml2::xml_children()

  target_df <- targets %>%
    attrs_to_df(stringsAsFactors = FALSE) %>%
    dplyr::rename("db_idx" = "db") %>%
    dplyr::mutate_at(c("length", "nsites", "db_idx"), as.integer) %>%
    dplyr::mutate(target_idx = (seq_len(nrow(.)) - 1))

  target_pfms <- targets %>%
    purrr::map(get_probability_matrix) %>%
    purrr::map(t)

  target_df$pfm <- target_pfms

  target_data <- target_df %>%
    dplyr::mutate(match_motif = purrr::pmap(list(.data$pfm, .data$id, .data$alt, .data$nsites), ~{
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

#' Return db/target/match data as merged dataframe
#'
#' @param tomtom_xml_data result from read_xml(tomtom_xml)
#'
#' @return hits lookup table or NULL if no matches detected (will pass message)
#' @noRd
get_tomtom_hits <- function(tomtom_xml_data){

  match_data <- tomtom_xml_data %>%
    get_tomtom_match_data()

  # TODO: revise this logic
  if (is.null(match_data)) {
    message("TomTom detected no matches")
    return(NULL)
  }

  target_db_lookup <- tomtom_xml_data %>%
    get_tomtom_db_data %>%
    dplyr::select("db_idx", "db_name")

  hits <- get_tomtom_target_data(tomtom_xml_data) %>%
    dplyr::left_join(target_db_lookup, by = "db_idx") %>%
    dplyr::left_join(match_data, by = "target_idx")

  return(hits)

}

#' Merge query data w/ hits data, nesting full tomtom results
#'
#' @param query tomtom query data
#' @param hits get_tomtom_hits() output
#'
#' @return
#' @importFrom rlang sym
#' @noRd
join_tomtom_tables <- function(query, hits){
  if (is.null(hits)){
    tomtom_results <- query %>%
      dplyr::mutate(best_match_name = NA_character_,
                    best_match_altname = NA_character_,
                    best_match_offset = NA_integer_,
                    best_match_pvalue = NA_real_,
                    best_match_evalue = NA_real_,
                    best_match_qvalue = NA_real_,
                    best_match_strand = NA_character_,
                    best_match_motif = NA,
                    tomtom = NA) %>%
      dplyr::select(-"query_idx", -"db_idx")
  } else {
    tomtom_results <- query %>%
      dplyr::left_join(hits, by = c("query_idx", "db_idx")) %>%
      # Rename columns for max compatibility with universalmotif
      dplyr::rename("match_name" = "match_id",
                    "match_altname" = "match_alt") %>%
      dplyr::select(-"query_idx", -"db_idx", -"target_idx") %>% 
      dplyr::arrange(!!rlang::sym("match_qvalue"),
                     !!rlang::sym("match_pvalue"))

  }

  if (!("altname" %in% names(tomtom_results))) {
    # Instantiate altname column w/ NA values if not exists
    tomtom_results %<>%
      dplyr::mutate(altname = NA_character_)
  }

  # Nest full match data & add best_match_ columns if hits exist
  if (!is.null(hits)){
    tomtom_results %<>%
      #nest_tomtom_results()
      nest_tomtom_results_best_top_row()
  }

  return(data.frame(tomtom_results))
}

#' Helper function selecting best match by lowest evalue
#'
#' @param df
#'
#' @return
#' @noRd
#' @importFrom utils head
#'
#' @examples
#' nest_tomtom_fun(tomtom_data, tomtom_best_match_min_evalue)
tomtom_best_match_min_evalue <- function(df){
    df %>%
      dplyr::filter("match_evalue" == min("match_evalue")) %>%
      utils::head(1) %>%
      dplyr::rename_all(~{paste0("best_", .x)})
}

#' Return tomtom best match by taking top column in tomtom df
#'
#' @param df
#'
#' @return
#' @noRd
tomtom_best_match_top_row <- function(df){
  df[1,] %>%
    dplyr::rename_all(~{paste0("best_", .x)})
}

#' Nest tomtom using `fun` to return best match
#'
#' @param tomtom_results unnested tomtom results
#' @param fun function to apply to tomtom best match
#'
#' @return tomtom_results data.frame with nested `tomtom` column
#' @noRd
#'
nest_tomtom_fun <- function(tomtom_results, fun){
  nest_cols <- c("match_name",
                 "match_altname",
                 "match_motif",
                 "db_name",
                 "match_offset",
                 "match_pvalue",
                 "match_evalue",
                 "match_qvalue",
                 "match_strand")

  # Need to remove "motif" S4 column & rejoin unique entries because `tibble` or
  # `tidyr` doesn't play nice with S4. S4 is allowed in the nested column, but
  # S4 is not allowed in a parent df column... *sigh*
  tomtom_motifs <- tomtom_results %>%
    dplyr::select(dplyr::any_of(c("name", "altname", "motif"))) %>%
    unique


  tomtom_results %>%
    dplyr::select(-.data$motif) %>%
    dplyr::group_by(.data$name, .data$altname) %>%
    tidyr::nest(data = (dplyr::any_of(nest_cols))) %>%
    dplyr::mutate(best_match_info = purrr::map(.data$data, fun)) %>%
    tidyr::unnest("best_match_info") %>%
    dplyr::rename("tomtom" = "data") %>%
    dplyr::mutate("tomtom" = purrr::map(.data$tomtom, data.frame)) %>%
    # Add back motif column
    dplyr::left_join(tomtom_motifs, by = c("name", "altname")) %>%
    # Reorder columns
    dplyr::select("name", "altname",
                  dplyr::matches("[^best_|^tomtom]"),
                  dplyr::matches("motif"),
                  dplyr::contains("best_"),
                  "tomtom")
}


#' Nest tomtom results & show only best match, store all others in `tomtom` list column
#'
#' @param tomtom_results
#'
#' @return data.frame with columns w/ all data for "best" match (defined by top
#'   hit, lowest pvalue). All other matches are nested into 'tomtom' column.
#'   Which is list of data.frames for each match too the given id.
#'
#' @importFrom tidyr nest
#' @importFrom rlang .data
#'
#' @noRd
nest_tomtom_results <- function(tomtom_results){
  nest_tomtom_fun(tomtom_results, tomtom_best_match_min_evalue)
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
  nest_tomtom_fun(tomtom_results, tomtom_best_match_top_row)
}


#' Return tomtom stat information and PWMs for matched motifs
#'
#' @param tomtom_xml_path path to tomtom.xml output
#' @param use_query_data whether to use the query motif metadata stored in the .xml file in tomtom_results
#'     used only when user-supplied motif metadata doesn't exist (ie during import from external .xml)
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
#' parseTomTom("inst/extdata/tomtom.xml")
#' @noRd
parseTomTom <- function(tomtom_xml, query_metadata = NULL){
  #tomtom_xml <- "inst/extdata/dreme_example/tomtom/tomtom.xml"
  tomtom_xml_data <- xml2::read_xml(tomtom_xml)

  hits <- get_tomtom_hits(tomtom_xml_data)

  if (is.null(query_metadata)){
    # TODO:
    # CHECK DATA TYPE?

    # get query as universalmotif_df since no external metadata
    query <- tomtom_query_motif_dfs(tomtom_xml_data)
  } else {
    # get query lookup table, then join back into the metadata so query metadata
    # will be propagated during downstream joins
    query <- get_tomtom_query_data(tomtom_xml_data) %>%
      add_query_metadata(query_metadata)
  }

  tomtom_results <- join_tomtom_tables(query, hits)

  return(tomtom_results)
}

