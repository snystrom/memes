handle_tomtom_database_path <- function(path = NULL){
  f <- dotargs::build_path_handler(environment_var = "TOMTOM_DB",
                              option_name = "tomtom_db")
  f(path = path)
}

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

  return(tomtom_out)
}

prepareTomTomFlags <- function(outdir, thresh, min_overlap, dist, evalue, ...){
  #TODO:
  # lookup table converting arguments with - to _ so
  # user doesn't have to escape flags
  argsDict = c("min_overlap" = "min-overlap",
               "outdir" = "oc")

  flags <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags(argsDict) %>%
    dotargs::crystallize_flags()

  return(flags)
}


readTomTom_txt <- function(txt){
  df <- readr::read_tsv(txt, col_names = c("query_id",
                                           "target_id",
                                           "optimal_offset",
                                           "pvalue",
                                           "evalue",
                                           "qvalue",
                                           "overlap",
                                           "query_consensus",
                                           "target_consensus",
                                           "orientation"),
                        skip = 1,
                        comment = "#",
                        col_types =
                        readr::cols(
                          query_id = readr::col_character(),
                          target_id = readr::col_character(),
                          optimal_offset = readr::col_double(),
                          pvalue = readr::col_double(),
                          evalue = readr::col_double(),
                          qvalue = readr::col_double(),
                          overlap = readr::col_integer(),
                          query_consensus = readr::col_character(),
                          target_consensus = readr::col_character(),
                          orientation = readr::col_character()
                        ))

  names(df) <- names(df) %>%
    gsub(" ", "_", .) %>%
    gsub("-", "", .) %>%
    tolower

  return(df)
}
