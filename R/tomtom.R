handle_tomtom_database_path <- function(path = NULL){
  f <- dotargs::build_path_handler(environment_var = "TOMTOM_DB",
                              option_name = "tomtom_db")
  f(path = path)
}


runTomTom <- function(input, database = NULL,
                      outdir = paste0(dirname(input), "/tomtom"),
                      meme_path = NULL, ...){

  command <- handle_meme_path(meme_path = meme_path, util = "tomtom")

  database <- handle_tomtom_database_path(path = database)

  flags <- prepareTomTomFlags(...)
  flags <- c(flags, input, database)

  processx::run(command, flags, spinner = T)

  tomtom_out <- dotargs::expected_outputs(c("tsv", "xml", "html"), "tomtom", outdir = outdir)

  tomtom_out %>%
    dotargs::check_files_exist()

  return(tomtom_out)
}

prepareTomTomFlags <- function(input, database, ...){
  #TODO:
  # lookup table converting arguments with - to _ so
  # user doesn't have to escape flags
  flags <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags() %>%
    dotargs::crystallize_flags()

  return(flags)
}
