#' Title
#'
#' @param input
#' @param control
#' @param outdir
#' @param database
#' @param meme_path
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
#'
#' @examples
runAme <- function(input,
       control = "shuffle",
       outdir = "auto",
       database = NULL,
       meme_path = NULL, ...){

  if (outdir == "auto") {outdir <- outdir_name(input, control)}

  # TODO: input either stringset or fasta path ~ dreme
  input <- dreme_input(input)
  # TODO: control either stringset, fasta path, or "shuffle"
  flags <- prepareAmeFlags(control, outdir, ...)
  # TODO: database can be path, or universalmotif list, (or vector: c(motifList, path)) (?)
  database <- handle_meme_database_path(database)
  command <- handle_meme_path(path = meme_path, util = "ame")

  flags <- c(flags, input$path, database)

  ps_out <- processx::run(command, flags, spinner = T, error_on_status = F)

  process_check_error(ps_out)
  # NOTE: sequences.tsv is only created when method == "fisher"
  # TODO: if `method` is unset or == "fisher", require sequences.tsv exists
  ame_out <- dotargs::expected_outputs(c("html", "tsv", "tsv"),
                                       c("ame", "ame", "sequences"),
                                       outdir = outdir) %>%
    purrr::set_names(c("html", "tsv", "sequences")) %T>%
    dotargs::check_files_exist()

}

prepareAmeFlags <- function(control, outdir, ...){

  argsDict <- c("outdir" = "oc")

  flagList <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags(argsDict) %>%
    purrr::set_names(~{gsub("_", "-", .x)})
  return(flagList)
  if (exists("control", flagList)) {
    if (flagList$control == "shuffle") {
      flagList$control <- "--shuffle--"
    }
  }

  flagList %>%
    dotargs::crystallize_flags(prefix = "--")
}

parseAme <- function(){
  # parsing this will be potentially complicated
  # might need to deploy switch statement for method since output varies depending on method type.

  # Strategey: build readr::cols() vector for each input type, the combine together using switch for import.
  # NOTE: need to test whether readr::col_* can be used in c() inside readr::cols()?

  # also maybe want flag for sequences = T if method = fisher to optionally also import sequences data,
  # relevant because data are very large and many people might not want it.
  # If you do implement this, store id in data.frame column, but provide helper function for converting to GRanges
  # if people use get_sequences, id will be the coordinates.
  # use this example in vignette.
  #
  # http://meme-suite.org/doc/ame-output-format.html

}
