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

  # TODO: input either stringset or fasta path
  # TODO: control either stringset, fasta path, or "shuffle"

  flags <- prepareAmeFlags(input, control, outdir, hit_lo_fraction, evalue_report_threshold, ...)
  # TODO: database can be path, or universalmotif list, (or vector: c(motifList, path)) (?)
  database <- handle_meme_database_path(database)
  command <- handle_meme_path(path = meme_path, util = "ame")

  # NOTE: sequences.tsv is only created when method == "fisher"
  # TODO: if `method` is unset or == "fisher", require sequences.tsv exists
  ame_out <- dotargs::expected_outputs(c("html", "tsv", "tsv"),
                                       c("ame", "ame", "sequences"),
                                       outdir = outdir) %>%
    purrr::set_names(c("html", "tsv", "sequences")) %T>%
    dotargs::check_files_exist()

}

prepareAmeFlags <- function(input, control, outdir, hit_lo_fraction, evalue_report_threshold, ...){

  argsDict <- c("outdir" = "oc")

  flagList <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags() %>%
    purrr::set_names(~{gsub("_", "-", .x)})

  if (control == "shuffle") {flagList$control <- "--shuffle--"}

  flagList %>%
    dotargs::crystallize_flags()
}

parseAme <- function(){
  # parsing this will be potentially complicated
  # might need to deploy switch statement for method since output varies depending on method type.

  # also maybe want flag for sequences = T if method = fisher to optionally also import sequences data,
  # relevant because data are very large and many people might not want it.
  # If you do implement this, store id in data.frame column, but provide helper function for converting to GRanges
  # if people use get_sequences, id will be the coordinates.
  # use this example in vignette.
  #
  # http://meme-suite.org/doc/ame-output-format.html

}
