#' Motif enrichment using AME
#'
#' @param input path to fasta, or DNAstringset (optional: DNAStringSet object
#'   names contain fasta score, required for partitioning mode)
#' @param control default: "shuffle", or set to
#'   DNAstringset or path to fasta file to use those sequences in discriminative
#'   mode. Set to `NA` for partitioning based on input fasta score (see
#'   get_sequence for assigning fasta score).
#' @param outdir default: auto
#' @param method default: fisher
#' @param database path to .meme format file, universalmotif list object, dreme
#'   results data.frame, or list() of multiple of these. If objects are assigned names in the list,
#'   that name will be used as the database id. It is highly recommended you set
#'   a name if not using a file path so the database name will be informative,
#'   otherwise the position in the list will be used as the database id. For
#'   example, if the input is: list("motifs.meme", list_of_motifs), the database
#'   id's will be: "motifs.meme" and "2". If the input is list("motifs.meme",
#'   "customMotifs" = list_of_motifs), the database id's will be "motifs.meme"
#'   and "customMotifs".
#' @param meme_path
#' @param sequences `logical(1)` add results from `sequences.tsv` to `sequences`
#'   list column to returned data.frame. Valid only if method = "fisher". See
#'   [AME outputs](http://alternate.meme-suite.org/doc/ame-output-format.html)
#'   webpage for more information.
#' @param silent whether to suppress stdout (default: TRUE), useful for debugging.
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
#'
#' @examples
#' \dontrun{
#' dreme_out <- runDreme("input.fa")
#' runAme("input.fa", database = "jasparmotifs.meme")
#' runAme("input.fa", database = list("jasparmotifs.meme", "my_dreme_motifs" = dreme_out))
#' }
runAme <- function(input,
       control = "shuffle",
       outdir = "auto",
       method = "fisher",
       database = NULL,
       meme_path = NULL,
       sequences = FALSE, silent = TRUE, ...){

  input <- sequence_input(input)

  if (!is.na(control)){
    control <- sequence_input(control)
  }

  # Autodetect outdir path from path names
  # note: this line must run after input&control are parsed to paths
  if (outdir == "auto") {outdir <- outdir_name(input, control)}

  user_flags <- prepareAmeFlags(control, outdir, method, ...)
  database <- handle_meme_database_path(database)
  command <- handle_meme_path(path = meme_path, util = "ame")

  flags <- c(user_flags, input, database)

  ps_out <- processx::run(command, flags, spinner = T, error_on_status = F)

  #uf <- dotargs::get_help_flag_names(user_flags) %>%
  #  grep("shuffle", ., invert = TRUE, value = TRUE)
  #return(uf)
  #return(user_flags)
  ps_out %>%
    process_check_error(help_fun = ~{ame_help(command)},
                        user_flags = dotargs::get_help_flag_names(user_flags) %>%
                          grep("shuffle", ., invert = TRUE, value = TRUE),
                        flags_fun = ~{gsub("-", "_", .x)})

  print_process_stdout(ps_out, silent = silent)
  print_process_stderr(ps_out, silent = silent)

  # NOTE: sequences.tsv is only created when method == "fisher"
  ame_out <- dotargs::expected_outputs(c("tsv", "html"), "ame", outdir)
  if (method == "fisher"){
    ame_seq <- dotargs::expected_outputs("tsv", "sequences", outdir)
    ame_out$sequences <- ame_seq[[1]]
  }

  ame_out %>%
    dotargs::check_files_exist()

  import_sequences <- FALSE
  if (method == "fisher" & sequences == TRUE){
    import_sequences <- ame_out$sequences
  }

  importAme(path = ame_out$tsv, method = method, sequences = import_sequences)
}

#' Returns ame help lines
#'
#' @param command path to ame. output of handle_meme_path(util = "ame")
#'
#' @return
#'
#' @noRd
ame_help <- function(command){
  processx::run(command, "--help", error_on_status = FALSE)$stderr
}

prepareAmeFlags <- function(control, outdir, method, ...){

  argsDict <- c("outdir" = "oc")

  flagList <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags(argsDict) %>%
    purrr::set_names(~{gsub("_", "-", .x)})

  if (exists("control", flagList)) {
    if (flagList$control == "shuffle") {
      flagList$control <- "--shuffle--"
    }
  }

  flagList %>%
    dotargs::crystallize_flags(prefix = "--")
}

#' Parse ame output
#'
#' @param path path to ame results file ("ame.tsv")
#' @param method ame run method used (one of: c("fisher", "ranksum", "dmhg3",
#'   "dmhg4", "pearson", "spearman")). Default: "fisher".
#' @param sequences FALSE or path to sequences file (only valid for method = "fisher")
#'
#' @return data.frame with method-specific results. See [AME
#'   results](http://meme-suite.org/doc/ame-output-format.html) webpage for more
#'   information. If sequences is set to a path to the sequences.tsv and method
#'   = "fisher", the list-column `sequences` will be appended to resulting
#'   data.frame.
#'
#' @importFrom magrittr %<>%
#' @export
#'
#' @family import
#'
#' @examples
importAme <- function(path, method = "fisher", sequences = FALSE){

  cols <- get_ame_coltypes(method)

  data <- readr::read_tsv(path,
                          col_names = names(cols$cols),
                          col_types = cols,
                          skip = 1,
                          comment = "#"
                          )

  if (nrow(data) == 0){
    message("AME detected no enrichment")
    return(NULL)
  }

  if (sequences == FALSE | method != "fisher"){return(data)}

  if (sequences != FALSE & method == "fisher"){
    seq <- importAmeSequences(sequences)

    if (is.null(seq)){
      return(data)
    }

    seq %<>%
      dplyr::group_by(motif_id, motif_db) %>%
      tidyr::nest() %>%
      dplyr::rename("sequences" = "data") %>%
      data.frame

    return(dplyr::left_join(data, seq, by = c("motif_id","motif_db")))
  }

}

#' Helper for combining readr::cols() objects
#'
#' @param col readr::cols() object
#' @param cols_list list of readr::cols() objects. **NOTE** order MATTERS!
#'
#' @return combined cols() object of all inputs
#' @noRd
combine_cols <- function(col, cols_list){
    # original cols to col
    # pass extra cols to cols as list

    out <- col

    purrr::walk(cols_list, ~{
      out$cols <<- c(out$cols, .x$cols)
    })
    return(out)
}

#' Import AME sequences information for method="fisher" runs.
#'
#' @param path path to sequences.tsv ame output file
#'
#' @return data.frame with columns:
#'  - motif_db: name of motif db the identified motif was found in
#'  - motif_id: name of identified motif (primary identifier)
#'  - seq_id (name of fasta entry)
#'  - label_[fasta|pwm]_score: score used for labeling positive (either fasta or pwm score)
#'  - class_[fasta|pwm]_score: score used for classifying positives (either fasta or pwm score)
#'  - class: whether a sequence was called true-positve ("tp") or false-positive ("fp")
#'
#' @importFrom magrittr %<>%
#'
#' @examples
#' \dontrun{
#' importAmeSequences("path/to/ame/sequences.tsv")
#' }
#'
#' @noRd
importAmeSequences <- function(path){

  sequences <- readr::read_tsv(path,
                               col_types = readr::cols("c", "c", "c", "d", "d", "c"),
                               col_names = T,
                               comment = "#")
  if (nrow(sequences) == 0){
    message("Sequences output is empty")
    return(NULL)
  }

  sequences %<>%
    dplyr::rename_all(tolower) %>%
    # positions 4 & 5 encode which score was used to "label" (4) vs "classify" (5)
    # can be either PWM for Fasta score, and can't predict which one easily, so
    # just prefix these two.
    dplyr::rename_at(4, function(x){paste0("label_", x)}) %>%
    dplyr::rename_at(5, function(x){paste0("class_", x)})

}

#' Generate columntypes/names for ame results.
#'
#' @param method ame run method used (one of: c("fisher", "ranksum", "dmhg3",
#'   "dmhg4", "pearson", "spearman")).
#'
#' @return readr::cols object w/ names & datatypes for given method
#'
#' @noRd
get_ame_coltypes <- function(method){
  # Strategey: build readr::cols() vector for each input type, the combine together using switch for import.
  # NOTE: need to test whether readr::col_* can be used in c() inside readr::cols()?

  cols_common <- readr::cols("rank" = "i",
                             "motif_db" = "c",
                             "motif_id" = "c",
                             "motif_alt_id" = "c",
                             "consensus" = "c",
                             "pvalue" = "d",
                             "adj.pvalue" = "d",
                             "evalue" = "d",
                             "tests" = "i"
                             )

  cols_fisher_ranksum_dmhg <- readr::cols("fasta_max" = "d",
                                          "pos" = "i",
                                          "neg" = "i"
                                          )
  cols_fisher <- readr::cols("pwm_min" = "d",
                             "tp" = "i",
                             "tp_percent" = "d",
                             "fp" = "i",
                             "fp_percent" = "d"
                             )

  cols_ranksum <- readr::cols("u" = "d",
                              "pleft" = "d",
                              "pright" = "d",
                              "pboth" = "d",
                              "adj.pleft" = "d",
                              "adj.pright" = "d",
                              "adj.both" = "d"
                              )

  cols_pearson <- readr::cols("pearson_cc" = "d",
                              "mean_squared_error" = "d",
                              "slope" = "d",
                              "intercept" = "d"
                              )

  cols_spearman <- readr::cols("spearman_cc" = "d")

  method <- gsub("[3,4]dmhg", "dmhg", method)
  cols <- switch(method,
         fisher = combine_cols(cols_common, list(cols_fisher_ranksum_dmhg, cols_fisher)),
         ranksum = combine_cols(cols_common, list(cols_fisher_ranksum_dmhg, cols_ranksum)),
         dmhg = combine_cols(cols_common, list(cols_fisher_ranksum_dmhg)),
         pearson = combine_cols(cols_common, list(cols_pearson)),
         spearman = combine_cols(cols_common, list(cols_spearman)),
         stop(paste0(method, " is not a valid method"))
         )

  return(cols)

}
