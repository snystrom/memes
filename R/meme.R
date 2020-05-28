#' Identify motifs with MEME
#'
#' NOTE: add note about the hsfrac issue
#'
#' @param input path to fasta, Biostrings::BStringSet list, or list of Biostrings::BStringSet (can generate using `get_sequence()`)
#' @param control any data type as in `input`, or a character vector of
#'   `names(input)` to use those regions as control sequences. Using sequences
#'   as background requires an alternative objective function. Users must pass a non-default value of
#'   `objfun` to `...` if using a non-NA control set (default: NA)
#' @param outdir (default: "auto")
#' @param alph one of c("dna", "rna", "protein") or path to alphabet file (default: "dna")
#' @param combined_sites `logical(1)` whether to return combined sites information (coerces output to list) (default: FALSE)
#' @param silent Whether to suppress printing stdout to terminal (default: TRUE)
#' @param meme_path path to "meme/bin/". If unset, will use default search
#'   behavior:
#'   1. `meme_path` setting in `options()`
#'   2. `MEME_PATH` settin in `.Renviron`
#' @param ...
#'
#' @return
#'
#' @details Additional arguments
#'
#' `runMeme()` accepts all valid arguments to meme as arguments passed to `...`.
#' For flags without values, pass them as `flag = TRUE`. The `dna`, `rna`, and
#' `protein` flags should instead be passed to the `alph` argument of
#' `runMeme()`.
#' The arguments passed to MEME often have many interactions with eachother, for
#' a detailed description of each argument see
#' the [MEME Help Page](meme-suite.org/doc/meme.html).
#'
#' For use with ChIP-seq data, see
#' [Using MEME with ChIP-Seq Data Tips](https://groups.google.com/forum/#%21topic/meme-suite/rIbjIHbcpAE)
#'
#' @details # Citation
#'
#' If you use `runMeme()` in your analysis, please cite:
#'
#' Timothy L. Bailey and Charles Elkan, "Fitting a mixture model by expectation
#' maximization to discover motifs in biopolymers", Proceedings of the Second
#' International Conference on Intelligent Systems for Molecular Biology, pp.
#' 28-36, AAAI Press, Menlo Park, California, 1994.
#' [pdf](https://tlbailey.bitbucket.io/papers/ismb94.pdf)
#'
#' @export
#'
#' @examples
runMeme <- function(input, control = NA, outdir = "auto", alph = "dna",
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){
  UseMethod("runMeme")
}

#' @export
#' @noRd
runMeme.list <- function(input, control = NA, outdir = "auto", alph = "dna",
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){

  x <- sequence_input_control_list(input, control)
  input <- x$input
  control <- x$control

  res <- purrr::map(input, runMeme.default,
             control = control,
             outdir = outdir,
             alph = alph,
             combined_sites = combined_sites,
             meme_path = meme_path,
             silent = silent,
             ...
             )
}

#' @export
#' @noRd
runMeme.BStringSetList <- function(input, control = NA, outdir = "auto", alph = "dna",
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){
  runMeme.list(as.list(input), control, outdir, alph, combined_sites, silent, meme_path, ...)
}

#' @export
#' @noRd
runMeme.default <- function(input, control = NA, outdir = "auto", alph = "dna",
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){

  input <- sequence_input(input)

  if (all(is.na(control))){
    control <- NA
  } else {
    control <- sequence_input(control)
  }

  if (outdir == "auto"){
    outdir <- outdir_name(input, control)
  }

  user_flags <- prepareMemeFlags(control, outdir,
                                 alph = alph, ...)

  command <- handle_meme_path(path = meme_path, util = "meme")

  ps_out <- processx::run(command, c(user_flags, input), error_on_status = FALSE, spinner = TRUE)

  ps_out %>%
    process_check_error(help_fun = ~{meme_help_flags(command)},
        user_flags = dotargs::get_help_flag_names(user_flags),
        default_help_fun = FALSE)

  print_process_stdout(ps_out, silent = silent)

  meme_out <- dotargs::expected_outputs(ext = c("txt", "xml", "html"), prefix = "meme", outdir = outdir)

  importMeme(meme_out$txt, combined_sites = combined_sites)
}

#' Return list of alphabet flag values
#'
#' @param alph one of c("dna", "rna", "protein") or a file path to alphabet file
#'
#' @return list w/ bool values of alphabet identity
#' @noRd
#'
meme_alph_to_args <- function(alph){

  flags <- list("dna" = FALSE,
       "rna" = FALSE,
       "protein" = FALSE,
       "alph" = FALSE)

  if (tolower(alph) %in% c("dna", "rna", "protein")){
    alph <- tolower(alph)

    flags[[alph]] <- TRUE

    return(flags)
  } else if (file.exists(alph)){
    flags[["alph"]] <- alph
    return(flags)
  } else {
    message(paste0(alph, " is not a valid file path on disk (file does not exist)"))
    error("alph value is invalid. Must be one of dna/rna/protein or a path to a valid file.")
  }
}

#' Convert user input flags into commandline flags for MEME
#'
#' @noRd
prepareMemeFlags <- function(control, outdir, alph, ...){
  argsDict <- c("outdir" = "oc",
                "control" = "neg")

  # handle alphabet assignment
  alph_flags <- meme_alph_to_args(alph) %>%
    dotargs::argsToFlags()

  flagsList <- dotargs::getAllArgs(drop = "alph") %>%
    dotargs::argsToFlags(argsDict)

  flags <- c(flagsList, alph_flags) %>%
    dotargs::crystallize_flags()

  return(flags)
}

#' Import MEME results
#'
#' @param meme_txt path to "meme.txt" output
#' @param parse_sequences whether to parse sequence headers into motif position
#'   information, only works if fasta files were written such that the sequence
#'   headers are in the form: "chr:start-end", or some variation of this form
#'   (delimiters can be any of: "[^[:alnum:]]+" (ie non-alphanumeric characters)).
#' @param combined_sites whether to add `combined_sites` output which contains coordinates of each sequence, the motif sequence
#'
#' @return
#' @export
#'
#' @examples
#' # If fasta headers do not have sequence information, parse_sequence must be set to FALSE
#' example_no_sequence <- system.file("extdata/meme_full.txt", package = "universalmotif", mustwork = TRUE)
#' importMeme(example_no_sequence, parse_sequences = FALSE)
#'
#' #TODO: Add example of file w/ sequence headers
importMeme <- function(meme_txt, parse_sequences = TRUE, combined_sites = FALSE){
  meme_res <- universalmotif::read_meme(meme_txt, readsites = TRUE, readsites.meta = TRUE)

  meme_dataframe <- meme_res$motifs %>%
    as_universalmotif_dataframe() %>%
    dplyr::mutate(width = nchar(consensus))

  if (parse_sequences) {
    if (class(meme_res$sites.meta) == "data.frame"){
      meme_res$sites.meta <- list(meme_res$sites.meta)
      names(meme_res$sites.meta) <- meme_dataframe$name
    }

    meme_sites_hits <- purrr::map2(meme_res$sites.meta,
                                   meme_dataframe$width,
                                   meme_sites_meta_to_granges
                                   ) %>%
      purrr::map(data.frame) %>%
      dplyr::bind_rows(.id = "name") %>%
      dplyr::group_by(name) %>%
      tidyr::nest() %>%
      dplyr::rename(sites_hits = data)

    meme_dataframe %<>%
      dplyr::left_join(meme_sites_hits, by = "name")

    if (!combined_sites){
      return(meme_dataframe)
    } else {

      meme_sites_combined <- meme_res$sites.meta.combined %>%
        meme_sites_meta_combined_to_granges()

      results_list <- list("meme_data" = meme_dataframe,
                           "combined_sites" = meme_sites_combined)
      return(results_list)

    }

  } else {
    return(meme_dataframe)
  }
}

#' Returns MEME help lines
#'
#' @param command path to meme. output of handle_meme_path(util = "meme")
#'
#' @return
#'
#' @noRd
meme_help <- function(command){
  processx::run(command, "-h", error_on_status = FALSE)$stderr
}

#' Get meme help flags as character vector
#'
#' Because meme commandline help is prefixed with [] and contains tabs, need a
#' custom implementation of flag parser.
#'
#' @param command path to meme. output of handle_meme_path(util = "meme")
#'
#' @return vector of flag arguments for meme
#' @export
#'
#' @noRd
meme_help_flags <- function(command){
  meme_help(command) %>%
    {strsplit(., "\n")[[1]]} %>%
    gsub("\t", " ", .) %>%
    gsub("\\[", "", .) %>%
    gsub("\\]", "", .) %>%
    dotargs::get_help_flag_names(processx = FALSE)
}

#' convert site metadata into motif positions in GRanges
#'
#' @param sites the .$sites.meta output of universalmotif::read_meme(readsites = T, readsites.meta = T)
#' @param motif_length length of the motif (used to shift position of coordinates)
#'
#' @return
#'
#' @noRd
meme_sites_meta_to_granges <- function(sites, motif_length){
  sites %>%
    as.data.frame %>%
    dplyr::rename_all(tolower) %>%
    tidyr::separate(sequence, c("seqnames", "start", "end")) %>%
    GenomicRanges::GRanges(.) %>%
    plyranges::anchor_5p() %>%
    plyranges::mutate(width = motif_length) %>%
    plyranges::shift_right(mcols(.)$position) %>%
    plyranges::select(-position)
}

#' return genomic positions of sites w/ combined pvalues & diagram information
#'
#' @param sites the .$sites.meta.combined output of universalmotif::read_meme(readsites = T, readsites.meta = T)
#' @noRd
meme_sites_meta_combined_to_granges <- function(sites){
  sites %>%
    as.data.frame %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename_all(~{gsub("\\.", "_", .x)}) %>%
    tidyr::separate(sequence, c("seqnames", "start", "end")) %>%
    GenomicRanges::GRanges(.)
}
