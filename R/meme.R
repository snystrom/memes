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
#' @param alph one of c("dna", "rna", "protein") or path to alphabet file (default: "dna").
#' @param parse_genomic_coord `logical(1)` whether to parse genomic coordinates
#'   from fasta headers. Requires headers are in the form: "chr:start-end", or
#'   will result in an error. Automatically set to `FALSE` if `alph =
#'   "protein"`. This setting only needs to be changed if using a custom-built
#'   fasta file without genomic coordinates in the header.
#' @param combined_sites `logical(1)` whether to return combined sites information (coerces output to list) (default: FALSE)
#' @param silent Whether to suppress printing stdout to terminal (default: TRUE)
#' @param meme_path path to "meme/bin/". If unset, will use default search
#'   behavior:
#'   1. `meme_path` setting in `options()`
#'   2. `MEME_PATH` setting in `.Renviron`
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
#' For use with ChIP-seq data, the suggestion is to set: `revcomp = TRUE, minw = 5, maxw = 20`.
#' See [using MEME with ChIP-Seq Data Tips](https://groups.google.com/forum/#%21topic/meme-suite/rIbjIHbcpAE)
#' for additional details.
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
#' @details ## Licensing
#' The MEME Suite is free for non-profit use, but for-profit users should purchase a
#' license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
#'
#' @export
#'
#' @examples
#' if (meme_is_installed()) {
#' seqs <- universalmotif::create_sequences("CCRAAAW", seqnum = 4)
#' names(seqs) <- 1:length(seqs)
#' runMeme(seqs, parse_genomic_coord = FALSE)
#'
#' }
#'
runMeme <- function(input, control = NA, outdir = "auto", alph = "dna", parse_genomic_coord = TRUE,
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){
  UseMethod("runMeme")
}

#' @export
#' @noRd
runMeme.list <- function(input, control = NA, outdir = "auto", alph = "dna", parse_genomic_coord = TRUE,
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
runMeme.BStringSetList <- function(input, control = NA, outdir = "auto", alph = "dna", parse_genomic_coord = TRUE,
                            combined_sites = FALSE, silent = TRUE, meme_path = NULL, ...){
  runMeme.list(as.list(input), control, outdir, alph, combined_sites, silent, meme_path, ...)
}

#' @export
#' @noRd
runMeme.default <- function(input, control = NA, outdir = "auto", alph = "dna", parse_genomic_coord = TRUE,
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

  command <- search_meme_path(path = meme_path, util = "meme")

  ps_out <- processx::run(command, c(user_flags, input), error_on_status = FALSE, spinner = TRUE)

  ps_out %>%
    process_check_error(help_fun = ~{meme_help_flags(command)},
        user_flags = cmdfun::cmd_help_parse_flags(user_flags),
        default_help_fun = FALSE)

  print_process_stdout(ps_out, silent = silent)

  meme_out <- cmdfun::cmd_file_expect(prefix = "meme", ext = c("txt", "xml", "html"), outdir = outdir)

  importMeme(meme_out$txt, parse_genomic_coord = alph_parse_coords(alph, parse_genomic_coord), combined_sites = combined_sites)
}

#' Override parse_genomic_coords setting if alph = protein
#' @noRd
alph_parse_coords <- function(alph, parse_coords = TRUE){
  if (alph %in% c("protein")) {
    return(FALSE)
  } else {
    return(parse_coords)
  }
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
    cmdfun::cmd_list_interp()

  flagsList <- cmdfun::cmd_args_all(drop = "alph") %>%
    cmdfun::cmd_list_interp(argsDict)

  flags <- c(flagsList, alph_flags) %>%
    cmdfun::cmd_list_to_flags()

  return(flags)
}

#' Import MEME results
#'
#' @param meme_txt path to "meme.txt" output
#' @param parse_genomic_coord whether to parse sequence headers into genomic
#'   coordinates for motif position information, only works if fasta files were
#'   written such that the sequence headers are in the form: "chr:start-end", or
#'   some variation of this form (delimiters can be any of: "[^[:alnum:]]+" (ie
#'   non-alphanumeric characters)), (default = FALSE).
#' @param combined_sites whether to add `combined_sites` output which contains coordinates of each sequence, the motif sequence
#'
#' @return
#' @export
#' @importFrom tidyr nest
#' @importFrom rlang .data
#'
#' @examples
#' example_meme_txt <- system.file("extdata/meme_full.txt", package = "universalmotif")
#' importMeme(example_meme_txt)
importMeme <- function(meme_txt, parse_genomic_coord = FALSE, combined_sites = FALSE){
  meme_res <- universalmotif::read_meme(meme_txt, readsites = TRUE, readsites.meta = TRUE)

  meme_dataframe <- meme_res$motifs %>%
    as_universalmotif_dataframe() %>%
    dplyr::mutate("width" = nchar(.data$consensus))

  ##
  # Add sites info as data.frame
  if (class(meme_res$sites.meta) == "data.frame"){
    meme_res$sites.meta <- list(meme_res$sites.meta)
    names(meme_res$sites.meta) <- meme_dataframe$name
  }

  meme_sites_hits <- purrr::map(meme_res$sites.meta,
                                 #meme_dataframe$width,
                                 meme_sites_meta_to_df
                                 ) %>%
    dplyr::bind_rows(.id = "name") %>%
    dplyr::group_by(.data$name) %>%
    tidyr::nest() %>%
    dplyr::rename("sites_hits" = "data")

  meme_dataframe %<>%
    dplyr::left_join(meme_sites_hits, by = "name")

  ##
  # Coerce sites hits info to granges
  # Currently, granges nested inside a data.frame causes printing issues,
  # so I convert these back to data.frame (*sigh*)
  if (parse_genomic_coord){
    meme_dataframe %<>%
      dplyr::mutate("sites_hits" = purrr::map2(.data$sites_hits,
                                               .data$width, ~{
                                               meme_sites_meta_to_granges(.x, .y) %>%
                    # temporary until come up with a fix for printing data.frames with nested Granges
                                                 data.frame
                                               })
                    )
  }

  if (!combined_sites){
    return(meme_dataframe)
  } else {

    if (!parse_genomic_coord){
      meme_sites_combined <- meme_res$sites.meta.combined %>%
        meme_sites_meta_combined_to_df()
    } else {
      meme_sites_combined <- meme_res$sites.meta.combined %>%
        meme_sites_meta_combined_to_df() %>%
        meme_sites_meta_combined_to_granges()
    }

    results_list <- list("meme_data" = meme_dataframe,
                         "combined_sites" = meme_sites_combined)
    return(results_list)

  }
}

#' Returns MEME help lines
#'
#' @param command path to meme. output of search_meme_path(util = "meme")
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
#' @param command path to meme. output of search_meme_path(util = "meme")
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
    cmdfun::cmd_help_parse_flags(split_newline = FALSE)
}

#' @param sites the .$sites.meta output of
#'   universalmotif::read_meme(readsites = T, readsites.meta = T)
#' @noRd
meme_sites_meta_to_df <- function(sites){
  sites %>%
    as.data.frame %>%
    dplyr::rename_all(tolower)
}

#' @param sites the .$sites.meta.combined output of
#'   universalmotif::read_meme(readsites = T, readsites.meta = T)
#' @noRd
meme_sites_meta_combined_to_df <- function(sites){
  sites %>%
    as.data.frame %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename_all(~{gsub("\\.", "_", .x)})
}

#' convert site metadata into motif positions in GRanges
#'
#' @param sites_df the .$sites.meta output of
#'   universalmotif::read_meme(readsites = T, readsites.meta = T) after pasing
#'   to `meme_sites_meta_to_df()`
#' @param motif_length length of the motif (used to shift position of coordinates)
#'
#' @return
#'
#' @importFrom tidyr separate
#'
#' @noRd
meme_sites_meta_to_granges <- function(sites_df, motif_length){
  sites_df %>%
    tidyr::separate(sequence, c("seqnames", "start", "end")) %>%
    GenomicRanges::GRanges(.) %>%
    plyranges::anchor_5p() %>%
    plyranges::mutate(width = motif_length) %>%
    plyranges::shift_right(mcols(.)$position) %>%
    plyranges::select(-position)
}

#' return genomic positions of sites w/ combined pvalues & diagram information
#'
#' @param sites_df the .$sites.meta.combined output of
#'   universalmotif::read_meme(readsites = T, readsites.meta = T) after passing
#'   to `meme_sites_meta_combined_to_df()`
#' @importFrom tidyr separate
#' @noRd
meme_sites_meta_combined_to_granges <- function(sites_df){
  sites_df %>%
    tidyr::separate(sequence, c("seqnames", "start", "end")) %>%
    GenomicRanges::GRanges(.)
}
