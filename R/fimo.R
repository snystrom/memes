#' Run FIMO on input sequences
#'
#' NOTE: runs best using `text = TRUE`, setting genomic coordinates in fasta
#' headers, and using `parse_genomic_coord = TRUE`.
#'
#' @param sequences path to fasta file, or stringset input.
#' @param motifs path to .meme format file, or universalmotif/universalmotif list input.
#' @param bfile path to background file, or special values: "motif", or
#'   "uniform" to use 0-order frequencies contained in the motif, or a uniform
#'   letter distribution. (default: "motif")
#' @param outdir output directory location. Only used if text = FALSE. Default:
#'   "auto" to autogenerate directory name. Note: if not using a fasta path as
#'   input, this will be a temporary location unless explicity set.
#' @param parse_genomic_coord `logical(1)` whether to parse genomic position
#'   from fasta name. Fasta names must be UCSC format positions (ie
#'   "chr:start-end"). If names of fasta entries are genomic coordinates and
#'   parse_genomic_coord == TRUE, results will contain positions of motif matches,
#'   otherwise FIMO will count positions from 1 to length of fasta file.
#' @param skip_matched_sequence `logical(1)` whether or not to include the DNA
#'   sequence of the match. Default: `FALSE`. Note: jobs will complete faster if
#'   set to FALSE. `add_sequence()` can be used to lookup the sequence after data import if
#'   `parse_genomic_coord` is `TRUE`, so setting this flag is not strictly needed.
#' @param max_strand if match is found on both strands, only report strand with
#'   best match.
#' @param text `logical(1)` (default: `TRUE`). No output files will be created on the filesystem.
#'   The results are unsorted and no q-values are computed. This setting allows
#'   fast searches on very large inputs.
#' @param meme_path path to `meme/bin/` (optional). Defaut: `NULL`, searches
#'   "MEME_PATH" environment variable or "meme_path" option.
#' @param silent `logical(1)` whether to suppress stdout/stderr printing to
#'   terminal (default: TRUE). NOTE: if `text = TRUE`, setting `silent = TRUE`
#'   will print all FIMO matches to terminal if the run is successful. Can be
#'   useful for troubleshooting.
#' @param ... additional commandline arguments to fimo see:
#'   \code{\link{http://meme-suite.org/doc/fimo.html?man_type=web}}
#'
#' @return GenomicRanges object contining positions of each match. Note: if
#'   `parse_genomic_coords = FALSE`, each `seqnames` entry will be the fasta
#'   header, and start/end will be the position within that sequence of the
#'   match. It is a good idea to use `parse_genomic_coords = TRUE` if using
#'   `get_sequence()` to create input sequences.
#' @export
#'
#' @details Additional arguments passed to `...`
#'
#'
#' | FIMO Flag         | allowed values | default | description                |
#' |:-----------------:|:--------------:|:-------:|:---------------------------|
#' | alpha             | `numeric(1)`   | 1       | alpha for calculating position-specific priors. Represents fraction of sites that are binding sites of TF of interest. Used in conjunction with `psp` |
#' | bfile             | "motif", "motif-file", "uniform", file path, | "motif" | If "motif" or "motif-file", use 0-order letter frequencies from motif. "uniform" sets uniform letter frequencies. |
#' | max_stored_scores | `integer(1)`   | NULL    | maximum number of scores to be stored for computing q-values. used when `text = FALSE`, see FIMO webpage for details |
#' | motif_pseudo      | `numeric(1)`   | 0.1     | pseudocount added to motif matrix |
#' | no_qvalue         | `logical(1)`   | FALSE   | only needed when `text = FALSE`, do not compute q-value for each p-value |
#' | norc              | `logical(1)`   | FALSE   | Do not score reverse complement strand |
#' | prior_dist        | file path      | NULL    | file containing binned distribution of priors |
#' | psp               | file path      | NULL    | file containing position specific priors. Requires `prior_dist` |
#' | qv_thresh         | `logical(1)`   | FALSE   | use q-values for the output threshold |
#' | thresh            | `numeric(1)`   | `1e-4`  | output threshold for returning a match, only matches with values less than `thresh` are returned. |
#'
#'
#' @details # Citation
#' If you use `runFimo()` in your analysis, please cite:
#'
#' Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, "FIMO:
#' Scanning for occurrences of a given motif", Bioinformatics, 27(7):1017-1018,
#' 2011. [full text](http://bioinformatics.oxfordjournals.org/content/early/2011/02/16/bioinformatics.btr064.full)
#'
#' @details ## Licensing
#' The MEME Suite is free for non-profit use, but for-profit users should purchase a
#' license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
#'
#' @md
#'
#' @examples
#' \dontrun{
#' seq <- universalmotif::create_sequences()
#' # sequences must have names in their fasta headers
#' names(seq) <- seq_along(seq)
#' motif <- universalmotif::create_motif()
#' runFimo(seq, motif, parse_genomic_coord = FALSE)
#' }
runFimo <- function(sequences, motifs, bfile = "motif",
                    outdir = "auto",
                    parse_genomic_coord = TRUE,
                    skip_matched_sequence = TRUE,
                    max_strand = TRUE,
                    text = TRUE,
                    meme_path = NULL,
                    silent = TRUE,
                    ...){

  sequences <- sequence_input(sequences)
  motifs <- motif_input(motifs)

  if (outdir == "auto") outdir <- paste0(outdir_name(sequences, motifs$path), "_fimo")

  user_flags <- prepareFimoFlags(bfile = bfile,
                            parse_genomic_coord = parse_genomic_coord,
                            skip_matched_sequence = skip_matched_sequence,
                            max_strand = max_strand,
                            text = text,
                            outdir = outdir,
                            ...)

  flags <- c(user_flags, motifs$path, sequences)

  command <- handle_meme_path(path = meme_path, util = "fimo")

  ps_out <- processx::run(command, flags, error_on_status = FALSE)

  ps_out %>%
    process_check_error(help_fun = ~{fimoHelp(command)},
                        user_flags = cmdfun::cmd_help_parse_flags(user_flags) %>%
                          # filter out special inputs to bfile
                          grep("--$", ., invert = TRUE, value = TRUE),
                        flags_fun = ~{gsub("-", "_", .x)}
                        )

  print_process_stdout(ps_out, silent = silent)
  print_process_stderr(ps_out, silent = silent)

  if (!is.na(text)) {
    if (text){
      if (ps_out$stdout == "") {
        message("No matches were detected")
        return(NULL)
      }
      fimo_res <- ps_out$stdout %>%
        parseFimo()
      return(fimo_res)
    }
  }


  fimo_out <- cmdfun::cmd_file_combn("fimo", "tsv", outdir = outdir)

  fimo_out$tsv %>%
    parseFimo() %>%
    return()
}

#' Parse flags for FIMO input
#'
#' @param bfile background file
#' @param parse_genomic_coord parse genomic coordinates from name
#' @param skip_matched_sequence skip matched
#' @param max_strand max strand
#' @param text text
#' @param outdir outdir
#' @param ... ...
#'
#' @return
#'
#' @noRd
prepareFimoFlags <- function(bfile, parse_genomic_coord, skip_matched_sequence, max_strand, text, outdir = outdir, ...){

  argsDict <- c("outdir" = "oc")

  flags <- cmdfun::cmd_args_all() %>%
    cmdfun::cmd_args_to_flags(argsDict) %>%
    purrr::set_names(~{gsub("_", "-", .x)})

  if (!is.null(bfile)){
    if (!file.exists(bfile) & bfile %in% c("motif", "uniform")){
      # if bfile isn't a path, user is probably inputting special keywords which
      # get wrapped in '--', but first drop all "-" in bfile input in case user
      # added '--' already.
      flags$bfile %<>% gsub("-", "", .) %>%
        gsub("(.+)", "--\\1--", .)
    }
  }

  flags %>%
    cmdfun::cmd_list_crystallize(prefix = "--")

}

#' Returns help string for FIMO
#'
#' @param command
#'
#' @return
#'
#' @noRd
fimoHelp <- function(command){
  processx::run(command, error_on_status = FALSE)$stderr
}

#' Import fimo matches as GRanges object.
#'
#' @param fimo_tsv path to fimo.tsv output
#'
#' @return GenomicRanges object for each match position. Note: if
#'   parse_genomic_coords == FALSE, each `seqnames` entry will be the fasta
#'   header, and start/end will be the position within that sequence of the
#'   match.
#'
#' @noRd
parseFimo <- function(fimo_tsv){

  fimo_matches <- readr::read_tsv(fimo_tsv, comment = "#") %>%
    dplyr::rename_all(~{gsub("-", "", .)}) %>%
    dplyr::rename("seqnames" = "sequence_name") %>%
    GenomicRanges::GRanges()
    # Compute q-value?
    #dplyr::group_by(motif_alt_id) %>%
    #dplyr::mutate(q.value = p.value/n())
}

#' Import FIMO results
#'
#' @param fimo_tsv
#'
#' @return
#' @export
#'
#' @examples
#'
#' \dontrun{
#' importFimo("fimo.tsv")
#' }
importFimo <- function(fimo_tsv){
  parseFimo(fimo_tsv)
}
