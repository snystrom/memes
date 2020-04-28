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
#'   set to FALSE. Other utilities can be used to lookup the sequence if
#'   `parse_genomic_coord` is `TRUE`.
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
#'   match. It is a good idea to use `parse_genomic_coords = TRUE`.
#' @export
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
                        user_flags = dotargs::get_help_flag_names(user_flags) %>%
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


  fimo_out <- dotargs::expected_outputs("tsv", "fimo", outdir = outdir)

  fimo_out$tsv %>%
    parseFimo() %>%
    return()
}

#' Parse flags for FIMO input
#'
#' @param bfile
#' @param parse_genomic_coord
#' @param skip_matched_sequence
#' @param max_strand
#' @param text
#' @param outdir
#' @param ...
#'
#' @return
#'
#' @noRd
prepareFimoFlags <- function(bfile, parse_genomic_coord, skip_matched_sequence, max_strand, text, outdir = outdir, ...){

  argsDict <- c("outdir" = "oc")

  flags <- dotargs::getAllArgs() %>%
    dotargs::argsToFlags(argsDict) %>%
    purrr::set_names(~{gsub("_", "-", .x)})

  if (!is.null(bfile)){
    if (!file.exists(bfile)){
      # if bfile isn't a path, user is probably inputting special keywords which
      # get wrapped in '--', but first drop all "-" in bfile input in case user
      # added '--' already.
      flags$bfile %<>% gsub("-", "", .) %>%
        gsub("(.+)", "--\\1--", .)
    }
  }

  flags %>%
    dotargs::crystallize_flags(prefix = "--")

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
