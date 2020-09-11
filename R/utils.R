#' See generic get_sequence docs for details on inputs
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom Biostrings getSeq
#' @export
#' @noRd
get_sequence.GRanges <- function(regions, genome, score_column = NULL, ...){

  chrNames <- seqnames(regions)
  # NOTE: parse_genomic_coords uses 1-based coordinates, so no need to shift here.
  startPos <- start(regions)
  endPos <- end(regions)

  feature_names <- paste0(chrNames,":",startPos,"-",endPos)

  if (!is.null(score_column)){
    stopifnot(score_column %in% names(mcols(regions)))

    feature_names <- paste(feature_names, mcols(regions)[[score_column]], sep = " ")
  }


  sequences <- Biostrings::getSeq(genome, regions, ...)
  names(sequences) <- feature_names
  return(sequences)
}

#' @export
#' @noRd
get_sequence.data.frame <- function(regions, genome, score_column = NULL, ...){
  regions <- tryCatch(GenomicRanges::GRanges(regions),
           error = function(e){stop(e)})

  get_sequence.GRanges(regions = regions, genome = genome, score_column = score_column, ...)

}

#' This and all functions below are simply to support all the different GRangesList types that exist.
#' they have no special inputs that differ from the other methods.
#' @export
#' @noRd
get_sequence.GenomicRangesList <- function(regions, genome, score_column = NULL, ...){
  sequences <- lapply(regions, function(x){
    get_sequence.GRanges(regions = x, genome = genome, score_column = score_column, ...)
  })

  Biostrings::BStringSetList(sequences)
}

#' @export
#' @noRd
get_sequence.GRangesList <- function(regions, genome, score_column = NULL, ...){
  get_sequence.GenomicRangesList(regions = regions, genome = genome, score_column = score_column, ...)
}

#' @export
#' @noRd
get_sequence.CompressedGRangesList <- function(regions, genome, score_column = NULL, ...){
  get_sequence.GenomicRangesList(regions = regions, genome = genome, score_column = score_column, ...)
}

#' @export
#' @noRd
get_sequence.SimpleGRangesList <- function(regions, genome, score_column = NULL, ...){
  get_sequence.GenomicRangesList(regions = regions, genome = genome, score_column = score_column, ...)
}

#' @export
#' @noRd
get_sequence.list <- function(regions, genome, score_column = NULL, ...){
  get_sequence.GenomicRangesList(regions = regions, genome = genome, score_column = score_column, ...)
}

#' @export
#' @noRd
get_sequence.character <- function(regions, genome, score_column = NULL, ...){
  get_sequence.data.frame(regions = regions, genome = genome, score_column = score_column, ...)
}

#' Add nucleic acid sequence of regions to metadata column
#'
#' @param ranges GRanges object
#' @param genome BSgenome object or any other valid input to `Biostrings::getSeq()` (Do `showMethods(Biostrings::getSeq)` to show valid types)
#' @param name name of metadata column to hold sequence information (default: "sequence")
#'
#' @return `ranges` with new metadata column named "sequence" (or another value
#'   passed to `name`) holding the DNA or RNA sequence from `genome`
#' @export
#'
#' @examples
#' data(example_peaks, package = "memes")
#' dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3
#' add_sequence(example_peaks, dm.genome)
add_sequence <- function(ranges, genome, name = "sequence"){
  name <- rlang::enquo(name)
  seq <- ranges %>%
    get_sequence(genome)

  seq_ranges <- sequence_as_granges(seq, name = !!name)

  # Ensure seq_ranges and ranges are in the same order
  #stopifnot(identical(granges(seq_ranges), granges(plyranges::set_strand(ranges, "*"))))

  mcols(ranges)[rlang::quo_name(name)] <- mcols(seq_ranges)[rlang::quo_name(name)]
  return(ranges)
}

#' Convert Biostring w/ names of genomic coords to GRanges w/ sequence column
#'
#' Not meant to be called directly, used to convert output from get_sequence() to GRanges.
#'
#' @param seq Biostrings::BStringSet object
#' @param name name of column to hold sequence information
#'
#' @return
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#'
#' @noRd
sequence_as_granges <- function(seq, name = "sequence"){
  UseMethod("sequence_as_granges")
}

#' Macro for building sequence_as_granges
#' @importFrom rlang :=
#' @noRd
build_sequence_as_granges <- function(){
  fun <- function(seq, name = "sequence") {
    name <- rlang::enquo(name)
    seq %>%
      as.data.frame %>%
      tibble::rownames_to_column("coords") %>%
      dplyr::rename(!!name := "x") %>%
      tidyr::separate("coords", c("seqnames", "start", "end")) %>%
      GenomicRanges::GRanges(.)

  }
  return(fun)
}

#' @noRd
sequence_as_granges.DNAStringSet <- build_sequence_as_granges()
#' @noRd
sequence_as_granges.RNAStringSet <- build_sequence_as_granges()
#' @noRd
sequence_as_granges.BStringSet <- build_sequence_as_granges()

#' Write fasta file from stringset
#'
#' @param seq `Biostrings::XStringSet`
#' @param path path of fasta file to write (default: temporary file)
#'
#' @return path to created fasta file
#' @export
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' seq <- universalmotif::create_sequences()
#' \dontrun{
#' write_fasta(seq)
#' }
write_fasta <- function(seq, path = tempfile(fileext = ".fa")){
  Biostrings::writeXStringSet(x = seq,
                              filepath = path,
                              append = FALSE,
                              compress = FALSE,
                              format = "fasta")

  if (!file.exists(path)) {
    stop(paste0(path, " not created"))
  }

  return(path)
}

#' Check user's MEME install
#'
#' Checks if "meme/bin" directory exists, returns green check if yes, or red X if no. Returns the path it searched for.
#'
#' If 'meme/bin' is detected, next checks if each utility is installed, green check if yes, red X if no.
#'
#' @param meme_path path to "meme/bin" (if unset will search `MEME_BIN` environment variable or `meme_bin` option)
#'
#' @return message indicating which MEME utilities are installed and their location on disk
#' @export
#'
#' @examples
#' check_meme_install()
check_meme_install <- function(meme_path = NULL){
  cmdfun::cmd_install_check(search_meme_path, path = meme_path)
}

#' Returns logical vector indicating valid MEME-Suite install status
#'
#' Checks for a valid meme install using same heirarchy as `search_meme_path()`.
#' Returns `TRUE` if all supported utilities are found in the meme install
#' location, `FALSE` if not.
#'
#' @param path optional path to "meme/bin/"
#'
#' @return `logical(1)` indicating whether meme is installed with all supported utilities
#' @export
#'
#' @seealso search_meme_path check_meme_install
#'
#' @examples
#' # Will return TRUE if "meme/bin/" is detected
#' meme_is_installed()
#' \dontrun{
#' # will throw error if path to meme is invalid
#' meme_is_installed("bad/path")
#' }
meme_is_installed <- cmdfun::cmd_install_is_valid(search_meme_path, util = TRUE)
