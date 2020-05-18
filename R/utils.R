#' See generic get_sequence docs for details on inputs
#' @importFrom GenomicRanges seqnames start end mcols
#' @export
#' @noRd
get_sequence.GRanges <- function(regions, genome, score_column = NULL, ...){

  chrNames <- seqnames(regions)
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

  return(sequences)
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

#' Write fasta file from stringset
#'
#' @param seq DNAstringset
#' @param path path of fasta file to write (default: temporary file)
#'
#' @return path to created fasta file
#' @export
#'
#' @examples
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
  dotargs::check_install(handle_meme_path, path = meme_path)
}
