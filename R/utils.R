#' Get sequence from GRanges
#'
#' A light wrapper around Biostrings::getSeq to return named DNAStringSets.
#'
#' @param regions GRanges object
#' @param genome object of any valid type in showMethods(Biostrings::getSeq).
#'   Commonly a BSgenome object, or fasta file. Used to lookup sequences in regions.
#' @param ... additional arguments passed to Biostrings::getSeq.
#'
#' @return Biostrings::DNAStringSet object with names corresponding to genomic coordinates
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using character string as input
#' genomeFasta <- "path/to/genome.fa"
#' get_sequence("chr2L:100-200", geonmeFasta)
#'
#' # using BSgenome object for genome
#' drosophila.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
#' get_sequence("chr2L:100-200", drosophila.genome)
#'
#' # using GRanges object for regions
#' regions <- GRanges(seqnames=Rle(c('chrX', 'chr2L', 'chr3R'), c(3, 3, 4)), IRanges(1:10, width=5))
#' get_sequence(regions, drosophila.genome)
#'
#' }
get_sequence <- function(regions, genome, ...){

  regions <- tryCatch(GenomicRanges::GRanges(regions),
                      error = function(e){return(e)}
           )

  chrNames <- seqnames(regions)
  startPos <- start(regions)
  endPos <- end(regions)

  feature_names <- paste0(chrNames,":",startPos,"-",endPos)

  sequences <- Biostrings::getSeq(genome, regions, ...)
  names(sequences) <- feature_names
  return(sequences)
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
