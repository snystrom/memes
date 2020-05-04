#' @noRd
sequence_input <- function(x) UseMethod("sequence_input")

#' @noRd
motif_input <- function(x, ...) UseMethod("motif_input")

#' Get sequence from GRanges
#'
#' A light wrapper around Biostrings::getSeq to return named DNAStringSets.
#'
#' @param regions GRanges, or GRangesList object. Will also accept a data.frame
#'   as long as it can be coerced to a GRanges object.
#' @param genome object of any valid type in showMethods(Biostrings::getSeq).
#'   Commonly a BSgenome object, or fasta file. Used to lookup sequences in regions.
#' @param score_column optional name of column (in mcols() of `regions`)
#'   containing a fasta score, used in AME in partitioning mode. (default: `NULL`)
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
#' get_sequence("chr2L:100-200", genomeFasta)
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
get_sequence <- function(regions, genome, score_column, ...) UseMethod("get_sequence")
