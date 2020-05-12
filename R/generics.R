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


#' Denovo motif discovery of target regions
#'
#' @param input regions to scan for motifs. Can be any of:
#'   - path to fasta file
#'   - DNAStringSet object (can be generated from GRanges using `get_sequence()`)
#'   - List of DNAStringSet objects (generated from `get_sequence()`)
#'   - *NOTE:* If you want to retain the raw dreme output files, you must use a
#'   path to fasta file as input, or specify an "outdir"
#' @param control regions to use as background for motif search. Can be any of:
#'   - path to fasta file
#'   - DNAStringSet object (can be generated from GRanges using get_sequence)
#'   - if `input` is a list of DNAStringSet objects, a character vector of names
#'   in `input` will use those sequences as background. runDreme will not scan
#'   those regions as input.
#'   - "shuffle" to use dreme's built-in dinucleotide shuffle feature (NOTE: if
#'   `input` is a list object with an entry named "shuffle", the list entry will
#'   be used instead).
#'   Optionally can also pass `s = <any number>` to `...` to use as the random
#'   seed during shuffling. If no seed is passed, dreme will use 1 as the random
#'   seed, so results will be reproducible if rerunning. **NOTE:** beware
#'   system-specific differences. As of v5, dreme will compile using the default
#'   python installation on a system (either python2.7 or python3). The random
#'   number generator changed between python2.7 and python3, so results will not
#'   be reproducible between systems using python2.7 vs 3.
#' @param outdir path to output directory of dreme files, or "auto" to autogenerate path. Default: location of
#'   input fasta in dir named "\<input\>_vs_\<control\>". If input is
#'   DNAstringset, will be temporary path. This means that if you want to save
#'   the raw output files, you must use fasta files as input or use an
#'   informative (and unique) outdir name. dremeR will **not check** if it
#'   overwrites files in a directory. Directories will be recursively created if needed.
#' @param meme_path optional, path to "meme/bin/" on your local machine.
#'   runDreme will search 3 places in order for meme if this flag is unset:
#'    1. the option "meme_bin" (set with options(meme_bin = "path/to/meme/bin"))
#'    2. the environment variable "MEME_PATH" (set in .Renviron)
#'    3. "~/meme/bin/" as the default location
#'    - If the user sets meme_path in the function call, this value will always be used
#'
#' @param silent whether to suppress printing dreme stdout as a message when
#'   finishing with no errors. Can be useful for troubleshooting in situations
#'   where no motifs are discovered, but command completes successfully.
#'   (default: TRUE)
#'
#' @param ... dreme flags can be passed as R function arguments to use
#'   non-default behavior. For a full list of valid arguments, run your local
#'   install of dreme -h, or visit the dreme documentation
#'   [website](http://meme-suite.org/doc/dreme.html?man_type=web).
#'
#' @details
#' In addition to allowing any valid flag of dreme to be passed to `...`, we
#' provide a few user-friendly aliases for common flags which will hopefully
#' make it things easier to remember. For example, m = 2 will search for 2
#' motifs. This is equivalent to setting nmotifs = 2.
#'
#' List of aliased values which can be passed to `...`
#'  - nmotifs = max number of motifs to search for
#'  - sec = max runtime in seconds
#'  - evalue = max evalue cutoff
#'  - seed = random seed if using "shuffle" as control
#'  - ngen = number of REs to generalize
#'
#' @return data.frame with statistics for each discovered motif. The `motif`
#'   column contains a universalmotif object representation in PCM format of
#'   each DREME motif. If no motifs are discovered, returns NULL.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @md
#'
#' @examples
#' \dontrun{
#' # Runs dreme with default settings, shuffles input as background
#' runDreme("input.fa", "shuffle")
#'
#' # Runs searching for max 2 motifs, e-value cutoff = 0.1, explicitly using the DNA alphabet
#' runDreme("input.fa", "shuffle", m = 2, e = 0.1, dna = T)
#' }
runDreme <- function(input, control, outdir = "auto", meme_path = NULL, silent = TRUE, ...) {
  UseMethod("runDreme")
}
