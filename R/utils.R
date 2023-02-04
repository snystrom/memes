#' See generic get_sequence docs for details on inputs
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom Biostrings getSeq
#' @export
#' @noRd
get_sequence.GRanges <- function(regions, genome, score_column = NULL, ...){

  chrNames <- seqnames(regions)
  # NOTE: parse_genomic_coord uses 1-based coordinates, so no need to shift here.
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
#' @param genome BSgenome object or any other valid input to
#'   `Biostrings::getSeq()` (Do `showMethods(Biostrings::getSeq)` to show valid
#'   types)
#' @param name name of metadata column to hold sequence information (default:
#'   "sequence"). Note, this will overwrite existing columns without warning if
#'   the name already exists.
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

#' @export
#' @noRd
sequence_as_granges.DNAStringSet <- build_sequence_as_granges()

#' @export
#' @noRd
sequence_as_granges.RNAStringSet <- build_sequence_as_granges()

#' @export
#' @noRd
sequence_as_granges.BStringSet <- build_sequence_as_granges()

#' Write fasta file from stringset
#'
#' @param seq a `Biostrings::XStringSet`
#' @param path path of fasta file to write (default: temporary file)
#'
#' @return path to created fasta file
#' @export
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' seq <- universalmotif::create_sequences()
#' \donttest{
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
#' In order to use the run* family of functions, memes must detect a local
#' install of the MEME Suite. MEME is installed in a directory named meme/bin/
#' which can be located anywhere on the filesystem, but is typically found in `~/meme/bin`.
#' If the MEME Suite is installed at `~/meme/bin`, memes can autodetect the install. However,
#' in the case that the MEME Suite is found at a nonstandard location, the user
#' may specify the location of their meme/bin in three ways:
#'
#'  1. provide the full path to `meme/bin` to the `meme_path` argument to each `run*` function.
#'  2. set the `meme_bin` option using `options(meme_bin = "path/to/meme/bin")` once per R session.
#'  3. set the `MEME_BIN` environment variable either in `.Renviron` or `~/.bashrc` with the path to `meme/bin`
#'
#' To aid the user in determining if memes can detect their `meme/bin` install,
#' `check_meme_install()` will search the aforementioned locations for a valid
#' `meme/bin`, returning green checks for each detected tool, or red X's for
#' undetected tools. Alternatively, users can run `meme_is_installed()` to get a
#' boolean value indicating whether their MEME Suite can be detected.
#'
#' `check_meme_install()` searches using the following heirarchy. This heirarchy
#' mimics how all `run*` functions search for `meme/bin`, thus the paths printed
#' from `check_meme_install()` will indicate the paths used by each `run*`
#' function. The heirarchy is as follows:
#'  1. the `meme_path` function argument if set
#'  2. the `meme_bin` option
#'  3. the `MEME_BIN` environment variable
#'  4. the default location at `~/meme/bin`
#'
#' @param meme_path path to "meme/bin" (if unset will search `MEME_BIN`
#'   environment variable or `meme_bin` option)
#'
#' @return message indicating which MEME utilities are installed and their
#'   location on disk
#'
#' @md
#'
#' @export
#'
#' @examples
#' check_meme_install()
check_meme_install <- function(meme_path = NULL){
  # TODO: temporary fix, revert to else block below after cmdfun update:
  if (is.null(meme_path)){
    x <- try(cmdfun::cmd_install_check(search_meme_path, path = meme_path), 
             silent = TRUE)
    if (is(x,"try-error")) {
      message("Cannot detect meme install")
      return(invisible(NULL))
    }
  } else {
    cmdfun::cmd_install_check(search_meme_path, path = meme_path)
  }
  
}

#' Returns logical vector indicating valid MEME-Suite install status
#'
#' Checks for a valid meme install using same heirarchy as `check_meme_install()`.
#' Returns `TRUE` if all supported utilities are found in the meme install
#' location, `FALSE` if not.
#'
#' The search heirarchy is as follows:
#'  1. the `meme_path` function argument if set
#'  2. the `meme_bin` option
#'  3. the `MEME_BIN` environment variable
#'  4. the default location at `~/meme/bin`
#'
#' @param path optional path to "meme/bin/". If unset, will follow the search
#'   heirarchy listed above.
#'
#' @return `logical(1)` indicating whether meme is installed with all supported utilities
#'
#' @md
#' @export
#'
#' @seealso [check_meme_install()]
#'
#' @examples
#' meme_is_installed()
meme_is_installed <- function(path = NULL){
  
  # temporary fix to catch if directory doesn't exist:
  # will eventually fix upstream in cmdfun
  # https://github.com/snystrom/memes/issues/65
  f <- cmdfun::cmd_install_is_valid(search_meme_path, util = TRUE)
  valid <- tryCatch(f(path = path), error = function(e) return(FALSE))
  return(valid)
}
