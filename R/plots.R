#' Compare top tomtom hits to original motif
#'
#' Although TomTom does a good job of matching unknown motifs to known motifs,
#' sometimes the top hit is not the correct assignment. It can be useful to
#' manually inspect the hits. This function provides a quick utility to compare
#' matches.
#'
#' @param results results data.frame from runTomTom
#' @param top_n number of matched motifs to return in plot (default: "all")
#'
#' @return plot of input motif vs the top n number of tomtom matched motifs
#' @export
#'
#' @examples
#' \dontrun{
#' results <- runTomTom(motifs)
#' # show top 3 hits
#' view_tomtom_hits(results, top_n = 3)
#'
#' }
view_tomtom_hits <- function(results, top_n = "all"){

  purrr::map2(results$motifs, results$tomtom, ~{

    if (top_n == "all") {select <- 1:length(.y$match_motif)}
    else if (is.numeric(top_n)) {select <- 1:top_n}
    else {
      stop("n must be either 'all' or a number.")
    }

    motifList <- c(list(.x), .y$match_motif[select]) %>%
      purrr::discard(is.null)

    universalmotif::view_motifs(motifList)
  })
}
