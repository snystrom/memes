#' Title
#'
#' @param results
#' @param n
#'
#' @return
#' @export
#'
#' @examples
view_tomtom_hits <- function(results, n = "all"){

  purrr::map2(results$motifs, results$tomtom, ~{

    if (n == "all") {select <- 1:length(.y$match_motif)}
    else if (is.numeric(n)) {select <- 1:n}
    else {
      stop("n must be either 'all' or a number.")
    }

    motifList <- c(list(.x), .y$match_motif[select]) %>%
      purrr::discard(is.null)

    universalmotif::view_motifs(motifList)
  })
}
