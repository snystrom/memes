#' Add title to complex plot with cowplot
#'
#' @param plot plot object to add title to
#' @param title text of title
#' @param ... passed to [cowplot::draw_text()]
#'
#' @return
#'
#' @noRd
#'
cowplot_title <- function(plot, title, ...){
  title <- cowplot::ggdraw() +
    cowplot::draw_text(title, ...)
  cowplot::plot_grid(plotlist = list(title, plot), ncol = 1, rel_heights = c(0.1, 1))
}

#' Compare top tomtom hits to original motif
#'
#' Although TomTom does a good job of matching unknown motifs to known motifs,
#' sometimes the top hit is not the correct assignment. It can be useful to
#' manually inspect the hits. This function provides a quick utility to compare
#' matches.
#'
#' This is intended to be a function used interactively and may not always be
#' the best tool for creating publication-quality figures. Results with matches
#' return ggseqlogo outputs which can be further manipulated using
#' [ggplot2::theme()] calls, but results containing no matches are static plots.
#'
#' @param results results data.frame from runTomTom
#' @param top_n number of matched motifs to return in plot (default: "all")
#' @param ... passed to [universalmotif::view_motifs()]
#'
#' @return plot of input motif vs the top n number of tomtom matched motifs. If
#'   no match found, will plot "No Match". Note: the "No Match" plots are not
#'   amenable to ggplot theme() manipulations, while all others are.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- runTomTom(motifs)
#' # show top 3 hits
#' view_tomtom_hits(results, top_n = 3)
#'
#' }
view_tomtom_hits <- function(results, top_n = "all", ...){
  # TODO: if tomtom is empty, return NONE as plot below??
  purrr::map2(results$motif, results$tomtom, ~{

    if (is.null(.y)) {
      return(view_tomtom_nomatch(.x, ...))
    }

    if (top_n == "all") {select <- 1:length(.y$match_motif)}
    else if (is.numeric(top_n)) {select <- 1:top_n}
    else {
      stop("n must be either 'all' or a number.")
    }

    motifList <- c(list(.x), .y$match_motif[select]) %>%
      purrr::discard(is.null)

    universalmotif::view_motifs(motifList, ...)
  })
}

#' Plot motif with "No Match" below.
#'
#' NOTE: this doesn't do any checking, and is not meant to be called directly by
#' user. Used internally in [view_tomtom_hits()]
#'
#' @param motif universalmotif
#' @param ... passed to [universalmotif::view_motifs()]
#'
#' @return
#'
#' @noRd
view_tomtom_nomatch <- function(motif, ...){
  plot <- universalmotif::view_motifs(motif, ...) %>%
    cowplot_title(motif@name, hjust = 0.25, size = 11)
  noMatch <- cowplot::ggdraw() +
    cowplot::draw_text("No Match", size = 100)

  cowplot::plot_grid(plot, noMatch, ncol = 1)
}

