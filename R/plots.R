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
view_tomtom_hits <- function(results, top_n = "all"){
  # TODO: if tomtom is empty, return NONE as plot below??
  purrr::map2(results$motif, results$tomtom, ~{

    if (is.null(.y)) {
      return(view_tomtom_nomatch(.x))
    }

    # Needed to handle when tomtom discovers no hits for any motifs, in which case
    # tomtom is NA instead of NULL so the column is kept in the dataframe
    if (all(is.na(.y))) {
      return(view_tomtom_nomatch(.x))
    }

    if (top_n == "all") {select <- seq_len(length(.y$match_motif))}
    else if (is.numeric(top_n)) {select <- seq_len(top_n)}
    else {
      stop("n must be either 'all' or a number.")
    }

    motifList <- c(list(.x), .y$match_motif[select]) %>%
      purrr::discard(is.null)

    universalmotif::view_motifs(motifList)
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
#' @importFrom ggplot2 ggtitle theme element_text
#'
#' @noRd
view_tomtom_nomatch <- function(motif){
  # Thanks, Hadley: http://r-pkgs.had.co.nz/description.html
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  motif_logo <- universalmotif::view_motifs(motif)

  motif_logo <- motif_logo +
    ggtitle(motif@name) +
    theme(plot.title = element_text(hjust = 0.5))

  nomatch_logo <- nomatch_logo()

  cowplot::plot_grid(motif_logo, nomatch_logo, ncol = 1, rel_heights = c(1,0.6))
}


#' Returns motif matrix of "NO MATCH"
#'
#' To use ggseqlogo to render the "NO MATCH" text,
#' I need to build a matrix with custom alphabet
#'
#' In matrix form NO MATCH is a matrix with diagonal all 1,
#' except at "space" position where all values are 1
#'
#' @return
#'
#' @noRd
nomatch_matrix <- function(){
  m <- matrix(0,
              nrow = 7,
              ncol = 7)
  diag(m) <- 1
  n <- c("N", "O", "M", "A", "T", "C", "H")
  rownames(m) <- c("N", "O", "M", "A", "T", "C", "H")
  no <- m[,c(1,2)]
  match <- m[,c(3:7)]
  space <- matrix(1, nrow = 7, ncol = 1)
  mat <- cbind(no, space, match)

  return(mat)
}

#' Returns "NO MATCH" ggseqlogo
#'
#' @return

#' @importFrom ggplot2 element_text
#' @importFrom ggseqlogo make_col_scheme ggseqlogo
#' @noRd
nomatch_logo <- function(){
  mat <- nomatch_matrix()
  alph <- rownames(mat)
  col <- ggseqlogo::make_col_scheme(chars = alph, cols = rep("#333333", length(alph)))

  ggseqlogo::ggseqlogo(mat,
                       namespace = alph,
                       method = "bits",
                       col_scheme = col) +
    ggplot2::theme(axis.text = element_text(color = "white"),
                   axis.text.x = element_text(color = "white"),
                   axis.text.y = element_text(color = "white"),
                   axis.title = element_text(color = "white"))
}
