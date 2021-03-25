#' Remove duplicated motif entries
#'
#' This function identifies motif matrices which are duplicated in a
#' universalmotif list or universalmotif_df and removes them. This operation
#' ignores motif metadata and operates by removing all entries whose motif
#' matrices are identical. The first instance of a duplicated motif in the input
#' list is the one returned.
#'
#' @param x a universalmotif list or universalmotif_df
#'
#' @return A deduplicated list or universalmotif_df
#' @export
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' duplicated <- c(motif, motif)
#' remove_duplicate_motifs(duplicated)
remove_duplicate_motifs <- function(x){
  UseMethod("remove_duplicate_motifs")
}

#' @export
remove_duplicate_motifs.universalmotif_df <- function(x){
  remove_duplicate_motifs.data.frame(x)
}

#' @export
remove_duplicate_motifs.data.frame <- function(x){
  x %>% 
    universalmotif::to_list() %>%
    remove_duplicate_motifs.list() %>% 
    universalmotif::to_df()
}

#' @export
remove_duplicate_motifs.list <- function(x){
  ids <- identify_duplicate_motifs(x)
  x[ids$unique]
}

#' Get positional index of dups & nondup entries
#'
#' @param x a universalmotif list
#'
#' @return a list with entries "unique" for each unique entry, and "dups" for each duplicate entry
#' @noRd
identify_duplicate_motifs <- function(x){
  cor <- universalmotif::compare_motifs(x, method = "PCC")

  # Get the first occurrence of nondup motifs
  # this is done by setting diagonal & bottom triangle
  # of the Pearson matrix to 0.
  # then any off-diagonal 1's correspond to duplicate motifs.
  diag(cor) <- 0
  cor[lower.tri(cor)] <- 0

  uniqs <- which(matrixStats::colMaxs(cor) != 1)
  dups <- which(matrixStats::colMaxs(cor) == 1)
  
  return(list("unique" = uniqs,
              "dups" = dups
              ))
  
}

#' Check for duplicated motif matrices
#'
#' This function identifies whether any motif matrices in the input
#' universalmotif list or universalmotif_df are identical to each other. Note:
#' this operation is slow on large motif lists
#'
#' @param x a universalmotif list or universalmotif_df
#'
#' @return logical value indicating presence or absence of duplicated motif matrices
#' @export
#'
#' @examples
#' motif <- universalmotif::create_motif()
#' duplicated <- c(motif, motif)
#' has_duplicate_motifs(duplicated)
has_duplicate_motifs <- function(x){
  UseMethod("has_duplicate_motifs")
}

#' @export
has_duplicate_motifs.universalmotif_df <- function(x){
  has_duplicate_motifs.data.frame(x)
}

#' @export
has_duplicate_motifs.data.frame <- function(x){
  x %>%
    universalmotif::to_list() %>%
    has_duplicate_motifs.list()
}

#' @export
has_duplicate_motifs.list <- function(x){
  ids <- identify_duplicate_motifs(x)
  length(ids$unique) != length(x)
}
