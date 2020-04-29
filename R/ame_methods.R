#' Order AME results by clusters
#'
#' Reorders ame results from one or more runs for plotting heatmap showing unique/shared motifs across groups.
#'
#' @param ame ame results data.frame
#' @param id column of motif ID's to use. (default: `motif_id`)
#' @param group column to group samples by. To control order, this column must be a factor.
#'
#' @return
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo
#' @importFrom rlang !!
#'
#' @noRd
ame_order_by_cluster <- function(ame, id = motif_id, group = NULL, name = NULL){
  # orders data in "order" column first by TFs unique to each type,
  # then by motifs shared between types such that tfs are shown by unique, pairwise, 3-wise, etc.
  # starting from the first type upwards.

  # consider a factor: F with 3 levels (j), and 6 rows (i)
  # for a heatmap of the following:
  #
  # i6 | .   2   3
  # i5 | 1   .   3
  # i4 | 1   2   .
  # i3 | .   .   3
  # i2 | .   2   .
  # i1 | 1   .   .
  #    |___________
  #      j1  j2  j3
  #
  # We can compute the rank
  # by counting the number of entries for each i (ex i1 = 1, i5 = 2)
  # we can rank order by j-pairwise clusters by ranking by
  # count(F), min(F), max(F) in that order.
  # (ex i4 = c(2,1,2), while i5 = c(2,1,3))
  #
  # where F is implemented below as type_rank
  #
  # NOTE: I'm pretty sure this problem could be better solved using a binary
  # matrix & heirarchical clustering, but this works fine in my tests. If it's
  # not performant, consider revising.

  group <- enquo(group)
  id <- enquo(id)

  if (rlang::quo_is_null(group)){
    # thank you:
    # https://rpubs.com/tjmahr/quo_is_missing
    if (is.null(name)){name <- "All Regions"}

    res <- ame %>%
      dplyr::mutate(type = factor(name)) %>%
      tibble::rowid_to_column("order")

    return(res)
  }

  ame %>%
    dplyr::mutate(type = factor(!!group),
                  type_rank = as.integer(type)) %>%
    dplyr::group_by(!!id) %>%
    dplyr::mutate(nType = dplyr::n(),
                  minType = min(type_rank),
                  maxType = max(type_rank)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(nType, minType, maxType) %>%
    tibble::rowid_to_column("order")

}

#' Plot AME heatmap clustered by similarity in detected motifs
#'
#' @param ame ame results data.frame
#' @param id column of motif ids to use (default: motif_id).
#' @param group grouping column if comparing across multiple ame runs (optional, default: NULL).
#' @param value value to display as heatmap intensity. Default: -log10(adj.pvalue). Takes function or column name as input.
#' @param group_name when group = NULL, name to use for input regions. Ignored if group is set.
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo
#' @importFrom rlang !!
#'
#' @examples
#' \dontrun{
#' results <- runAme()
#' ame_plot_heatmap(results)
#' }
ame_plot_heatmap <- function(ame, id = motif_id, group = NULL, value = -log10(adj.pvalue), group_name = NULL){
  id <- enquo(id)
  group <- enquo(group)
  value <- enquo(value)

  if (rlang::quo_is_null(group)){
    res <- ame %>%
      ame_order_by_cluster(id = id, group = NULL, name = group_name)
  } else {
    res <- ame %>%
      ame_order_by_cluster(id = !!id, group = !!group)
  }

  res %>%
    ggplot(aes(reorder(!!id, order), as.factor(type))) +
      geom_tile(aes(fill = !!value), color = 'black', size = 0.3) +
      theme_bw() +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text = element_text(color = "black",
                                     size = 34 / .pt)) +
      labs(x = substitute(id),
           y = substitute(group),
           fill = substitute(value)) +
      scale_fill_gradient2(low = "white",
                           high = "firebrick")

}
