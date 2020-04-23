###
# AME

ame_get_tfid <- function(res){
  res %>%
    dplyr::mutate(tfid = stringr::str_split_fixed(motif_ID, "_", 2) %>% .[,1])
}

ame_get_best_tfid <- function(ame){
  # get TfId with lowest p-value by group
  testthat::expect_true(c("tfid") %in% names(ame),
                        info = "run ame_get_tfid() before this function to make tfid column")
  ame %>%
    dplyr::group_by(tfid, type) %>%
    dplyr::filter(adj_p.value == min(adj_p.value)) %>%
    dplyr::ungroup()

}

ame_order_by_cluster <- function(ame){
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

  ame %>%
    dplyr::mutate(type = factor(type),
                  type_rank = as.integer(type)) %>%
    dplyr::group_by(tfid) %>%
    dplyr::mutate(nType = dplyr::n(),
                  minType = min(type_rank),
                  maxType = max(type_rank)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(nType, minType, maxType) %>%
    tibble::rowid_to_column("order")

}

ame_plot_ordered_heatmap <- function(ame){
  # TODO:
  # user-configure x-axis variable, fill-variable.
  ame %>%
    ame_get_tfid() %>%
    ame_get_best_tfid() %>%
    ame_order_by_cluster() %>%
    ggplot(aes(reorder(tfid, order), type)) +
      geom_tile(aes(fill = -log10(adj_p.value)), color = 'black', size = 0.3) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(color = "black",
                                     size = 34 / .pt)) +
      labs(x = NULL,
           y = NULL,
           fill = "-log10(p.adj)") +
      scale_fill_gradient2(low = "white",
                           high = "firebrick")

}
