skip_if(T, "Testing plots must view manually.")
# TODO: User vdiffr
# https://github.com/r-lib/vdiffr



# check order_by_cluster
ame_analysis %>%
  dplyr::mutate(grp = c("1", "2"), grp = factor(grp, levels = c("2", "1"))) %>%
  ame_order_by_cluster(id = motif_id, group = grp)

ame_analysis %>%
  dplyr::mutate(grp = c("1", "2")) %>%
  ame_order_by_cluster(id = motif_id, group = grp)

ame_analysis %>%
  dplyr::mutate(grp = c("1", "2")) %>%
  ame_order_by_cluster(id = motif_id, group = NULL)



# check ame_plot_heatmap
## test each input type
ame_analysis %>%
  dplyr::mutate(grp = c("1", "2"), grp = factor(grp, levels = c("2", "1"))) %>%
  ame_plot_heatmap()

# y axis labeled "grp" with values 1 & 2.
ame_analysis %>%
  dplyr::mutate(grp = c("1", "2"), grp = factor(grp, levels = c("2", "1"))) %>%
  ame_plot_heatmap(id = motif_id, group = grp)

ame_analysis %>%
  dplyr::mutate(grp = c("1", "2"), grp = factor(grp, levels = c("2", "1"))) %>%
  ame_plot_heatmap(group = grp, value = -log10(evalue))

ame_analysis %>%
  ame_plot_heatmap()

# 1 y-axis value labeled 'test'
ame_analysis %>%
  ame_plot_heatmap(group_name = "test")
