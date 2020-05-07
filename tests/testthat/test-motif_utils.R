
d <- mutate_motif(.data = dreme_out, "name" = "alt")
dreme_out %>%
  mutate_motif("name" = "altname",
               "altname" = "seq") %>%
  .$motif

motif_analysis %>%
  mutate_motif("name" = "best_match_name") %>%
  {universalmotif::view_motifs(.$motif)}

motif_analysis %>%
  dplyr::mutate(id = seq) %>%
  update_motifs() -> z
#
#
#dots <- mutate_motif(.data = dreme_out, name = id)
#mutate_motif(.data = dreme_out, name = id)
