
test_that("motif_utils work", {
  seq1 <- "CCRAAAW"
  seq2 <- "TTAAGGC"
  name1 <- "seq1"
  name2 <- "seq2"

  motifs <- list(universalmotif::create_motif(seq1, name = name1),
                 universalmotif::create_motif(seq2, name = name2))

  df <- as_universalmotif_dataframe(motifs)

  expect_true(is_universalmotif_dataframe(df))

  df %<>%
    dplyr::mutate(altname = name)

  expect_equal(update_motifs(df)$motif[[1]]@altname, name1)

  df %<>%
    dplyr::mutate(strand = c("-", "+"))

  expect_equal(update_motifs(df)$motif[[1]]@strand, "-")
  expect_equal(update_motifs(df)$motif[[2]]@strand, "+")

  # The following changes to df_bad should be ignored
  df_bad <- df %>%
    dplyr::mutate(consensus = c("A", "B"))

  expect_equal(update_motifs(df_bad)$motif[[1]]@consensus, seq1)
  expect_equal(update_motifs(df_bad)$motif[[2]]@consensus, seq2)

  df_bad %<>%
    dplyr::mutate(icscore = -1)

  expect_equal(update_motifs(df_bad)$motif[[1]]@icscore, motifs[[1]]@icscore)
  expect_equal(update_motifs(df_bad)$motif[[2]]@icscore, motifs[[2]]@icscore)

  df_bad %<>%
    dplyr::mutate(alphabet = -1)

  expect_equal(update_motifs(df_bad)$motif[[1]]@alphabet, motifs[[1]]@alphabet)
  expect_equal(update_motifs(df_bad)$motif[[2]]@alphabet, motifs[[2]]@alphabet)

  # Error if motif column cannot be found
  df_bad %<>%
    dplyr::rename(not_motif = motif)
  expect_error(update_motifs(df_bad))

  # error if motif column isn't universalmotif list
  df_bad <- df %>%
    dplyr::rename("name" = "motif",
                  "motif" = "name")

  expect_error(update_motifs(df_bad))

})
