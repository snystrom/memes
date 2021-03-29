test_that("Duplicate motif removal works", {

  m <- universalmotif::create_motif("CCAA")
  duplicated <- c(m,m,m)

  # universalmotif list
  expect_length(remove_duplicate_motifs(duplicated), 1)
  expect_identical(remove_duplicate_motifs(duplicated), list(m))
  # universalmotif_df
  expect_equal(nrow(remove_duplicate_motifs(universalmotif::to_df(duplicated))), 1)
  expect_identical(remove_duplicate_motifs(universalmotif::to_df(duplicated)), universalmotif::to_df(m))

  # returned motif is the first in the list
  m2 <- universalmotif::create_motif("CCAA", name = "first")
  expect_equal(remove_duplicate_motifs(c(m2, m))[[1]]["name"], "first")
  
  expect_true(has_duplicate_motifs(duplicated))
  expect_true(has_duplicate_motifs(universalmotif::to_df(duplicated)))
  expect_false(has_duplicate_motifs(c(m)))
  expect_false(has_duplicate_motifs(universalmotif::to_df(m)))

  no_dups <- c(universalmotif::create_motif("AAT"),
               universalmotif::create_motif("TTA"))

  expect_false(has_duplicate_motifs(no_dups))
  expect_false(has_duplicate_motifs(universalmotif::to_df(no_dups)))

})
