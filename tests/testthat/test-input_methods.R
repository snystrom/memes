context("input methods & meme_database_path search_function")

#####
# MOVE THIS TO HELPER
fa <- memes:::duplicate_file("inst/extdata/fasta_ex/fa1.fa")
dreme_out <- runDreme(fa, "shuffle", e = 39)
####

test_that("database works with user-input settings", {
  fb <- "inst/extdata/db/fly_factor_survey_id.meme"
  expect_equal(basename(handle_meme_database_path(fb)), basename(fb))

  path_motif_named <- c(fb, "motif" = universalmotif::create_motif())
  expect_equal(basename(handle_meme_database_path(path_motif_named)),
               c(basename(fb), "motif"))

  noNames <- c(fb, universalmotif::create_motif())
  expect_equal(basename(handle_meme_database_path(noNames)),
               c("fly_factor_survey_id.meme", "2"))

  both_named <- c('flyfactor' = fb, "motif" = universalmotif::create_motif())
  expect_equal(basename(handle_meme_database_path(both_named)),
               c("flyfactor", "motif"))
})

test_that("dreme-res works as database",{
  # uses name
  expect_equal(handle_meme_database_path(list("test" = dreme_out)) %>% basename, "test")
  # uses order as default name
  expect_equal(handle_meme_database_path(list(dreme_out)) %>% basename, "1")
})

test_that("Type checking errors correctly", {
  expect_error(handle_meme_database_path(dreme_out), "data.frame")
  expect_error(handle_meme_database_path(2), "path must be")
  expect_error(handle_meme_database_path(factor(1,2,3)), "path must be")
})

test_that("motif input method dispatch works",{
  fb <- "inst/extdata/db/fly_factor_survey_id.meme"
  # character(1)
  expect_equal(motif_input(fb), list(metadata = NULL, path = fb))
  expect_error(motif_input(c(fb, fb)), "length == 1")

  # dreme-results
  path <- tempfile()
  expect_equal(motif_input(dreme_out, path), list(metadata = dreme_out, path = path))

  # universalmotif list
  motifList <- purrr::map(1:2, universalmotif::create_motif)
  expect_equal(motif_input(motifList, path), list(metadata = as_universalmotif_dataframe(motifList), path = path))

  # universalmotif
  motif <- universalmotif::create_motif()
  expect_equal(motif_input(motif, path), list(metadata = as_universalmotif_dataframe(motif), path = path))
})
