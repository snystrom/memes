context("input methods & meme_database_path search_function")

#####
# MOVE THIS TO HELPER
#fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes") %>%
#  Biostrings::readDNAStringSet()
#dreme_out <- runDreme(fa, "shuffle", e = 39)
dreme_out <- system.file("extdata/dreme_example/dreme.xml", package = "memes") %>%
  importDremeXML()

fb <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes")
fb_name <- basename(fb)
####

test_that("database works with user-input settings", {
  expect_equal(basename(search_meme_database_path(fb)), basename(fb))

  path_motif_named <- c(fb, "motif" = universalmotif::create_motif())
  expect_equal(basename(search_meme_database_path(path_motif_named)),
               c(basename(fb), "motif"))

  noNames <- c(fb, universalmotif::create_motif())
  expect_equal(basename(search_meme_database_path(noNames)),
               c(fb_name, "2"))

  both_named <- c('flyfactor' = fb, "motif" = universalmotif::create_motif())
  expect_equal(basename(search_meme_database_path(both_named)),
               c("flyfactor", "motif"))
})

test_that("dreme-res works as database",{
  # uses name
  expect_equal(search_meme_database_path(list("test" = dreme_out)) %>% basename, "test")
  # uses order as default name
  expect_equal(search_meme_database_path(list(dreme_out)) %>% basename, "1")
})

test_that("Type checking errors correctly", {
  expect_error(search_meme_database_path(dreme_out), "data.frame")
  expect_error(search_meme_database_path(2), "path must be")
  expect_error(search_meme_database_path(factor(1,2,3)), "path must be")
})

test_that("motif input method dispatch works",{
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
