skip_if(!meme_is_installed(), "MEME is not installed")

setup({

  db <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "dremeR")
  fa <- system.file("extdata/fasta_ex/fa1.fa", package = "dremeR")
  dreme_out <- runDreme(fa, "shuffle", e = 39, outdir = tempdir())

})

teardown({

  rm(db)
  rm(fa)
  rm(dreme_out)
  rm(tt_out)

})

test_that("tomtom target PWM and target metadata correctly assigned to eachother", {
  tt_out <<- runTomTom(dreme_out, database = db)
  expect_equal(tt_out$best_match_motif[[2]]@name, tt_out$best_match_name[[2]])
})

test_that("tomtom error checking suggests alternatives", {
  expect_error(
    suppressMessages(runTomTom(dreme_out, intternal = T)),
    "internal", class = "error")
  expect_error(
    suppressMessages(runTomTom(dreme_out, incomplete_score = T)),
    "incomplete_scores", class = "error")
})

# Test return vals
## best_match_motif list has names?
## motif list has names
test_that("tomtom returns correct data types", {
  motifs <- dreme_out$motif
  ## Test that returns NA cols if all motifs no match
  expect_message(nomatch <- runTomTom(motifs[[1]], database = db))
  expect_true(is.na(nomatch$tomtom))
  expect_true(is.na(nomatch$best_match_motif))

  expect_message(nomatch_multi <- runTomTom(motifs[c(1,3)], database = db))
  expect_true(all(is.na(nomatch_multi$tomtom)))
  expect_true(all(is.na(nomatch_multi$best_match_motif)))

  ## Test that returns NULL + real value if some motifs no match
  classes <- unique(vapply(tt_out$best_match_motif, class, character(1)))
  expect_true(identical(c("NULL", "universalmotif"), classes))
})

test_that("view_tomtom_hits works", {
  # correctly returns "noMatch" instead of error msg
  # Just check for list output, since it's easier...
  expect_type(view_tomtom_hits(tt_out), "list")

  # warning when no motifs match but input is >1 length
  x <- runTomTom(dreme_out$motif[c(1,3)], database = db)
  # should succeed making "no match" plots as well
  expect_type(view_tomtom_hits(x), "list")

})

# Input tests
test_that("all input types are accepted", {
  expect_success(tt_out <- runTomTom(dreme_out, database = db))
  expect_success(tt_meme_file <- runTomTom(system.file("extdata/example.meme", package = "dremeR"), database = db))
  # universalmotif
  expect_success(tt_um <- runTomTom(create_motif("CCAAAA", altname = "alt"), database = db))
  # universalmotif (no alt name)
  expect_success(tt_um_noalt <- runTomTom(create_motif("CCAAAA"), database = db))
  motifs <- c(create_motif("CCAAAA"), create_motif("TTTAAAA"))
  expect_success(runTomTom(motifs, database = db))
  # NEED TO USE MEME data also
  expect_true(FALSE)
})
