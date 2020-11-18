skip_if(!meme_is_installed(), "MEME is not installed")

test_that("tomtom target PWM and target metadata correctly assigned to eachother", {
  tt_out <<- runTomTom(dreme_out, database = db)
  expect_equal(tt_out$best_match_motif[[2]]@name, tt_out$best_match_name[[2]])
  
  # Ensure that results are correctly sorted by descending p-value
  expect_true(all(tt_out$tomtom[[2]]$match_pvalue == sort(tt_out$tomtom[[2]]$match_pvalue)))
  # Ensure that results are correctly sorted by descending q-value
  expect_true(all(tt_out$tomtom[[2]]$match_qvalue == sort(tt_out$tomtom[[2]]$match_qvalue)))
  
})

test_that("tomtom error checking suggests alternatives", {
  expect_error(
    suppressMessages(runTomTom(dreme_out, database = db, intternal = TRUE)),
    "internal", class = "error")
  expect_error(
    suppressMessages(runTomTom(dreme_out, database = db, incomplete_score = TRUE)),
    "incomplete_scores", class = "error")
})

# Test return vals
## best_match_motif list has names?
## motif list has names
test_that("tomtom returns correct data types", {
  motifs <- dreme_out$motif
  ## Test that returns NA cols if all motifs no match
  expect_message(nomatch <- runTomTom(motifs[[1]], database = db))
  nomatch <- runTomTom(motifs[[1]], database = db)
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
  expect_success(runTomTom(system.file("extdata/example.meme", package = "memes", mustWork = TRUE), database = db))
  # universalmotif
  expect_success(tt_um <- runTomTom(universalmotif::create_motif("CCAAAA", altname = "alt"), database = db))
  # universalmotif (no alt name)
  expect_success(tt_um_noalt <- runTomTom(universalmotif::create_motif("CCAAAA"), database = db))
  motifs <- c(universalmotif::create_motif("CCAAAA"), universalmotif::create_motif("TTTAAAA"))
  expect_success(runTomTom(motifs, database = db))
  # runMeme output as input
  expect_success(runTomTom(meme_out, database = db))
})

test_that("tomtom works w/ nonstandard db inputs", {
  db1 <- universalmotif::create_motif("CRAW", name = "motif_1", altname = "1")
  db2 <- universalmotif::create_motif("CCRAAAW", name = "motif_2", altname = "2")
  motif <- universalmotif::create_motif("CCAAAAW", name = "test_motif")
  # expect warning that db is too small:
  expect_message(runTomTom(motif, database = list(db1, db2)), "database size too small")
  # expect warning re <50 motifs inaccurate p-value
  expect_message(runTomTom(motif, database = list(db1, db2)), "at least 50")
  # "too small" warning should increment as entries are added
  expect_message(runTomTom(motif, database = list(db1, db2)), "(2)")
  expect_message(runTomTom(motif, database = list(db1, db2, db2)), "(3)")
  expect_message(runTomTom(motif, database = purrr::map(1:51, ~{db1})), NA)
 
  # TODO: Fix tests below
  # expect warning re discarded motifs 
  # this command shows expected output: 
  # runTomTom(motif, database = "inst/extdata/db/fly_factor_survey_id.meme", silent = FALSE)
  expect_message(runTomTom(motif, database = "inst/extdata/db/fly_factor_survey_id.meme"))
  # db w/ no altname
  # > 1 db with identical values
  expect_true(FALSE)
})
