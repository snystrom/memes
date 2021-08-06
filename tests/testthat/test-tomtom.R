skip_if(!meme_is_installed(), "MEME is not installed")

# Setup
withr::local_options(list("meme_db" = NULL))

db <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE)
fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes", mustWork = TRUE)
dreme_out <- runDreme(fa, "shuffle", e = 39, outdir = tempdir())

fa2 <- universalmotif::create_sequences(seqnum = 2, seqlen = 100, rng.seed = 123)
names(fa2) <- seq_along(fa2)
meme_out <- runMeme(fa2, parse_genomic_coord = FALSE)
####
test_that("tomtom target PWM and target metadata correctly assigned to eachother", {
  # NOTE: don't change dist from "pearson" this is needed downstream for the NULL match class test
  suppressMessages(tt_out <<- runTomTom(dreme_out, database = db, dist = "pearson"))
  expect_equal(tt_out$best_match_motif[[2]]@name, tt_out$best_match_name[[2]])
  
  # Ensure that results are correctly sorted by descending p-value
  expect_true(all(tt_out$tomtom[[2]]$match_pval == sort(tt_out$tomtom[[2]]$match_pval)))
  # Ensure that results are correctly sorted by descending q-value
  expect_true(all(tt_out$tomtom[[2]]$match_qval == sort(tt_out$tomtom[[2]]$match_qval)))
  
})

test_that("tomtom error checking suggests alternatives", {
  expect_error(
    suppressMessages(runTomTom(dreme_out, database = db, intternal = TRUE)),
    "internal", class = "error")
  expect_error(
    suppressMessages(runTomTom(dreme_out, database = db, incomplete_score = TRUE)),
    "incomplete_scores", class = "error")
  
  # Make sure norc = TRUE imports work correctly
  expect_error(
    runTomTom(system.file("extdata/example.meme", package = "memes", mustWork = TRUE), 
              database = db, norc = TRUE), 
    NA)
})

# Test return vals
## best_match_motif list has names?
## motif list has names
test_that("tomtom returns correct data types", {
  motifs <- universalmotif::to_list(dreme_out, extrainfo = TRUE)
  ## Test that returns NA cols if all motifs no match
  ## Set impossibly low threshold to guarantee no match
  suppressMessages(expect_message(nomatch <<- runTomTom(motifs[[1]], database = db, thresh = 1e-08), "detected no matches"))
  expect_true(is.na(nomatch$tomtom))
  expect_true(is.na(nomatch$best_match_motif))

  suppressMessages(expect_message(nomatch_multi <- runTomTom(motifs[c(1,3)], database = db, thresh = 1e-08), 
                                  "detected no matches"))
  expect_true(all(is.na(nomatch_multi$tomtom)))
  expect_true(all(is.na(nomatch_multi$best_match_motif)))

  ## Test that returns NULL + real value if some motifs no match
  classes <- unique(vapply(tt_out$best_match_motif, 
                           function(x) {class(x)}, 
                           character(1)))
  expect_true(identical(c("NULL", "universalmotif"), sort(classes)))
})

test_that("view_tomtom_hits works", {
  # correctly returns "noMatch" instead of error msg
  # Just check for list output, since it's easier...
  expect_type(view_tomtom_hits(tt_out), "list")

  # warning when no motifs match but input is >1 length
  suppressMessages(x <- runTomTom(dreme_out[c(1,3),], database = db, thresh = 1e-8))
  # should succeed making "no match" plots as well
  expect_s3_class(view_tomtom_hits(x)[[1]], c("gg", "ggplot"))
  expect_s3_class(view_tomtom_hits(x)[[2]], c("gg", "ggplot"))
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
  suppressMessages(expect_message(runTomTom(motif, database = list(db1, db2)), "database size too small"))
  # expect warning re <50 motifs inaccurate p-value
  suppressMessages(expect_message(runTomTom(motif, database = list(db1, db2)), "at least 50"))
  # "too small" warning should increment as entries are added
  suppressMessages(expect_message(runTomTom(motif, database = list(db1, db2)), "(2)"))
  suppressMessages(expect_message(runTomTom(motif, database = list(db1, db2, db2)), "(3)"))
  expect_message(runTomTom(motif, database = purrr::map(1:51, ~{db1})), NA)
 
  # TODO: Fix tests below
  # expect warning re discarded motifs 
  dup_db <- system.file("extdata", "flyFactor_dups.meme", package = "memes", mustWork = TRUE)
  suppressMessages(expect_message(runTomTom(motif, database = dup_db), "duplicated in the database"))
 
  # TODO: 
  # runTomTom doesn't throw duplicate motif error
  # when using universalmotif db as input...
  #flyFactor_data %>% 
  #  dplyr::filter(consensus == "MMCACCTGYYV") %>% 
  #  to_list() %>% 
  #  {runTomTom(motif, database = ., silent = F)} -> x
  
})

test_that("Unusual db formats work", {
  skip_if(TRUE, "Need to implement tests below here")
  db1 <- universalmotif::create_motif("CRAW", name = "motif_1", altname = "1")
  db2 <- universalmotif::create_motif("CCRAAAW", name = "motif_2", altname = "2")
  motif <- universalmotif::create_motif("CCAAAAW", name = "test_motif")
  
  # db w/ no altname
  #db1["altname"] <- NA_character_
  #runTomTom(motif, database = list(db1, db1), silent = F)
  #expect_message(runTomTom(motif, database = list(db1, db1), silent = F), "duplicated in the database")
  # Goal here is warn if there is no altname? Does tomtom warn about this?
  # I think this is needed to ensure that altname exists for join operation to parse data
  # What I really need to do is fix the tomtom import functions to create an empty altname for the merge,
  # or only use altname column if it exists for the join.
  # I'll leave this failing test as a reminder.
  expect_message(runTomTom(motif, database = list(motif)), 
                 "some inputs are missing an altname. Using the motif name as altname.")
  
  # TODO:
  # warn if > 1 db with identical values Not sure what the error message will be
  # here Also not sure if this is needed. The idea would be to help users figure
  # out which db contributed the duplicated motifs. However, it looks like
  # tomtom won't catch this so I'd have to check for it before writing out the
  # data, which doesn't quite work with the current framework.
  #expect_true(FALSE)
})
