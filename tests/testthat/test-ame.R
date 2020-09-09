skip_if_not(meme_is_installed())

testthat::setup(
  options(meme_db = system.file("extdata/db/fly_factor_survey_id.meme", package = "memes"))
)

testthat::teardown(
  options(meme_db = "")
)

# Setup
prepare_seq <- function(n, seed){
  if (n > 26) warning("using > 26 sequences requires rewriting this function")

  seq <- universalmotif::create_sequences(seqnum = n, rng.seed = seed)
  names(seq) <- letters[1:n]
  return(seq)
}

db <- system.file("extdata/db/fly_factor_survey_id.meme", package = "memes")

n_seq <- 26
seq <- prepare_seq(n_seq, 123)

seq_score <- seq
names(seq_score) <- paste(names(seq_score), n_seq:1)

seq2 <- prepare_seq(n_seq, 321)

fa1 <- system.file("extdata/fasta_ex/fa1.fa", package = "memes")


test_that("AME Internal functions work",{
  expect_equal(prepareAmeFlags("shuffle", "outdir"), c("--control", "--shuffle--", "--oc", "outdir"))
  expect_equal(prepareAmeFlags("test.fa", "outdir"), c("--control", "test.fa", "--oc", "outdir"))
  expect_equal(prepareAmeFlags(NULL, "outdir"), c("--oc", "outdir"))
  expect_equal(prepareAmeFlags(NA, "outdir"), c("--oc", "outdir"))
})

test_that("AME Runs", {
  expect_s3_class(ame1 <- runAme(seq, NULL, evalue_report_threshold = 100), "data.frame")

  # expect error: file must have at least 2 entries
  expect_error(suppressMessages(runAme("inst/extdata/fasta_ex/fa1.fa", NULL, database = db)))

  # expect error, partitioning mode, fasta requires fasta scores
  expect_identical(runAme(seq, NA, evalue_report_threshold = 100, database = db), ame1)
  # expect success, partitioning mode, fasta requires fasta scores
  # the remaining tests may break if RNG changes, but this is easier than
  # shipping a bunch of extra test data
  # Also allows detecting if changes happen to RNG to alert users
  ame2 <- runAme(seq_score, NA, evalue_report_threshold = 100, database = db)
  expect_identical(ame2$motif_id, c("dsf_SANGER_5", "br-Z3_FlyReg"))
  # expect success, discriminative mode
  ame_discriminative <<- runAme(seq, seq2, database = db)
  expect_s3_class(ame_discriminative, "data.frame")
  expect_equal(ame_discriminative$motif_id, "Caup_Cell")
  # expect error

})

test_that("list input similar runDreme tests", {
  skip_if(TRUE, "Need to write list input tests")
  # Check list input
  expect_equal("NO", "Need to write test for list input to runAme similar to runDreme tests")
})

test_that("List input works", {
  seqList <- list("seq" = seq,
                  "seq2" = seq2)

  # use discriminative background from list
  expect_identical(runAme(seqList, "seq2", database = db)$seq, ame_discriminative)

  seqList2 <- list("A" = prepare_seq(5, 123),
                   "B" = prepare_seq(5, 321))
  # use shuffled background
  suppressMessages(ame_list_shuf <- runAme(seqList2, "shuffle", database = db))
  expect_named(ame_list_shuf)
  expect_null(unlist(ame_list_shuf))

  # Use invalid background
  expect_error(runAme(seqList2, "c", database = db), "input names: c")
})

test_that("input error checking works", {
  expect_error(
    suppressMessages(
    runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa", evalue_reportt_htreshold = 1, database = db),
               "\"evalue_report_threshold\" instead of: \" evalue_reportt_htreshold\""),
    class = "error"
  )

  expect_error(
    suppressMessages(
    runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa", evalue_reportt_htreshold = 1, database = db)
    ), class = "usethis_error"
  )
})


## Need functions to test import of all methods
test_that("ame import methods work", {
  skip_if(TRUE, "import tests not done")
  # Need to make some tests/testdata/ ame output tsv files
  expect_equal("NO", "You need to write import tests for AME")
})
