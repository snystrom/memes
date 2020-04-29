testthat::setup(
  options(meme_db = "inst/extdata/db/fly_factor_survey_id.meme")
)
testthat::teardown(
  options(meme_db = "")
)
test_that("AME Internal functions work",{
  expect_equal(prepareAmeFlags("shuffle", "outdir"), c("--control", "--shuffle--", "--oc", "outdir"))
  expect_equal(prepareAmeFlags("test.fa", "outdir"), c("--control", "test.fa", "--oc", "outdir"))
  expect_equal(prepareAmeFlags(NULL, "outdir"), c("--oc", "outdir"))
  expect_equal(prepareAmeFlags(NA, "outdir"), c("--oc", "outdir"))
})

# expect error
runAme("inst/extdata/fasta_ex/fa1.fa", NULL)
# expect error, partitioning mode, fasta requires fasta scores
runAme("inst/extdata/fasta_ex/fa1.fa", NA)
# expect success, partitioning mode, fasta requires fasta scores
# TODO: use fasta file with scores
runAme("inst/extdata/fasta_ex/fa1.fa", NA)
# expect success, discriminative mode
runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa")
# expect error


test_that("input error checking works", {
  expect_error(
    suppressMessages(
    runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa", evalue_reportt_htreshold = 1),
               "\"evalue_report_threshold\" instead of: \" evalue_reportt_htreshold\""),
    class = "error"
  )

  expect_error(
    suppressMessages(
    runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa", evalue_reportt_htreshold = 1)
    ), class = "usethis_error"
  )
})


## Need functions to test import of all methods

expect_equal("NO", "You need to write import tests for AME")

# Different versions of ame w/ sequence import.
# maybe useful for checking import sequences function?
#ame_analysis_seq <- peaks %>%
#  resize(200, "center") %>%
#  get_sequence(dm.genome) %>%
#  runAme(evalue_report_threshold = 30, sequences = TRUE)
#
## ame fails
#ame_analysis_seq2 <- peaks %>%
#  resize(200, "center") %>%
#  get_sequence(dm.genome) %>%
#  runAme(control = fa, evalue_report_threshold = 100, sequences = T, silent = F)
## ame succeeds paritiotning mode
#ame_analysis_seq3 <- peaks %>%
#  resize(200, "center") %>%
#  get_sequence(dm.genome) %>%
#  runAme(control = NA, evalue_report_threshold = 100, sequences = T)
#
## ame discriminative fisher
#bg <- universalmotif::create_sequences()
#names(bg) <- rep("chr2L:100-200", length(bg))
#ame_analysis_seq4 <- peaks %>%
#  resize(200, "center") %>%
#  get_sequence(dm.genome) %>%
#  runAme(control = bg, evalue_report_threshold = 100, sequences = T, silent = F)
