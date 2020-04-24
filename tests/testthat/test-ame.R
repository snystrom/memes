testthat("AME Internal functions work",{
  expect_equal(prepareAmeFlags("shuffle", "outdir"), c("--control", "--shuffle--", "--oc", "outdir"))
  expect_equal(prepareAmeFlags("test.fa", "outdir"), c("--control", "test.fa", "--oc", "outdir"))
  expect_equal(prepareAmeFlags(NULL, "outdir"), c("--oc", "outdir"))
  expect_equal(prepareAmeFlags(NA, "outdir"), c("--oc", "outdir"))
})

# expect error
runAme("inst/extdata/fasta_ex/fa1.fa", NULL)
# expect success, partitioning mode, fasta requires fasta scores
runAme("inst/extdata/fasta_ex/fa1.fa", NA)
# expect success, discriminative mode
runAme("inst/extdata/fasta_ex/fa1.fa", "inst/extdata/fasta_ex/fa2.fa")


## Need functions to test import of all methods
