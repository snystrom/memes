context("XML import works")

test_that("DREME XML import",{
  expect_s3_class(importDremeXML("inst/extdata/dreme_example/dreme.xml"), "data.frame")
  expect_true(is_dreme_results(importDremeXML("inst/extdata/dreme_example/dreme.xml")))
})

# testthat NULL tomtom result returns empty columns
test_that("TomTom Import works",{
  expect_s3_class(importTomTomXML("inst/extdata/tomtom_ex/tomtom_good.xml"), "data.frame")
  # Should return query motifs w/ NA values for each column.
  expect_s3_class(suppressWarnings(importTomTomXML("inst/extdata/tomtom_ex/tomtom_bad.xml")), "data.frame")
})
