context("XML import works")

test_that("DREME XML import",{
  dreme_xml <- system.file("extdata/dreme_example/dreme.xml", package = "memes")
  expect_s3_class(dreme_xml_data <<- importDremeXML(dreme_xml), "data.frame")
  expect_true(is_dreme_results(dreme_xml_data))
})

# testthat NULL tomtom result returns empty columns
test_that("TomTom Import works",{
  # good has matches
  tomtom_good_xml <- system.file("extdata/tomtom_ex/tomtom_good.xml", package = "memes")
  # bad has no detected matches
  tomtom_bad_xml <- system.file("extdata/tomtom_ex/tomtom_bad.xml", package = "memes")

  expect_s3_class(tomtom_good <<- importTomTomXML(tomtom_good_xml), "data.frame")
  # Should return query motifs w/ NA values for each column.
  expect_s3_class(suppressMessages(tomtom_bad <<- importTomTomXML(tomtom_bad_xml)), "data.frame")
  expect_message(importTomTomXML(tomtom_bad_xml), "TomTom detected no matches")

  # Ensure "query_idx", "target_idx", or "db_idx" cols don't appear in results
  expect_false(grepl("_idx", names(tomtom_good)))
  expect_false(grepl("_idx", names(tomtom_bad)))

  # TODO:
  # Check valid tomtom objects
})
