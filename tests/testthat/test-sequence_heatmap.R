
test_that("sequence heatmap functions work", {
  # length error works
  expect_error(sequence_to_df(c("AAA", "AA")), "must be equal length")
  expect_error(sequence_to_df(c("AAA", "AAA")), NA)
  
  good_df <- sequence_to_df(c("AAA", "ATA"))
  expect_equal(c(1,1,1,2,2,2), good_df$id)
  expect_equal(c(1,2,3,1,2,3), good_df$position)
  expect_equal(c("A", "A", "A", "A", "T", "A"), good_df$letters)
  
  # alph detection works
  expect_error(sequence_to_heatmap(c("AAA"), "DNA"), NA)
  expect_error(sequence_to_heatmap(c("AAA"), "RNA"), NA)
  expect_error(sequence_to_heatmap(c("AAA"), "AA"), NA)
  expect_error(sequence_to_heatmap(c("AAA"), "AAA"))
  
  #TODO: deal with named sequences? BStringsList?
  # set_names(named_seq, NULL)?
})
  