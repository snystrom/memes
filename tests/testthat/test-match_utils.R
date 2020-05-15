data("example_tomtom")

# test force_best_match
expect_error(
  force_best_match(example_tomtom, c("bad"= "bad")),
  "invalid names: bad"
)

expect_error(
  force_best_match(example_tomtom, c("example_motif"= "bad")),
  "bad is not found within "
)

expect_equal(
  force_best_match(example_tomtom, c("example_motif"= "Lag1_Cell"))$best_match_name,
  "Lag1_Cell"
)
