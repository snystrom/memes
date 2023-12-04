meme.full <- system.file("extdata", "meme_full.txt", package = "universalmotif")

test_that("MEME Import errors correctly", {
  suppressWarnings(
    expect_error(importMeme(meme.full, parse_genomic_coord = TRUE), 
                 "Problem parsing genomic coordinates"),
  )
  
})

test_that("prepareMemeFlags warns about objfun", {
  # TODO:
  expect_warning(prepareMemeFlags(control = NA), "cannot use classic")
  expect_warning(prepareMemeFlags(control = NA, objfun = "classic"), "cannot use classic")
  expect_warning(prepareMemeFlags(control = NA, objfun = "se"), NA)
})
