test_that("Dreme flags correctly parsed", ~{
  expect_equal(prepareDremeFlags(input = "input.fa", control = "shuffle"),
               c("-p", "input.fa"))
  expect_equal(prepareDremeFlags(input = "input.fa", control = "background.fa"),
               c("-p", "input.fa", "-n", "background.fa"))
})

