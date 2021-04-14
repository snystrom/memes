test_that("streme works", {
  withr::local_options(meme_bin = "~/meme/bin")
  seqs <- universalmotif::create_sequences("CCRAAAW", rng.seed = 123)
  expect_s3_class(runStreme(seqs, "shuffle"), "universalmotif_df")
  expect_s3_class(runStreme(seqs, control = NA, objfun = "cd"), "universalmotif_df")
  # Error if control is set when objfun = "cd"
  expect_error(runStreme(seqs, control = "shuffle", objfun = "cd"), "must be NA")
  
  # Test that error suggestion works
  expect_error(suppressMessages(runStreme(seqs, "shuffle", versionn = TRUE)), 
               "\"version\" instead of:")
  expect_error(suppressMessages(runStreme(seqs, "shuffle", seeed = 123)), 
               "\"seed\" instead of:")
 
  few_seqs <- universalmotif::create_sequences(seqnum = 5, rng.seed = 321) 
  suppressMessages(expect_message(runStreme(few_seqs, "shuffle"), 
               "Warning: No hold-out set"))
  suppressMessages(expect_message(runStreme(few_seqs, "shuffle"), 
               "Warning: Ignoring <pvt>"))
  
  
  # Test list input
  seqlist <- list("one" = seqs,
       "two" = universalmotif::create_sequences("AATAATT", rng.seed = 321)
       )
  
  expect_type(runStreme.list(seqlist, "two"), "list")
  expect_type(runStreme.list(seqlist, "shuffle"), "list")
  expect_named(runStreme.list(seqlist, "shuffle"), c("one", "two"))
})
