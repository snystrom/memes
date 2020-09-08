skip_if_not(meme_is_installed(), "Only works with dreme install")

# Setup
fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes") %>%
  Biostrings::readDNAStringSet()
fa2 <- system.file("extdata/fasta_ex/fa2.fa", package = "memes") %>%
  Biostrings::readDNAStringSet()

dreme_out <- runDreme(fa, "shuffle", e = 39)

test_that("runDreme works", {
  # won't detect @ normal threshold
  expect_null(runDreme(fa, "shuffle"))

  expect_true(is_dreme_results(dreme_out))

  # test dreme catch error if fasta empty
  empty_fa <- tempfile()
  file.create(empty_fa)
  expect_error(suppressMessages(runDreme(empty_fa, "shuffle"), "No sequences"))
  file.remove(empty_fa)

})

test_that("Dreme flags correctly parsed", ~{
  expect_equal(prepareDremeFlags(input = "input.fa", control = "shuffle"),
               c("-p", "input.fa"))
  expect_equal(prepareDremeFlags(input = "input.fa", control = "background.fa"),
               c("-p", "input.fa", "-n", "background.fa"))
})

test_that("dreme_results & universalmotif list validators work", {
  # Test dreme_results validators
  spec <- new_dreme_results()
  bad <- spec[,1:5]

  expect_true(is_dreme_results(dreme_out))
  expect_false(is_dreme_results(bad))

  expect_null(error_dreme_results(dreme_out))
  expect_error(error_dreme_results(spec), "motif column is empty")
  expect_error(error_dreme_results(bad), "Missing columns")

  # test universalmotif validators
  expect_true(is_universalmotif_list(list(universalmotif::create_motif())))
  expect_false(is_universalmotif_list(list(universalmotif::create_motif(), "no")))
})

test_that("runDreme dispatch works", {
  # test multiple dispatch
  peaks <- system.file("extdata/peaks/peaks.tsv", package = "memes") %>%
    readr::read_tsv(col_types = c("seqnames" = "c",
                                  "start" = "i",
                                  "end" = "i")) %>%
    GenomicRanges::GRanges()
  dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

  seq <- get_sequence(peaks, dm.genome)
  seq <- get_sequence(peaks, dm.genome)
  fa <- write_fasta(seq)

  ## test fasta path input
  expect_equal(runDreme(fa, "shuffle", e = 49), runDreme(seq, "shuffle", e = 49))

  ## test stringset list input
  seq_by_type <- peaks %>%
    data.frame %>%
    dplyr::mutate(type = c(rep("A", nrow(.) / 2), rep("B", nrow(.) / 2))) %>%
    GenomicRanges::GRanges() %>%
    split(mcols(.)$type) %>%
    # genome is dm3
    get_sequence(dm.genome)

  # Use all input w/ shuffled background sequence
  expect_named(runDreme(seq_by_type, "shuffle", e = 70), c("A", "B"))
  # use "A" as background
  expect_named(runDreme(seq_by_type, "A", e = 70), "B")
  # Use invalid background
  expect_error(runDreme(seq_by_type, "d", e = 70), "input names: d")

  # test that control names are all inside names of input
  # ie that if name isn't in name, throws error
  expect_error(runDreme(seq_by_type, c("B", "C"), e = 70), "names passed to control do not exist")
  expect_error(runDreme(seq_by_type, c("shuffle", "C"), e = 70), "names passed to control do not exist")
  expect_error(runDreme(seq_by_type, c("shuffle", "B"), e = 70), "names passed to control do not exist")

})

test_that("runDreme error checking works", {
  ## test error checking
  expect_error(suppressMessages(runDreme(fa, "shuffle", et = 39)),
               "\"e\" instead of", class = "error")
  ## Ensure error catch works with list input also
  expect_error(suppressMessages(runDreme(seq_by_type, "A", et = 70)),
               "\"e\" instead of", class = "error")


})

