skip_if_not(meme_is_installed())

test_that("runFimo works", {
  peaks <- system.file("extdata/peaks/peaks.tsv", package = "memes") %>%
    readr::read_tsv(col_types = c("seqnames" = "c", "start" = "i", "end" = "i")) %>%
    GenomicRanges::GRanges()

  dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

  e93_motif <- universalmotif::create_motif("CCRAAAW", name = "E93", altname = "Eip93F")

  regions <- peaks %>%
    GenomicRanges::resize(200, fix = "center")

  fimo_res <- regions %>%
    get_sequence(dm.genome) %>%
    runFimo(e93_motif, thresh = 1e-3, skip_matched_sequence = FALSE) %>%
    add_sequence(dm.genome)

  expect_s4_class(fimo_res, "GenomicRanges")

  # Ensure add_sequence matches the detected sequences
  fimo_res %<>% add_sequence(dm.genome)
  expect_equal(fimo_res$matched_sequence, fimo_res$sequence)

  # errror no names in header
  seq <- universalmotif::create_sequences()
  motif <- universalmotif::create_motif()
  expect_message(try(runFimo(seq, motif, parse_genomic_coord = FALSE), silent = TRUE), "header")
  expect_error(suppressMessages(runFimo(seq, motif, parse_genomic_coord = FALSE)))

  # run, but message "no matches detected"
  names(seq) <- seq_along(seq)
  expect_null(suppressMessages(runFimo(seq, motif, thresh = 1e-10)))
  expect_message(runFimo(seq, motif, thresh = 1e-10), "No matches were detected")

  # Suggest psp instead of psspp
  expect_error(suppressMessages(runFimo(seq, motif, psspp = 'x')), class = "usethis_error", regexp = "Invalid flags")
  expect_message(try(runFimo(seq, motif, psspp = 'x'), silent = TRUE), "Error processing command line options")

})
