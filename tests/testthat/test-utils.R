withr::local_options(list(meme_bin = NULL))
withr::local_envvar(list(MEME_BIN = ""))
test_that("get_sequence works", {

  peaks <- system.file("extdata/peaks/peaks.tsv", package = "memes", mustWork = TRUE) %>%
    readr::read_tsv(col_types = c("seqnames" = "c", "start" = "i", "end" = "i")) %>%
    GenomicRanges::GRanges()

  dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

  regions <- peaks %>%
    GenomicRanges::resize(10, "center") %>%
    data.frame %>%
    dplyr::mutate(type = c(rep("A", nrow(.)/2), rep("B", nrow(.)/2)),
                  score = sample(1:10, nrow(.))) %>%
    GenomicRanges::GRanges()

  regions_list <- regions %>%
    split(mcols(.)$type)

  seqs <- get_sequence(regions, dm.genome)
  expect_equal(unique(Biostrings::width(seqs)), 10)
  expect_s4_class(seqs, "DNAStringSet")

  # score is added correctly
  seqs_score <- get_sequence(regions, dm.genome, score_column = "score")

  expect_true(seqs[1][[1]] == seqs_score[1][[1]])
  expect_equal(names(seqs_score)[1], paste(names(seqs)[1], regions[1]$score))

  # list dispatch works
  seqs_list <- get_sequence(regions_list, dm.genome)

  expect_true(all(seqs_list$A == seqs[seq(1, length(seqs)/2)]))
  expect_true(all(names(seqs_list$A) %in% names(seqs)))
  expect_true(all(names(seqs_list$B) %in% names(seqs)))

  # add_sequence works
  regions_with_seq <- regions %>%
    add_sequence(dm.genome)

  expect_true(all(regions_with_seq$sequence == seqs))
})

test_that("meme_is_installed works", {
  # catches if fails when dir not exist
  expect_false(meme_is_installed("/this/dir/doesnot/exist"))
})
