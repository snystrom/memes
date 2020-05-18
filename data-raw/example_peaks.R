## code to prepare `example_peaks` dataset goes here
example_peaks <- system.file("extdata/peaks/peaks.tsv", package = "dremeR") %>%
  readr::read_tsv() %>%
  GenomicRanges::GRanges()

usethis::use_data(example_peaks)
