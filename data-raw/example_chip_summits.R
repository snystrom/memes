library(dplyr)
library(plyranges)

example_chip_summits <- system.file("extdata/peaks/e93_chr3.csv", package = "memes") %>%
  readr::read_csv() %>%
  GRanges() %>%
  anchor_start() %>%
  mutate(width = 1) %>%
  shift_right(.$summit_position) %>%
  GRanges %>%
  select(id, peak_binding_description, e93_sensitive_behavior)

usethis::use_data(example_chip_summits)
