library(dremeR)
library(GenomicRanges)
library(magrittr)

peaks <- system.file("extdata/peaks/e93_chr3.csv", package = "dremeR") %>%
  readr::read_csv() %>%
  GRanges

# These data use the dm3 reference genome
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

# compute summits using the summit_position column
summits <- peaks %>%
  plyranges::anchor_start() %>%
  plyranges::mutate(width = 1) %>%
  plyranges::shift_right(mcols(.)$summit_position)

# Get sequences in a 100bp window around the peak summit
summit_flank <- summits %>%
  plyranges::anchor_center() %>%
  plyranges::mutate(width = 100)

# Use E93 motif from motifdb
e93_motif <- MotifDb::MotifDb %>%
  MotifDb::query("Eip93F") %>%
  universalmotif::convert_motifs() %>%
  .[[1]]

# Rename the motif for simplicity
e93_motif["name"] <- "E93"

example_fimo <- summit_flank %>%
  get_sequence(dm.genome) %>%
  runFimo(motifs = e93_motif, thresh = 1e-3)

usethis::use_data(example_fimo)
