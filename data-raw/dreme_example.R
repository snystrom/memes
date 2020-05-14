library(dremeR)

peaks <- system.file("extdata/peaks/e93_chr3.csv", package = "dremeR") %>%
  readr::read_csv() %>%
  GenomicRanges::GRanges

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

# split by response to E93 binding
by_sens <- summit_flank %>%
  split(mcols(.)$e93_sensitive_behavior) %>%
  get_sequence(dm.genome)

dreme_by_sens_vs_static <- runDreme(by_sens, "Static")

dreme_example <- dreme_by_sens_vs_static$Decreasing

#db <- system.file("extdata/db/fly_factor_survey_id.meme")
#tomtom_example <- runTomTom(dreme_example, database = db)
usethis::use_data(dreme_example)

