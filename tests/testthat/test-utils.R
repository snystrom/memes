skip_if(T, "delete when done building tests")

peaks <- "inst/extdata/peaks/peaks.tsv" %>%
  readr::read_tsv() %>%
  GRanges

dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

regions <- peaks %>%
  resize(200, "center") %>%
  data.frame %>%
  dplyr::mutate(type = c(rep("A", nrow(.)/2), rep("B", nrow(.)/2)),
                score = sample(1:10, nrow(.))) %>%
  GRanges

regions_list <- regions %>%
  split(mcols(.)$type)

get_sequence(regions, dm.genome)
get_sequence(regions, dm.genome, score_column = "score")

get_sequence(regions_list, dm.genome)
get_sequence(regions_list, dm.genome, score_column = "score") -> z
names(z$A)
