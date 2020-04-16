skip_if(T, "delete when done building tests")

peaks <- "inst/extdata/peaks/peaks.tsv" %>%
  readr::read_tsv() %>%
  GRanges
get_sequence(peaks, dm.genome)
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

x <- Biostrings::getSeq(dm.genome, peaks)
x

get_sequence(peaks, dm.genome)
