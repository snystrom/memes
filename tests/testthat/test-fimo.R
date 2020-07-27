ffdb <- universalmotif::read_meme("inst/extdata/db/fly_factor_survey_id.meme")
df <- memes:::universalmotif_to_meme_df(ffdb)
df %<>%
  dplyr::filter(grepl("Eip93F", id))
e93_motif <- ffdb %>%
  universalmotif::filter_motifs(name = df$id)

universalmotif::write_meme(e93_motif, "inst/extdata/db/e93_motif.meme")

suppressPackageStartupMessages(library(GenomicRanges))

peaks <- "inst/extdata/peaks/peaks.tsv" %>%
  readr::read_tsv() %>%
  GRanges



e93_motif <- universalmotif::read_meme("inst/extdata/db/e93_motif.meme")
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

fimo_res <- peaks %>%
  resize(200, fix = "center") %>%
  get_sequence(dm.genome) %>%
  runFimo(e93_motif)

fimo_res$stdout[1] %>%
  readr::write_lines("fimo_test.tsv", sep = "")

#motif_input(e93_motif) -> x

readr::read_tsv(fimo_res$stdout) %>%
  dplyr::rename_all(~{gsub("-", "", .)}) %>%
  dplyr::rename("seqnames" = "sequence_name") %>%
  GRanges


####

# errror no names in header
seq <- universalmotif::create_sequences()
motif <- universalmotif::create_motif()
runFimo(seq, motif, parse_genomic_coord = FALSE)

# run, but message "no matches detected"
names(seq) <- seq_along(seq)
runFimo(seq, motif, parse_genomic_coord = FALSE)
runFimo(seq, motif)

# Suggest psp instead of psspp
runFimo(seq, motif, psspp = 'x')
