
db <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE)
fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes", mustWork = TRUE)
dreme_out <- runDreme(fa, "shuffle", e = 39, outdir = tempdir())

fa2 <- universalmotif::create_sequences(seqnum = 2, seqlen = 100, rng.seed = 123)
names(fa2) <- seq_along(fa2)
meme_out <- runMeme(fa2, parse_genomic_coord = FALSE)
