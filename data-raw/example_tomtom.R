library(dremeR)
db <- system.file("extdata/db/fly_factor_survey_id.meme", package = "dremeR")

motif <- universalmotif::create_motif("CCRAAAW")
motif["name"] <- "example_motif"
motif["altname"] <- 'motif_ccraaaw'
example_tomtom <- runTomTom(motif, database = db, silent = F)

usethis::use_data(example_tomtom)
