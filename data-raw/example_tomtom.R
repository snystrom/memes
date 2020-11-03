library(memes)
# Use FlyFactor database
db <- here::here("inst/extdata/flyFactorSurvey_cleaned.meme")

motif <- universalmotif::create_motif("CCRAAAW")
motif["name"] <- "example_motif"
motif["altname"] <- 'motif_ccraaaw'
example_tomtom <- runTomTom(motif, database = db, silent = FALSE)

usethis::use_data(example_tomtom)
