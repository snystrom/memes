library(dremeR)
data("example_dreme")
db <- system.file("extdata/db/fly_factor_survey_id.meme", package = "dremeR")
example_dreme_tomtom <- runTomTom(example_dreme, database = db, silent = F)

usethis::use_data(example_dreme_tomtom)
