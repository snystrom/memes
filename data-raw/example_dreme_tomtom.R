library(dremeR)
data("dreme_example")
db <- system.file("extdata/db/fly_factor_survey_id.meme", package = "dremeR")
example_dreme_tomtom <- runTomTom(dreme_example, database = db, silent = F)

usethis::use_data(example_dreme_tomtom)
