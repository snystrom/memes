library(memes)
data("example_dreme")
# Use FlyFactor database
db <- here::here("inst/extdata/flyFactorSurvey_cleaned.meme")
example_dreme_tomtom <- runTomTom(example_dreme, database = db, silent = FALSE)

usethis::use_data(example_dreme_tomtom)
