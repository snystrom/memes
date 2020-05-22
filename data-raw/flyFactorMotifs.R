library(dremeR)
library(universalmotif)

flyFactorDb <- MotifDb::MotifDb %>%
  MotifDb::query("FlyFactorSurvey")

flyFactorMotifs <- flyFactorDb %>%
  convert_motifs()

flyFactor_data <- flyFactorMotifs %>%
  as_universalmotif_dataframe()


flyFactor_data %<>%
  dplyr::rename("altname" = "name",
                "name" = "altname")

flyFactor_data %<>%
  # Critical to set remove = FALSE to keep the `name` column
  tidyr::separate(name, c("tfid"), remove = FALSE, extra = "drop") %>%
  # Only use the tfid if the altname contains an FBgn
  dplyr::mutate(altname = ifelse(grepl("^FBgn", altname), tfid, altname)) %>%
  mutate(name = gsub("_FBgn\\d+", "", name))

flyFactorMotifs <- as_universalmotif(flyFactor_data)

write_meme(flyFactorMotifs, "inst/extdata/flyFactorSurvey_cleaned.meme")
