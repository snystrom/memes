library(dremeR)
library(universalmotif)
library(magrittr)

flyFactorDb <- MotifDb::MotifDb %>%
  MotifDb::query("FlyFactorSurvey")

flyFactorMotifs <- flyFactorDb %>%
  convert_motifs()

flyFactor_data <- flyFactorMotifs %>%
  as_universalmotif_dataframe()

#############
# Add motifs that are missing from MotifDb FlyFactor Survey data
meme_ffs <- read_meme("inst/extdata/db/fly_factor_survey_id.meme") %>%
  as_universalmotif_dataframe()

flyFactor_data %<>%
  dplyr::mutate(id = gsub("_FBgn\\d+", "", altname))

missing_ff <- meme_ffs %>%
  dplyr::filter(!(name %in% flyFactor_data$id)) %>%
  dplyr::mutate(symbol = gsub("(.+)-?_.+", "\\1", name)) %>%
  dplyr::mutate(symbol = gsub("-.+", "", symbol)) %>%
  dplyr::mutate(symbol = ifelse(symbol == "Blimp", "Blimp-1", symbol)) %>%
  dplyr::mutate(isoform = gsub("_.+", "", name)) %>%
  dplyr::mutate(fbgn = gsub("_.+", "", altname)) %>%
  tidyr::unite(name, c("name", "fbgn"))

missing_ffs_noentries <- missing_ff %>%
  dplyr::filter(!(symbol %in% flyFactor_data$name)) %>%
  dplyr::mutate(fbgn = gsub("_.+", "", altname)) %>%
  dplyr::select(-altname) %>%
  dplyr::rename(altname = symbol) %>%
  dplyr::select(-isoform)

missing_ffs_motifs <- as_universalmotif(missing_ffs_noentries)
# Add datasource entry
missing_ffs_motifs <- purrr::map(missing_ffs_motifs, ~{
  .x["extrainfo"] <- c("dataSource" = "FlyFactorSurvey")
  return(.x)
})

##########

# Tidy motifdb data
flyFactor_data %<>%
  dplyr::rename("altname" = "name",
                "name" = "altname")

flyFactor_data %<>%
  # Critical to set remove = FALSE to keep the `name` column
  tidyr::separate(name, c("tfid"), remove = FALSE, extra = "drop") %>%
  # Only use the tfid if the altname contains an FBgn
  dplyr::mutate(altname = ifelse(grepl("^FBgn", altname), tfid, altname)) %>%
  dplyr::mutate(name = gsub("_FBgn\\d+", "", name))


flyFactorMotifs <- as_universalmotif(flyFactor_data)
# Combine data together
flyFactorMotifs <- c(flyFactorMotifs, missing_ffs_motifs)


write_meme(flyFactorMotifs, "inst/extdata/flyFactorSurvey_cleaned.meme")
