library(memes)
library(magrittr)
library(universalmotif)

flyFactorDb <- MotifDb::MotifDb %>%
  MotifDb::query("FlyFactorSurvey")

flyFactorMotifs <- flyFactorDb %>%
  convert_motifs()

flyFactor_data <- flyFactorMotifs %>%
  to_df()

flyFactor_data %<>%
  dplyr::rename("altname" = "name",
                "name" = "altname") %>% 
  # Critical to set remove = FALSE to keep the `name` column
  tidyr::separate(name, c("tfid"), remove = FALSE, extra = "drop") %>%
  # Only use the tfid if the altname contains an FBgn
  dplyr::mutate(altname = ifelse(grepl("^FBgn", altname), tfid, altname)) %>% 
  dplyr::mutate(name = gsub("_FBgn\\d+", "", name)) %>%
  # rename all "da" instances using their tfid value instead
  dplyr::mutate(altname = ifelse(altname == "da", tfid, altname))

swap_alt_id <- c("CG6272", "Clk", "Max", "Mnt", "Jra")
remove <- "Bgb"

flyFactor_data %<>%
  dplyr::mutate(altname = ifelse(altname %in% swap_alt_id, tfid, altname)) %>%
  dplyr::filter(!(altname %in% remove))

# This operation takes a while to run on large motif lists
flyFactor_dedup <- remove_duplicate_motifs(flyFactor_data)
flyFactorMotifs_final <- to_list(flyFactor_dedup, extrainfo = FALSE)

write_meme(flyFactorMotifs_final, "inst/extdata/flyFactorSurvey_cleaned.meme")

flyFactor_data %>% 
  dplyr::filter(consensus == "MMCACCTGYYV") %>% 
  to_list() %>% 
  write_meme("inst/extdata/flyFactor_dups.meme")
