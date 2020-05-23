library(dremeR)
library(dplyr)

link <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97956&format=file&file=GSE97956%5FgeneFPKMs%5FL3%5F24hr%5F44hr%2Exlsx"
path <- curl::curl_download(link, destfile = tempfile())

full_rnaseq <- readxl::read_excel(path) %>%
  tidyr::pivot_longer(-matches("gene_symbol"), names_to = "time", values_to = "fpkm") %>%
  dplyr::filter(time != "44hr") %>%
  dplyr::mutate(time = case_when(time == "L3" ~ "Early",
                                 time == "24hr" ~ "Late")) %>%
  dplyr::rename(symbol = gene_symbol) %>%
  dplyr::group_by(symbol) %>%
  dplyr::filter(max(fpkm) != 0)

data("flyFactorMotifs", package = "dremeR")

ff_df <- flyFactorMotifs %>%
  as_universalmotif_dataframe()

tf_rnaseq <- full_rnaseq %>%
  dplyr::mutate(symbol = case_when(symbol == "Bap" ~ "bap",
                                   T ~ symbol)) %>%
  dplyr::filter(symbol %in% ff_df$altname)

set.seed(123)
no_tf_rnaseq <- full_rnaseq %>%
  dplyr::filter(!(symbol %in% ff_df$altname)) %>%
  dplyr::slice(sample(1:nrow(.), size = 200, replace = F))

example_rnaseq <- tf_rnaseq %>%
  dplyr::bind_rows(no_tf_rnaseq)

usethis::use_data(example_rnaseq, overwrite = T)
