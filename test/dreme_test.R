# dreme test
txt_path <- "inst/extdata/dreme_example/dreme.txt"

do <- universalmotif::read_meme("inst/extdata/dreme_example/dreme.txt")
do2 <- motifStack::importMatrix("inst/extdata/dreme_example/dreme.txt", format = "meme")

parse_dreme_txt <- function(path){
  #system2("tail", "-n+8", stdin = path, stdout = "./MYtest.meme")
  system2("tail", "-n+8", stdin = path, stdout = "./MYtest.meme_pregrep")
  system2("grep", c("-v \\#", "./MYtest.meme_pregrep"), stdout = "./MYtest.meme")
}

parse_dreme_txt(txt_path)

out_meme <- "MYtest.meme"

universalmotif::read_meme(out_meme) -> do
motifStack::importMatrix(out_meme) -> do


####################
do <- universalmotif::read_meme("inst/extdata/dreme_example/m01_DREME-1_DGCARC.meme")

system2("cat", "inst/extdata/dreme_example/*.meme", stdout = "allMotifs.meme")

do2 <- universalmotif::read_meme("allMotifs.meme")

do2[[1]]["name"]


do2 %>%
  {purrr::set_names(., purrr::map_chr(., ~{.["name"]}))}
