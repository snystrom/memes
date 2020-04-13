skip_if(T, "Only works with dreme install")

fa <- "inst/extdata/fasta_ex/fa1.fa"
fa2 <- "inst/extdata/fasta_ex/fa2.fa"

dreme_out <- runDreme(fa, "shuffle", e = 39)

###
#test tomtom
dreme_file <- dotargs::expected_outputs(c("txt", "xml", "html"), "dreme", "inst/extdata/fasta_ex/fa1_vs_shuffle/")

tt_files <- runTomTom(dreme_file$txt, database = "inst/extdata/db/fly_factor_survey_id.meme", thresh = 10)

# dev tt_result_parse
# want data.frame w/ results & motif_match col w/ universalmotif object

#tomtom_xml_path <- tt_files$xml
tomtom_xml_path <- "inst/extdata/dreme_example/tomtom/tomtom.xml"
tt_xml <- xml2::read_xml(tomtom_xml_path)

queries <- xml2::xml_find_all(tt_xml, "//queries") %>%
  xml2::xml_children() %>%
  attrs_to_df() %>%
  dplyr::mutate(idx = 1:nrow(.))

matches <- xml2::xml_find_all(tt_xml, "//matches") %>%
  xml2::xml_children()

targets <- xml2::xml_find_all(tt_xml, "//targets") %>%
  xml2::xml_children()

target_df <- targets %>%
  attrs_to_df(stringsAsFactors = F) %>%
  dplyr::mutate(length = as.integer(length),
                nsites = as.integer(length))

target_pfms <- targets %>%
  purrr::map(get_probability_matrix) %>%
  purrr::map(t)

target_list <- target_df %>%
  split(.$id)

pfmList <- purrr::map2(target_list, target_pfms, ~{
  universalmotif::create_motif(.y,
                               type = "PCM",
                               name = .x$id,
                               altname = .x$alt,
                               nsites = .x$nsites)

})

target_df$motif_match <- pfmList
#####
stop("dev below")
runDreme(fa, fa2)
