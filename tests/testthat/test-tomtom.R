dreme_file <- dotargs::expected_outputs(c("txt", "xml", "html"), "dreme", "inst/extdata/fasta_ex/fa1_vs_shuffle/")

tt_files <- runTomTom(dreme_file$txt, database = "inst/extdata/db/fly_factor_survey_id.meme", thresh = 10)
#########
tomtom_xml_path <- "inst/extdata/dreme_example/tomtom/tomtom.xml"

tt <- parseTomTom(tomtom_xml_path)

tt %>%
  dplyr::group_by(id, alt) %>%
  tidyr::nest() %>%
  dplyr::mutate(best_match = purrr::map(data, ~{
    .x %>%
      dplyr::filter(match_evalue == min(match_evalue)) %>%
      .$match_motif
  })) -> r

r %>%
  tidyr::unnest(best_match)

get_best_match_motif <- function(res, by = "match_evalue", motif_colname = "match_motif"){
  by = rlang::quo_name(by)

  res %>%
    dplyr::filter(!!by == min(!!by)) %>%
    .[motif_colname]
}

tt %>%
  dplyr::group_by(id, alt) %>%
  tidyr::nest() %>%
  dplyr::mutate(best_match = purrr::map(data, get_best_match_motif)) %>%
  tidyr::unnest(best_match) -> rr

