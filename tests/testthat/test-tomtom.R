skip_if(T, "not developed yet")
dreme_file <- dotargs::expected_outputs(c("txt", "xml", "html"), "dreme", "inst/extdata/fasta_ex/fa1_vs_shuffle/")

tt_files <- runTomTom(dreme_file$txt, database = "inst/extdata/db/fly_factor_survey_id.meme", thresh = 10)

###########




dreme_file <- dotargs::expected_outputs(c("txt"), "dreme", "inst/extdata/fasta_ex/fa1_vs_shuffle/")
dreme_file <- duplicate_file(dreme_file$txt)
runTomTom(dreme_file, "inst/extdata/db/fly_factor_survey_id.meme")
runTomTom(dreme_out, "inst/extdata/db/fly_factor_survey_id.meme")

mList <- dreme_out$motifs
m <- mList[[1]]
options(tomtom_db = "inst/extdata/db/fly_factor_survey_id.meme")
runTomTom(mList)

# FAILS because no motif detected, Need method to skip if NA
runTomTom(m)
# Works because m2 is the only one that gets matched
runTomTom(mList[[2]])

#ume <- as.environment("package:universalmotif")
#eval(my_write_meme(dreme_out$motifs, "myMeme.meme"), envir = ume)
#universalmotif::write_meme(dreme_out$motifs, "hisMeme.meme")

###
# dealing w/ empty tomtom results (ie no match ever)
fa <- duplicate_file("inst/extdata/fasta_ex/fa1.fa")
dreme_out <- runDreme(fa, "shuffle", e = 39)
# this first motif matches nothing
m <- dremeR:::write_meme_list(dreme_out$motifs[[1]])

options(tomtom_db = "inst/extdata/db/fly_factor_survey_id.meme")
runTomTom(m, outdir = "tt_nomatch_dev")

nomatch_xml <- "tt_nomatch_dev/tomtom.xml" %>%
  xml2::read_xml()

###



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

