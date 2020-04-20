
tomtom_query_motif_dfs <- function(tt_xml){
  tt_motif_list <- tt_xml %>%
    xml2::xml_find_all("//queries") %>%
    xml2::xml_children() %>%
    purrr::map(tomtom_xml_motif_to_dataframe, tt_xml)

  tt_motifs <- purrr::map(tt_motif_list, "motif")
  tt_data <- purrr::map(tt_motif_list, ~{dplyr::select(.x, -"motif")}) %>%
    dplyr::bind_rows()

  tt_data$motif <- tt_motifs

  tt_data %<>%
    dplyr::mutate(query_idx = (1:nrow(.) - 1))
    dplyr::rename("db_idx" = "db")
  return(tt_data)
}

tomtom_xml_motif_to_dataframe <- function(entry, tt_xml){
  df <- attrs_to_df(entry, stringsAsFactors = FALSE) %>%
    dplyr::mutate_at(c("db", "length", "nsites"), as.integer) %>%
    dplyr::mutate_at(c("evalue"), as.double)

  pfm <- t(get_probability_matrix(entry))

  # Background frequency
  tt_run_info <- xml2::xml_children(tt_xml)[1]
  bkg <- dreme_get_background_freq(tt_run_info)

  motif <- universalmotif::create_motif(pfm,
                               name = df$id,
                               altname = check_col(df, "alt", character(0)),
                               bkg = bkg,
                               pval = check_col(df, "pvalue"),
                               nsites = check_col(df, "nsites"),
                               eval = check_col(df, "evalue")
                               )

  df$motif <- motif

  return(df)
}

check_col <- function(df, col, type = numeric(0)){
  ifelse(!is.null(df[[col]]), df[[col]], type)
}
