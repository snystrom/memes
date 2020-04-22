
#' Convert tomtom query entries to meme_df format
#'
#' @param tt_xml tomtom xml2 data
#'
#' @return
#'
#' @noRd
tomtom_query_motif_dfs <- function(tt_xml){
  tt_motif_list <- tt_xml %>%
    xml2::xml_find_all("//queries") %>%
    xml2::xml_children() %>%
    purrr::map(tomtom_xml_motif_to_universalmotif, tt_xml)

  data <- universalmotif_to_meme_df(tt_motif_list) %>%
    dplyr::mutate(query_idx = (1:nrow(.) - 1),
      db_idx = purrr::map_int(tt_motif_list, ~{
      .x@extrainfo["db"] %>%
        as.integer()
    }))

  return(data)
}

#' Title
#'
#' @param entry motif XML entry from tomtom.xml
#' @param tt_xml tomtom xml2 data structure (used to grab metadata)
#'
#' @return universalmotif object w/ metadata of entry
#'
#' @noRd
tomtom_xml_motif_to_universalmotif <- function(entry, tt_xml){
  data <- attrs_to_df(entry, stringsAsFactors = FALSE) %>%
    dplyr::mutate_at(c("length", "nsites"), as.integer) %>%
    dplyr::mutate_at(c("evalue"), as.double)

  pfm <- t(get_probability_matrix(entry))

  # Background frequency
  tt_run_info <- xml2::xml_children(tt_xml)[1]
  bkg <- dreme_get_background_freq(tt_run_info)

  motif <- universalmotif::create_motif(pfm,
                               name = data$id,
                               altname = check_col(data, "alt", character(0)),
                               bkg = bkg,
                               pval = check_col(data, "pvalue"),
                               nsites = check_col(data, "nsites"),
                               eval = check_col(data, "evalue"),
                               extrainfo = c("db" = data$db)
                               )


  return(motif)
}

#' Return valid data types
#'
#' Return empty data type if value is missing, else return value.
#'
#' Useful for passing values to a class definition with type-checking if you
#' can't predict which values may be missing. values of `col` will
#'
#' @param df data.frame
#' @param col column name to check. Values of `col` will be coerced to type in `type`
#' @param type data type to return (typically one of `character(0)`, `integer(0)`, etc. default: `numeric(0)`)
#'
#' @return value if defined, empty data type if undefined
#'
#' @examples
#' \dontrun{
#' df_undef <- data.frame(a = NULL)
#' check_col(df_undef$a)
#' df_def <- data.frame(a = 1)
#' check_col(df_undef$a)
#' }
check_col <- function(df, col, type = numeric(0)){
  val <- ifelse(!is.null(df[[col]]), df[[col]], type)
  class(val) <- class(type)
  return(val)
}
