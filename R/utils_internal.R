#' Returns valid path to meme bin directory or supported meme utility
#'
#' @param path user override path to meme bin (optional)
#' @param util one of c(NULL,"dreme","ame","fimo","tomtom")
#'
#' @return valid path to meme/bin or meme utility
#'
#' @examples
#'
#' @noRd
handle_meme_path <- function(path = NULL, util = NULL){
  f <- dotargs::build_path_handler(environment_var = "MEME_PATH",
                                       option_name = "meme_bin",
                                       default_path = "~/meme/bin",
                                       utils = c("dreme", "ame", "fimo", "tomtom")
                                       )
  f(path, util)
}


#' Converts xml attrs to data-frame
#'
#' @param xml xml object
#'
#' @return
#'
#' @examples
#'
#' @noRd
attrs_to_df <- function(xml, ...) {
  # parsing DREME output
  # converts xml attributes to dataframe
  # where each column is an attribute
  xml2::xml_attrs(xml) %>%
    data.frame() %>%
    t() %>%
    data.frame(row.names = NULL, ...)
}
