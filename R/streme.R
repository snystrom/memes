#' Title
#'
#' @param input 
#' @param control 
#' @param outdir output directory of results. (default: "auto" uses temp directory)
#' @param objfun one of c("de", "cd"). Default: "de" for differential
#'   enrichment. "cd" for central distance (control must be set to NA for "cd").
#' @param alph one of c("dna", "rna", "protein") or a path to a MEME format alph file.
#' @param meme_path path to "meme/bin"
#' @param silent Whether to print stdout to console (default: TRUE)
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
runStreme <- function(input, control, outdir = "auto", objfun = "de", 
                      alph = "dna", meme_path = NULL, silent = TRUE, ...){

  # Handle multiple input types by multiple dispatch
  # input & control will be coerced to file paths
  input <- sequence_input(input)
  
  if (all(is.na(control))) {
    control <- NA
  } else {
    control <- sequence_input(control)
  }
  
  # no control allowed if objfun = "cd"
  if(objfun == "cd" & !is.na(control)) {
    stop("control must be NA if objfun = 'cd'", call. = FALSE)
  }
  
  if (outdir == "auto") {outdir <- outdir_name(input, control)}

  flags <- prepareStremeFlags(input = input, control = control, 
                              outdir = outdir, alph = alph, ...)
  
  command <- search_meme_path(meme_path, "streme")
  ps_out <- processx::run(command, flags, spinner = TRUE, 
                          error_on_status = FALSE)
  # TODO: remove
  #ps_out$outdir <- outdir
  #return(ps_out)
  
  #TODO: process_check_error
  #message(ps_out$stdout)
  #message(ps_out$stderr)
  #TODO: return NULL if no motifs found
  
  streme_out <- cmdfun::cmd_file_expect("streme", c("txt", "html", "xml"), 
                                        outdir = outdir)

  streme_results <- parseStreme(streme_out$xml)

  return(streme_results)
}

prepareStremeFlags <- function(input, control, outdir, alph, ...){
  argsDict <- c("input" = "p",
               "control" = "n",
               "outdir" = "oc")
  
  # handle alphabet assignment
  alph_flags <- meme_alph_to_args(alph) %>%
    cmdfun::cmd_list_interp()

  flagsList <- cmdfun::cmd_args_all(drop = "alph") %>%
    cmdfun::cmd_list_interp(argsDict) %>% 
    cmdfun::cmd_list_drop(c("n" = "shuffle"))

  flags <- c(flagsList, alph_flags) %>%
    cmdfun::cmd_list_to_flags(prefix = "--")
  return(flags)
}

#' Import streme results from xml
#'
#' @param xml path to streme xml
#'
#' @return universalmotif_df
parseStreme <- function(xml){
  streme_stats <- streme_motif_stats(xml)

  pfms <- streme_to_pfm(xml)

  streme_stats$motif <- pfms
  # Convert to motif_df format
  # suppressing messages about adding empty motif slots
  suppressMessages(
    universalmotif::update_motifs(streme_stats, extrainfo = TRUE)
    )
}

#' Get universalmotifs from streme xml
#'
#' @param streme_xml_path 
#'
#' @return
#' @noRd
streme_to_pfm <- function(streme_xml_path){
  streme_xml <- xml2::read_xml(streme_xml_path)

  streme_run_info <- xml2::xml_children(streme_xml)[1] %>%
    xml2::xml_children()

  background_freq <- dreme_get_background_freq(streme_run_info)

  motifs_matrix <- xml2::xml_children(streme_xml)[2] %>%
    xml2::xml_children() %>% 
    purrr::map(get_probability_matrix) %>% 
    purrr::map(t)

  motif_stats_list <- streme_motif_stats(streme_xml_path) %>%
    split(.$name)

  pfmList <- purrr::map2(motif_stats_list, motifs_matrix, ~{
    universalmotif::create_motif(.y,
                                 type = "PCM",
                                 name = .x$name,
                                 altname = .x$altname,
                                 bkg = background_freq,
                                 pval = .x$pval,
                                 nsites = .x$nsites)

  })

  return(pfmList)
}

#' Parse stats for Streme results
#'
#' @param streme_xml_path path to streme.xml
#'
#' @return data.frame with full statistics for each streme motif
#' @noRd
streme_motif_stats <- function(streme_xml_path){
  int_cols <- c("width", "initial_width")
  dbl_cols <- c("score_threshold", "train_pos_count", "train_neg_count", 
                "train_log_pvalue", "train_pvalue", "train_dtc", 
                "train_bernoulli", "test_pos_count", "test_neg_count", 
                "test_log_pvalue", "test_pvalue",  "test_dtc", "test_bernoulli", 
                "elapsed_time")
  
  streme_xml <- xml2::read_xml(streme_xml_path)
  
  xml2::xml_children(streme_xml)[2] %>%
    xml2::xml_children() %>% 
    attrs_to_df() %>% 
    dplyr::mutate("id" = gsub("-", "_", id) %>% 
                    gsub("(\\d)", "m\\1", .)) %>% 
    dplyr::rename("name" = "id",
                  "altname" = "alt") %>% 
    dplyr::mutate_at(int_cols, as.integer) %>% 
    dplyr::mutate_at(dbl_cols, as.double) %>% 
    dplyr::mutate(nsites = test_pos_count + train_pos_count) %>% 
    # Put test variables first (do before rename)
    dplyr::select("name", "altname", "width", 
                  "initial_width", "seed", "nsites", "score_threshold",
                  tidyselect::starts_with("test_"), 
                  tidyselect::starts_with("train_"), 
                  tidyselect::everything()) %>% 
    # Change the "test_" columns for compatibility w/ universalmotif_df
    dplyr::rename_with(~{gsub("test_", "", .)},
                     tidyselect::starts_with("test_"), 
                     ) %>% 
    # pvalue -> pval for universalmotif compatability
    dplyr::rename_with(~{gsub("pvalue", "pval", .)},
                     tidyselect::contains("pvalue"), 
                     ) 
}
