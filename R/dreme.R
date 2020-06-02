#' @export
#' @noRd
runDreme.list <- function(input, control, outdir = "auto", meme_path = NULL, silent = TRUE, ...){

  x <- sequence_input_control_list(input, control)
  input <- x$input
  control <- x$control

  res <- purrr::map(input, runDreme.default,
             control = control,
             outdir = outdir,
             meme_path = meme_path,
             silent = silent,
             ...
             )
  return(res)
}

#' @export
#' @noRd
runDreme.BStringSetList <- function(input, control, outdir = "auto", meme_path = NULL, silent = TRUE, ...){
  runDreme.list(as.list(input), control, outdir, meme_path, silent, ...)
}

#' @export
#' @noRd
runDreme.default <- function(input, control, outdir = "auto", meme_path = NULL, silent = TRUE, ...){

  # Handle multiple input types by multiple dispatch
  # input & control will be coerced to file paths
  input <- sequence_input(input)
  control <- sequence_input(control)

  if (outdir == "auto") {outdir <- outdir_name(input, control)}

  flags <- prepareDremeFlags(input = input, control = control, outdir = outdir, ...)

  command <- handle_meme_path(path = meme_path, util = "dreme")
  ps_out <- processx::run(command, flags, spinner = T, error_on_status = F)

  ps_out %>%
    process_check_error(help_fun = ~{dreme_help(command)},
                        user_flags = cmdlr::get_help_flag_names(flags))

  print_process_stdout(ps_out, silent = silent)

  n_motifs <- dreme_nmotifs_found(ps_out)

  if (n_motifs == 0) {return(NULL)}

  dreme_out <- cmdlr::expected_outputs(c("txt", "html", "xml"), "dreme", outdir = outdir)

  dreme_out %>%
    cmdlr::cmd_check_files_exist()

  dreme_results <- parseDreme(dreme_out$xml)

  return(dreme_results)
}

#' Prepare flags for Dreme
#'
#' @param input input file path
#' @param control control file path
#' @param outdir output directory
#' @param ... additional flags defined by dreme. Pass flag and resulting value
#'   (if any) as arg = value. If the flag has no value, set arg = TRUE.
#'
#' @return vector of flags for system2 or processx
#'
#' @importFrom magrittr %>%
#'
#' @noRd
prepareDremeFlags <- function(input, control, outdir, ...){
  argDict <- c(nmotifs = "m",
               sec = "t",
               evalue = "e",
               seed = "s",
               input = "p",
               control = "n",
               outdir = "oc",
               ngen = "g")

  flags <- cmdlr::getAllArgs() %>%
    cmdlr::cmd_args_to_flags(argDict) %>%
    cmdlr::cmd_drop_flags(c("n" = "shuffle")) %>%
    cmdlr::cmd_crystallize_flags()

  return(flags)

}

#' Import dreme output to R
#'
#'
#' @param xml path to dreme_out/dreme.xml
#'
#' @return data.frame with `motif` column containing universalmotif object representation of each DREME motif.
#'
#' parseDreme("dreme_out/dreme.xml")
#' @noRd
parseDreme <- function(xml){
  dreme_stats <- dreme_motif_stats(xml)

  pfms <- dreme_to_pfm(xml)

  dreme_stats$motif <- pfms
  return(dreme_stats)
}

#' Returns Dreme help lines
#'
#' @param command path to ame. output of handle_meme_path(util = "ame")
#'
#' @return
#'
#' @noRd
dreme_help <- function(command){
  processx::run(command, "-h", error_on_status = FALSE)$stderr
}

#' Return statistics of DREME results
#'
#' @param dreme_xml_path path to dreme.xml
#'
#' @return data.frame of dreme result statistics
#'
#' Columns:
#'  - id: motif id. Number represents order of discovery.
#'  - alt: alternative id. Number represents order of discovery.
#'  - seq: IUPAC sequence of matched motif
#'  - length: basepair length of discovered motif
#'  - nsites: number of times motif is discovered in reference sequence (can be more than 1 per entry)
#'  - positive_hits: number of positive sequences with at least 1 site
#'  - negative_hits: number of negative sequences with at least 1 site
#'  - pvalue: pvalue
#'  - evalue: evalue
#'  - unerased_evalue
#'  - positive_total: total number of sequences in positives (ie number of fasta entries)
#'  - negative_total: total number of sequences in negatives
#'  - pos_frac/neg_frac: fraction of positive or negative sites with a match
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @noRd
dreme_motif_stats <- function(dreme_xml_path) {
  # Extract metadata about number of matches for each motif, etc.

  dreme_xml <- xml2::read_xml(dreme_xml_path)

  # Get information about # of positive & # negative regions
  pos_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_find_all("//positives") %>%
    attrs_to_df() %>%
    dplyr::mutate(count = as.numeric(as.character(count)))

  neg_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_find_all("//negatives") %>%
    attrs_to_df() %>%
    dplyr::mutate(count = as.numeric(as.character(count)))

  # Extract statistics & motif counts for each dreme motif
  motif_stats <- xml2::xml_children(dreme_xml)[2] %>%
    xml2::xml_children() %>%
    attrs_to_df()
  dbl_cols <- grep("[^id|alt|seq]", names(motif_stats), value = T)
  motif_stats %<>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate_at(dbl_cols, as.numeric) %>%
    dplyr::mutate_at(c("length", "nsites",
                       "p", "n"), as.integer)

  # append info about positive / negative regions
  # compute some useful statistics
  motif_stats_final <- motif_stats %>%
    dplyr::rename("positive_hits" = "p",
                  "negative_hits" = "n") %>%
    dplyr::mutate("positive_total" = pos_info$count %>% as.integer,
                  "negative_total" = neg_info$count %>% as.integer,
                  "pos_frac" = positive_hits/positive_total,
                  "neg_frac" = negative_hits/negative_total) %>%
    dplyr::mutate(rank = gsub("^m", "", id) %>% as.integer(),
                  id = paste0(id, "_", seq)) %>%
    dplyr::select(rank, dplyr::everything()) %>%
    # Finally, change id and alt to "name" and "altname"
    # for compatibility with universalmotif
    dplyr::rename("name" = "id",
                  "altname" = "alt")

  return(motif_stats_final)
}


#' Return named vector of letter frequencies from DREME run
#'
#' @param dreme_run_info `xml_nodeset` of the dreme run info
#'
#' @return named vector of letter frequencies, suitable as input to `bkg` in
#'   universalmotif::create_motif
#'
#' @importFrom magrittr %>%
#'
#' @noRd
dreme_get_background_freq <- function(dreme_run_info){
  background_entry <- xml2::xml_find_all(dreme_run_info, "//background")

  background_df <- background_entry %>%
    attrs_to_df() %>%
    dplyr::select(-from) %>%
    lapply(function(x) as.character(x) %>% as.numeric) %>%
    dplyr::bind_cols()

  background_freq <- as.numeric(background_df) %>%
    magrittr::set_names(names(background_df))

  return(background_freq)
}

#' return list of universalmotif PCM objects
#'
#' @param dreme_xml_path path to dreme.xml output
#'
#' @return list of PCM output in Universalmotif format. Appended with metadata
#'   from DREME output.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
dreme_to_pfm <- function(dreme_xml_path){
  dreme_xml <- xml2::read_xml(dreme_xml_path)

  dreme_run_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_children()

  background_freq <- dreme_get_background_freq(dreme_run_info)

  motifs_matrix <- xml2::xml_children(dreme_xml)[2] %>%
    xml2::xml_children() %>%
    purrr::map(get_probability_matrix) %>%
    purrr::map(t)

  motif_stats_list <- dreme_motif_stats(dreme_xml_path) %>%
    split(.$name)

  pfmList <- purrr::map2(motif_stats_list, motifs_matrix, ~{
    universalmotif::create_motif(.y,
                                 type = "PCM",
                                 name = .x$name,
                                 altname = .x$altname,
                                 bkg = background_freq,
                                 pval = .x$pvalue,
                                 nsites = .x$nsites,
                                 bkgsites = .x$negative_total,
                                 eval = .x$evalue)

  })

  return(pfmList)
}

#' Return probability matrix for dreme motif
#'
#' Takes a <motif></motif> XML entry to return the probability matrix
#'
#' @param motif_xml_entry
#'
#' @return position probability matrix
#'
#' @noRd
get_probability_matrix <- function(motif_xml_entry){
  # takes a <motif></motif> XML entry to return the probability matrix
  # WARNING: matrix is a character matrix (NOT NUMERIC)
  # need to do the lapply(matrix, function(x)
  # as.character(x) %>% as.numeric()) %>% bind_cols(.)
  # trick for numeric matrix
  motif_attr <- attrs_to_df(motif_xml_entry, stringsAsFactors = F)

  nsites <- motif_attr$length %>%
     as.character() %>%
     as.integer()

  freqs <- motif_xml_entry %>%
    xml2::xml_children(.) %>%
    .[1:nsites]

  freq_table <- lapply(freqs, attrs_to_df, stringsAsFactors = F) %>%
    dplyr::bind_rows()

  freq_matrix <- lapply(freq_table, function(x) as.character(x) %>% as.numeric) %>%
                    dplyr::bind_rows(.) %>%
                    as.matrix(.)

  return(freq_matrix)
}

#' Return line reporting number of motifs passing cutoff in DREME stdout
#'
#' @param stdout stdout from processx
#'
#' @return
#'
#' @noRd
dreme_nmotifs_line <- function(stdout){
  lines <- strsplit(stdout, "\n") %>%
    .[[1]]

  matchLine <- grep("\\d motifs with E-value <", lines, value = T)
  return(matchLine)
}

#' Grab number of discovered motifs from stdout line
#'
#' @param line output of dreme_nmotifs_line
#'
#' @return
#'
#' @noRd
dreme_nmotifs <- function(line){
  nmotifs <- gsub("(^\\d+).+", "\\1", line)
  return(as.integer(nmotifs))
}

#' Return number of discovered DREME motifs
#'
#' @param processx_out output of processx run
#'
#' @return `integer(1)` of number of motifs passing threshold
#' @export
#'
#' @noRd
dreme_nmotifs_found <- function(processx_out){
  processx_out$stdout %>%
    dreme_nmotifs_line() %>%
    dreme_nmotifs()
}
