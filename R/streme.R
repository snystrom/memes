#' Denovo motif discovery of target regions using STREME
#'
#' STREME discovers short, ungapped, *de-novo* motifs that are enriched or
#' relatively enriched relative to a control set of sequences. STREME can be run
#' to discover motifs relative to a shuffled set of input sequences, against a
#' separately provided set of "control" sequences, or to determine whether
#' motifs are centrally enriched within input sequences.
#'
#' Properly setting the `control` parameter is key to discovering biologically
#' relevant motifs. Often, using `control = "shuffle"` will produce a suboptimal
#' set of motifs; however, some discriminative analysis designs don't have
#' proper "control" regions other than to shuffle.
#' 
#'
#' If you have fewer than 50 sequences, consider using [runMeme()] instead.
#'
#' @param input regions to scan for motifs. If using `objfun = "cd"` to test for
#'   centrally enriched motifs, be sure to include sufficient flanking sequence
#'   (e.g. +/- 500bp) for an accurate estimate. Can be any of:
#'   - path to fasta file
#'   - DNAStringSet object (can be generated from GRanges using `get_sequence()`)
#'   - List of DNAStringSet objects (generated from `get_sequence()`)
#'   - *NOTE:* if using StringSet inputs, each entry must be named (set with `names()`).
#'   - *NOTE:* If you want to retain the raw streme output files, you must use a
#'   path to fasta file as input, or specify an "outdir"
#' @param control regions to use as background for motif search. These should
#'   have a similar length distribution as the input sequences. Can be any of:
#'   - path to fasta file
#'   - DNAStringSet object (can be generated from GRanges using get_sequence)
#'   - A Biostrings::BStringSetList (generated using `get_sequence`), in which
#'   case all sequences in the list will be combined as the control set.
#'   - if `input` is a list of DNAStringSet objects, a character vector of names
#'   in `input` will use those sequences as background. runstreme will not scan
#'   those regions as input.
#'   - "shuffle" to use streme's built-in dinucleotide shuffle feature (NOTE: if
#'   `input` is a list object with an entry named "shuffle", the list entry will
#'   be used instead).
#'   Optionally can also pass `seed = <any number>` to `...` to use as the random
#'   seed during shuffling. If no seed is passed, streme will use 0 as the random
#'   seed, so results will be reproducible if rerunning. 
#' @param outdir path to output directory of streme files, or "auto" to
#'   autogenerate path. Default: location of input fasta in dir named
#'   "\<input\>_vs_\<control\>". If input is DNAstringset, will be temporary
#'   path. This means that if you want to save the raw output files, you must
#'   use fasta files as input or use an informative (and unique) outdir name.
#'   memes will **not check** if it overwrites files in a directory. Directories
#'   will be recursively created if needed. (default: "auto")
#' @param objfun one of c("de", "cd"). Default: "de" for differential
#'   enrichment. "cd" for central distance (control must be set to NA for "cd").
#' @param alph one of c("dna", "rna", "protein") or a path to a MEME format alph file. (default: "dna")
#' @param meme_path path to "meme/bin"
#' @param silent Whether to suppress printing stdout & stderr to console
#'   (default: TRUE). Warnings are always printed regardless of this setting.
#' @param ... pass any commandline options as R function arguments. For a
#'   complete list of STREME options, see 
#'   [the STREME manual](https://meme-suite.org/meme/doc/streme.html).
#'
#' @return a `universalmotif_df` of STREME Motifs
#' @seealso `?universalmotif::tidy-motifs`
#' @details 
#' 
#'  # Citation
#'
#'  If you use `runStreme()` in your analysis, please cite:
#'
#'  Timothy L. Bailey, "STREME: Accurate and versatile sequence motif
#'  discovery", Bioinformatics, 2021.
#'  https://doi.org/10.1093/bioinformatics/btab203
#'  
#'  # Licensing
#'  The MEME Suite is free for non-profit use, but for-profit users should purchase a
#'  license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
#'
#' @export
#'
#' @examples
runStreme <- function(input, control, outdir = "auto", objfun = "de", 
                      alph = "dna", meme_path = NULL, silent = TRUE, ...){
  UseMethod("runStreme")
}

#' @export
#' @noRd
runStreme.list <- function(input, control, outdir = "auto", objfun = "de", 
                      alph = "dna", meme_path = NULL, silent = TRUE, ...){
  
  x <- sequence_input_control_list(input, control)
  input <- x$input
  control <- x$control

  res <- purrr::map(input, runStreme.default,
             control = control,
             outdir = outdir,
             objfun = objfun,
             alph = alph,
             meme_path = meme_path,
             silent = silent,
             ...
             )
  return(res)
}

#' @export
#' @noRd
runStreme.BStringSetList <- function(input, control, outdir = "auto", objfun = "de", 
                      alph = "dna", meme_path = NULL, silent = TRUE, ...){
  runStreme.list(as.list(input), control, outdir, objfun, alph, meme_path, silent, ...)
}

#' @export
#' @noRd
runStreme.BStringSet <- function(input, control, outdir = "auto", objfun = "de", 
                      alph = "dna", meme_path = NULL, silent = TRUE, ...){
  runStreme.default(input, control, outdir, objfun, alph, meme_path, silent, ...)
}

#' @export
#' @noRd
runStreme.default <- function(input, control, outdir = "auto", objfun = "de", 
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
  process_check_error(ps_out, 
                      help_fun = ~{streme_help_flags(command)}, 
                      user_flags = cmdfun::cmd_help_parse_flags(flags), 
                      default_help_fun = FALSE)
  if (!silent){
    message(ps_out$stdout)
    message(ps_out$stderr)
  }
  
  #TODO: return NULL if no motifs found?
  # - Not sure here, seems like streme will always find motifs
  # - unless patience = 0
  
  print_streme_messages(ps_out$stderr)
  
  streme_out <- cmdfun::cmd_file_expect("streme", c("txt", "html", "xml"), 
                                        outdir = outdir)

  streme_results <- parseStreme(streme_out$xml)

  return(streme_results)
}

#' Prepare flags vector input
#'
#' @param input 
#' @param control 
#' @param outdir 
#' @param alph 
#' @param ... 
#'
#' @return
#' @noRd
prepareStremeFlags <- function(input, control, outdir, alph, ...){
  argsDict <- c("input" = "p",
               "control" = "n",
               "outdir" = "oc")
  
  # handle alphabet assignment
  # allows setting alph = "dna" instead of dna = TRUE
  alph_flags <- meme_alph_to_args(alph) %>%
    cmdfun::cmd_list_interp()

  flagsList <- cmdfun::cmd_args_all(drop = "alph") %>%
    cmdfun::cmd_list_interp(argsDict) %>% 
    cmdfun::cmd_list_drop(c("n" = "shuffle"))

  flags <- c(flagsList, alph_flags) %>%
    cmdfun::cmd_list_to_flags(prefix = "--")
  return(flags)
}

#' return streme help text
#'
#' @param command path to streme
#'
#' @return
#' @noRd
streme_help <- function(command){
  processx::run(command, "--help", error_on_status = FALSE)$stderr
}

#' return vector of streme flag names
#'
#' @param command path to streme
#'
#' @return
#' @noRd
streme_help_flags <- function(command){
  streme_help(command) %>% 
    {strsplit(., "\n")[[1]]} %>% 
    gsub("\t", " ", .) %>% 
    gsub("\\[", "", .) %>%
    gsub("\\]", "", .) %>% 
    cmdfun::cmd_help_parse_flags(split_newline = FALSE) %>% 
    gsub(";", "", .) %>% 
    unique()
}

#' Print warnings from Streme
#'
#' @param stderr stderr text
#'
#' @return
#' @noRd
print_streme_messages <- function(stderr){
  stderr %>% 
    strsplit(., "\n") %>% 
    .[[1]] %>% 
    {
      x <- .
      grep("Warning|would have had fewer", x, value = TRUE) %>% 
        gsub("^# ", "", .) %>% 
        gsub("$", "\n", .) %>% 
        purrr::walk(message)
    }
}
#' Import streme results from xml
#'
#' @param xml path to streme xml
#'
#' @return universalmotif_df
#' @noRd
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
  streme_version <- xml2::xml_attr(streme_xml, "version")

  streme_run_info <- xml2::xml_children(streme_xml)[1] %>%
    xml2::xml_children()

  background_freq <- streme_get_background_freq(streme_run_info, streme_version)

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

#' Return named vector of letter frequencies from STREME run
#'
#' @param streme_run_info `xml_nodeset` of the dreme run info
#' @param version STREME version number
#'
#' @return named vector of letter frequencies, suitable as input to `bkg` in
#'   universalmotif::create_motif
#'
#' @importFrom magrittr %>%
#' @importFrom utils compareVersion
#'
#' @noRd
streme_get_background_freq <- function(streme_run_info, version){

  # STREME changed how background frequencies were encoded in XML around 5.4.1
  # unsure when the change was made because it is undocumented in release notes.
  if (compareVersion(version, "5.4.1") >= 0) {
    return(streme_get_background_freq_new(streme_run_info))
  } else {
    return(streme_get_background_freq_old(streme_run_info))
  }
}

#' Return named vector of letter frequencies from STREME run
#'
#' @param streme_run_info `xml_nodeset` of the dreme run info
#'
#' @return named vector of letter frequencies, suitable as input to `bkg` in
#'   universalmotif::create_motif
#'
#' @importFrom magrittr %>%
#'
#' @noRd
streme_get_background_freq_new <- function(streme_run_info){

  background_entry <- xml2::xml_find_all(streme_run_info, "//background_frequencies")

  alphabet_array <- background_entry %>%
    xml2::xml_find_all("//alphabet_array") %>%
    xml2::xml_children()

  background_freq <- xml2::xml_double(alphabet_array) %>%
    magrittr::set_names(xml2::xml_attrs(alphabet_array) %>% unlist())

  return(background_freq)
}

#' Return named vector of letter frequencies from STREME run
#'
#' @param streme_run_info `xml_nodeset` of the dreme run info
#'
#' @return named vector of letter frequencies, suitable as input to `bkg` in
#'   universalmotif::create_motif
#'
#' @importFrom magrittr %>%
#'
#' @noRd
streme_get_background_freq_old <- function(streme_run_info){
  background_entry <- xml2::xml_find_all(streme_run_info, "//background")

  background_df <- background_entry %>%
    attrs_to_df() %>%
    dplyr::select(-"from") %>%
    lapply(function(x) as.character(x) %>% as.numeric) %>%
    dplyr::bind_cols()

  background_freq <- as.numeric(background_df) %>%
    magrittr::set_names(names(background_df))

  return(background_freq)
}


#' Parse motif IDs from streme output
#' 
#' 0-pad numeric ids to preserve sort order, add 'm' prefix to match other MEME
#' Suite tool naming convention
#' 
#' @param ids the raw ids from the STREME XML (e.g. "1-CACAC")
#' @return A vector of cleaned ids (e.g. "m01-CACAC")
#' 
#' @noRd
clean_streme_ids <- function(ids) {
  
  raw_numeric_ids <- strsplit(ids, "-", fixed = TRUE) %>% 
    purrr::map_chr(head, 1)
  
  # The max number of 0's to pad to each id 
  numeric_id_lengths <- nchar(raw_numeric_ids)
  
  # The minimum length of a cleaned numeric ID after 0 padding
  # e.g. 2 for "1" => "01"
  min_numeric_id_length <- 2
  
  # Pad so that all numeric IDs are equal length (match length of longest)
  if (max(numeric_id_lengths) > min_numeric_id_length) {
    min_numeric_id_length <- max(numeric_id_lengths)
  }
  
  zeros_to_add <- min_numeric_id_length - numeric_id_lengths
 
  # Pad zeros to correct length 
  purrr::map2_chr(ids, zeros_to_add, ~{
    # e.x. 
    # 1-ATGC => m01-ATGC
    paste0("m", rep(0, length.out = .y), .x)
  }) %>% 
    gsub("-", "_", .)
  
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
    dplyr::mutate("id" = clean_streme_ids(id)) %>% 
    dplyr::rename("name" = "id",
                  "altname" = "alt") %>% 
    dplyr::mutate_at(int_cols, as.integer) %>% 
    dplyr::mutate_at(dbl_cols, as.double) %>% 
    dplyr::mutate(nsites = test_pos_count + train_pos_count) %>% 
    # Put test variables first (do before rename)
    dplyr::select("name", "altname", "width", 
                  "initial_width", "seed", "nsites", "score_threshold",
                  dplyr::starts_with("test_"), 
                  dplyr::starts_with("train_"), 
                  dplyr::everything()) %>% 
    # Change the "test_" columns for compatibility w/ universalmotif_df
    dplyr::rename_with(~{gsub("test_", "", .)},
                     dplyr::starts_with("test_"), 
                     ) %>% 
    # pvalue -> pval for universalmotif compatability
    dplyr::rename_with(~{gsub("pvalue", "pval", .)},
                     dplyr::contains("pvalue"), 
                     )
}
