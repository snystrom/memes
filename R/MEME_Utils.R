
#' Run DREME from fasta files
#'
#' @param input
#' @param control
#' @param evalue
#' @param n_motifs
#' @param sec
#' @param outdir
#' @param meme_path
#' @param mink
#' @param maxk
#' @param ngen number of regular expressions to generalize
#' @param alph one of c("dna", "rna", "protein") or a file path to custom alphabet file. Note: "protein" is not recommended by DREME authors.
#'
#' @return
#' @export
#'
#' @examples
runDreme <- function(input, control, evalue = 0.05, n_motifs = 10, sec = 18000, alph = c("dna", "rna", "protein"),
                     ngen = 100,
                     outdir = paste0(dirname(input), "/", basename(tools::file_path_sans_ext(input)),  "_vs_", basename(tools::file_path_sans_ext(control))),
                     meme_path = Sys.getenv("MEME_DIR"),
                     mink = 3,
                     maxk = 8){
  # TODO: warn if alph = "protein" | double check MEME documentation for valid inputs
  # # NOTE: currently does not support custom alphabet file
  # NOTE: control if shuffle will not currently use dreme shuffling algorithm. add this functionality.

  if (file.exists(alph)) {
    alph <- glue::glue("-alph {alph}")
  } else {
    alph <- tolower(alph)
    assertthat::assert_that(alph %in% c("dna", "rna", "protein"),
                            msg = glue::glue("alph must be one of 'dna', 'rna', 'protein', or a path to custom alphabet."))

    if (alph == "protein") warning("Setting `alph` to `protein` is not recommended by DREME authors")

    alph <- glue::glue("-{alph}")

  }

  # Requires setting Sys.setenv(MEME_DIR = "path/to/meme")
  dreme_command <- glue::glue("{meme_path}/bin/dreme {alph} -oc {outdir} \\
                              -g {ngen} \\
                              -p {input} -n {control} -t {sec} \\
                              -e {evalue} -m {n_motifs} \\
                              -mink {mink} -maxk {maxk}")

  print(dreme_command)

  system(dreme_command)

  dreme_out <- list(
    html = paste0(outdir, "/dreme.html"),
    xml = paste0(outdir, "/dreme.xml"),
    txt = paste0(outdir, "/dreme.txt"))

  return(dreme_out)
}

#' Title
#'
#' @param input
#' @param database
#' @param min_overlap
#' @param dist
#' @param thresh
#' @param outdir
#' @param meme_path
#'
#' @return
#' @export
#'
#' @examples
runTomTom <- function(input, database = Sys.getenv("TOMTOM_DB"), min_overlap = 5, dist = "pearson", thresh = 10.0,
                      outdir = paste0(dirname(input), "/tomtom"),
                      meme_path = Sys.getenv("MEME_DIR")){
  # Optional set Sys.setenv(TOMTOM_DB = "path/to/db.meme")
  tomtom_command <- glue::glue("{meme_path}/bin/tomtom -no-ssc -oc {outdir} \\
                               -min-overlap {min_overlap} \\
                               -dist {dist} \\
                               -evalue \\
                               -thresh {thresh} \\
                               {input} \\
                               {database}")
  system(tomtom_command)

  tomtom_out <- list(
    html = paste0(outdir, "/tomtom.html"),
    xml = paste0(outdir, "/tomtom.xml"),
    txt = paste0(outdir, "/tomtom.tsv"))

  return(tomtom_out)
}

#' Title
#'
#' @param input
#' @param control
#' @param outdir
#' @param scoring
#' @param method
#' @param bgformat
#' @param pvalue_threshold
#' @param verbose
#' @param database
#' @param meme_path
#'
#' @return
#' @export
#'
#' @examples
runAme <- function(input, control = "shuffle",
                   outdir = paste0(dirname(input), "/", basename(input), "_vs_", basename(control), "/ame"),
                   #inc = "\\*",
                   scoring = "avg", method = "ranksum", bgformat = 1, pvalue_threshold = 0.05, verbose = 1,
                   database = Sys.getenv("TOMTOM_DB"),
                   meme_path = Sys.getenv("MEME_DIR")){

  # inc = regex matching motif names. By default will use all motifs, but regex will restrict it to only matches.
  #       Note: must escape wildcards so not evaluated by shell

  dir.create(outdir, recursive = T)

  if(control == "shuffle") {
    print("Using shuffled input as background")
    shuffle_fa <- shuffle_fasta_seq(input)

    shuffle_fa_path <- glue::glue("{input}.shuffle")
    Biostrings::writeXStringSet(shuffle_fa, filepath = shuffle_fa_path, format = "fasta")
    control <- shuffle_fa_path
  }

  #ame_command <- glue::glue("{meme_path}/bin/ame --verbose {verbose} --bgformat {bgformat} \\
  #                          --inc {inc} \\
  #                          --scoring {scoring} --method {method} \\
  #                          --pvalue-report-threshold {pvalue_threshold} \\
  #                          --o {outdir} --control {control} {input} {database}")
  ame_command <- glue::glue("{meme_path}/bin/ame --verbose {verbose} --bgformat {bgformat} \\
                            --scoring {scoring} --method {method} \\
                            --pvalue-report-threshold {pvalue_threshold} \\
                            --o {outdir} --control {control} {input} {database}")

  system(ame_command)

  ame_out <- list(
    html = paste0(outdir, "/ame.html"),
    txt = paste0(outdir, "/ame.txt"))

  return(ame_out)
}

#' Title
#'
#' @param motifs
#' @param genome
#' @param background
#' @param thresh
#' @param meme_path
#'
#' @return
#' @export
#'
#' @examples
runFimoGenome <- function(motifs, genome, background, thresh = 0.01,
                   meme_path = Sys.getenv("MEME_DIR")){

  # motifs = path to meme format motif file
  # genome = path to genome fasta (ie dm3.fa)
  # background = path to background frequency file generated using meme tools

  # Runs fimo on whole genome, imports all matches to motif & computes qvalue for each with benjamini-hochberg correction.
  # Takes meme format file as input.

  # create temporary output file to store results on disk
  timestamp <- Sys.time() %>% gsub(" ", "_", .)
  outfile <- glue::glue("fimo_out_{timestamp}.tmp")

  fimo_command <- glue::glue("{meme_path}/bin/fimo --bgfile {background} --thresh {thresh} \\
                             --max-strand \\
                             --text --skip-matched-sequence --verbosity 2 \\
                             {motifs} {genome} > {outfile}")

  system(fimo_command)

  fimo_matches <- readr::read_tsv(outfile) %>%
    dplyr::rename_all(function(x) gsub("-", ".", x)) %>%
    dplyr::rename_all(function(x) gsub("#", "", x)) %>%
    dplyr::rename_all(function(x) gsub(" ", "", x)) %>%
    dplyr::rename(seqnames = sequence_name) %>%
    dplyr::group_by(motif_alt_id) %>%
    dplyr::mutate(q.value = p.value/n())

  unlink(outfile)

  return(fimo_matches)

}


#' Title
#'
#' @param xml
#'
#' @return
#'
#' @examples
attrs_to_df <- function(xml) {
  # parsing DREME output
  # converts xml attributes to dataframe
  # where each column is an attribute
  xml2::xml_attrs(xml) %>%
    data.frame() %>%
    t() %>%
    data.frame(., row.names = NULL)
}

get_probability_matrix <- function(motif_xml_entry){
  # takes a <motif></motif> XML entry to return the probability matrix
  # WARNING: matrix is a character matrix (NOT NUMERIC)
  # need to do the lapply(matrix, function(x) as.character(x) %>% as.numeric()) %>% bind_cols(.) trick for num
  # this should be an internal function
  motif_attr <- attrs_to_df(motif_xml_entry)

  nsites <- motif_attr$length %>%
     as.character() %>%
     as.integer()

  freqs <- motif_xml_entry %>%
    xml2::xml_children(.) %>%
    .[1:nsites]

  freq_table <- lapply(freqs, attrs_to_df) %>%
    dplyr::bind_rows()

  return(freq_table)
}

dreme_to_pfm <- function(dreme_xml_path){
  # input: dreme xml output
  # output: pfm matrix list
  # NOTE: see `write_dreme_motifs` for comments on things
  # work in progress... maybe just write your own object type....
  dreme_xml <- xml2::read_xml(dreme_xml_path)

  dreme_attr <- attrs_to_df(dreme_xml)

  dreme_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_children()

  background_entry <- xml2::xml_find_all(dreme_info, "//background")
  bg <- background_entry %>%
    attrs_to_df() %>%
    dplyr::select(-from) %>%
    lapply(., function(x) as.character(x) %>% as.numeric) %>%
    dplyr::bind_cols(.)

  bg_params <- c(A=bg$A, C=bg$C, G=bg$G, T=bg$T)

  motifs <- xml2::xml_children(dreme_xml)[2] %>%
    xml2::xml_children()

  # parse & write-out each motif entry to MEME format
  pfmList <- lapply(motifs, function(motif_entry){
    motif_attr <- attrs_to_df(motif_entry)

    # parse filename:
    #outfile <- paste0(outdir, motif_attr$id, "_", motif_attr$alt, "_", motif_attr$seq, ".meme")


    # write probability matrix info pfm
    prob_info_line <- paste0("letter-probability matrix: alength= 4 w= ", motif_attr$length,
           " nsites= ", motif_attr$nsites ,
           " E= ", motif_attr$evalue)
    # write probability matrix to pwm output:
    freq_table <- get_probability_matrix(motif_entry)
    freq_matrix <- lapply(freq_table, function(x) as.character(x) %>% as.numeric) %>%
                      dplyr::bind_rows(.) %>%
                      as.matrix(.)
    freq_matrix
    #pfm <- TFBSTools::PFMatrix(ID=as.character(motif_attr$id),
    #                           name=as.character(motif_attr$alt),
    #                           tags = list(type = "DREME"),
    #                           strand = "+",
    #                           bg=bg_params,
    #                           profileMatrix=t(freq_matrix))
  })
  return(pfmList)
}

write_dreme_motifs <- function(dreme_xml_path, outdir = gsub(basename(dreme_xml_path), "", dreme_xml_path), name_type = "seq") {
  # takes dreme.xml file as input
  # will write out *.meme format motif files for each de-novo motif called

  # where name_type is one of c("seq", "id", "alt"), such that
  # seq = MOTIF AACA
  # and alt = MOTIF DREME-1
  # and id = MOTIF m01
  # in other words, seq will name the motif as its consensus sequence, while id will name is by its dreme ranking
  # NOTE: only works for DNA output. Will set strand = both (+ -), though this could
  # theoretically be parsed, I need to figure out a better way to do it.
  # DOES NOT work with custom alphabet and will not parse it (forces ACGT).
  # ASSUMES DNA and strand = both and will not check otherwise
  # TODO:
    # could have it parse the source in "background_line1" but for now just says "from DREME"
    # need more rigorous testing perhaps of parameters/validity, but this should do for now...
    # check that name_type is either c("seq", "id", "alt")

  if (stringr::str_sub(outdir, -1) != "/") {outdir <- paste0(outdir, "/")}
  R.utils::mkdirs(outdir)

  dreme_xml <- xml2::read_xml(dreme_xml_path)

  dreme_attr <- attrs_to_df(dreme_xml)

  dreme_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_children()

  # for future reference, these are the positions of relevant data types in the xml info:
  #alphabet <- dreme_info[4]
  #strand <- dreme_info[5]
  #background_entry <- dreme_info[6]

  # parse background frequency
  # NOTE: this could be a source of breakage if there are multiple entries containing "background".
  # I'm pretty positive this won't grep any commandline arguments--should only grab xml <background> tag.
  # I opted for this instead of the dreme_info[6] positional argument calling
  # because I think this is more future-proof but I could be wrong... If background entry stuff breaks later, try this first.
  background_entry <- xml2::xml_find_all(dreme_info, "//background")
  bg <- background_entry %>%
    attrs_to_df() %>%
    dplyr::select(-from)

  # write MEME-format motif file header:
  meme_version <- paste0("MEME version ", dreme_attr$version, "\n")
  alphabet <- "ALPHABET= ACGT\n"
  strands <- "strands: + -\n" # not sure how to parse otherwise. field is optional, could omit.
  background_line1 <- "Background letter frequencies (from DREME):"
  background_line2 <- paste0(paste("A", bg$A, "C", bg$C, "G", bg$G, "T", bg$T), "\n")

  meme_header <- c(meme_version, alphabet, strands, background_line1, background_line2)

  # list of each motif entry
  motifs <- xml2::xml_children(dreme_xml)[2] %>%
    xml2::xml_children()

  # parse & write-out each motif entry to MEME format
  outfiles <- lapply(motifs, function(motif_entry){
    motif_attr <- attrs_to_df(motif_entry)
    # parse filename:
    outfile <- paste0(outdir, motif_attr$id, "_", motif_attr$alt, "_", motif_attr$seq, ".meme")
    # write motif line
    motif_line <- switch(name_type,
                         seq = paste0("MOTIF ", motif_attr$seq, "\n"),
                         id = paste0("MOTIF ", motif_attr$id, "\n"),
                         alt = paste0("MOTIF ", motif_attr$alt, "\n"))

    # write probability matrix info to MEME output file
    prob_info_line <- paste0("letter-probability matrix: alength= 4 w= ", motif_attr$length,
           " nsites= ", motif_attr$nsites ,
           " E= ", motif_attr$evalue)
    write(c(meme_header, motif_line, prob_info_line), file = outfile)
    # write probability matrix to output:
    freq_table <- get_probability_matrix(motif_entry)
    freq_table$A <- paste0("  ", freq_table$A)

    write.table(freq_table, file = outfile, append = T,
                quote = FALSE,
                sep = "  \t",
                row.names = FALSE, col.names = FALSE)

    # add trailing newline: very important for parsing to some programs (*sigh*)
    # worst file format ever
    write("\n", file = outfile, append = T)
    return(outfile)
  })
  return(unlist(outfiles))
}

dreme_motif_stats <- function(dreme_xml_path, ...) {
  # Extract metadata about number of matches for each motif, etc.

  dreme_xml <- xml2::read_xml(dreme_xml_path)
  dreme_attrs <- attrs_to_df(dreme_xml)

  # Get information about # of positive & # negative regions
  pos_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_find_all(., "//positives") %>%
    attrs_to_df() %>%
    dplyr::mutate(count = as.numeric(as.character(count)))

  neg_info <- xml2::xml_children(dreme_xml)[1] %>%
    xml2::xml_find_all(., "//negatives") %>%
    attrs_to_df() %>%
    dplyr::mutate(count = as.numeric(as.character(count)))

  # Extract statistics & motif counts for each dreme motif
  motif_stats <- xml2::xml_children(dreme_xml)[2] %>%
    xml2::xml_children() %>%
    attrs_to_df(.)
  dbl_cols <- grep("[^id|alt|seq]", names(motif_stats), value = T)
  motif_stats %<>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate_at(dbl_cols, as.numeric)

  # append info about positive / negative regions
  # compute some useful statistics
  motif_stats_final <- motif_stats %>%
    dplyr::rename(positive_hits = p,
                  negative_hits = n) %>%
    dplyr::mutate(positive_total = pos_info$count,
                  negative_total = neg_info$count,
                  pos_frac = positive_hits/positive_total,
                  neg_frac = negative_hits/negative_total)

  return(motif_stats_final)
}


centrimoSiteCounts <- function(site_count_path){
  siteCounts_tsv <- readr::read_tsv(site_count_path, col_names = F)

  # extract each motif name & info by grabbing text entries
  grepString <- siteCounts_tsv[1,1] # not sure if best approach if string is not always the same across samples (could build regex to drop numbers)

  (motifInfo <- siteCounts_tsv %>%
    dplyr::filter(grepl(grepString, X1)) %>%
    dplyr::mutate(rownum = grep(grepString, siteCounts_tsv$X1)) %>%
    dplyr::rename(id = "X2", alt = "X3"))

  # for each motif
  # get position-specific info from range k:j
  # where i is rownum of first motif + 1
  # and j is rownum of second motif - 1

  # make motif_range as vector, map to motifInfo
  # this will extract the position-specific data while excluding the header info
  pos_start <- motifInfo$rownum + 1
  pos_end <- motifInfo$rownum - 1
  pos_end <- c(pos_end[-1], nrow(siteCounts_tsv))

  motifInfo$motif_range <- map(1:nrow(motifInfo), ~{pos_start[.]:pos_end[.]})

  # extract annotated centrimo output:
  centrimo_out <- motifInfo %>%
    dplyr::mutate(info = map(motif_range, ~{siteCounts_tsv[.,]})) %>%
    dplyr::mutate(info = map(info, ~{
      # need to rename columns to reasonable names, and convert to integers
      df <- .
      df %<>%
        dplyr::rename(position = "X1", matches = "X2", negative = "X3")

      for (n in names(df)) {
        df[[n]] %<>% as.integer()
      }

      return(df)
    })) %>%
    tidyr::unnest(info, .drop = T) %>%
    dplyr::select(-rownum)
  return(centrimo_out)
}

centrimo_txt_parse <- function(centrimo_txt_path){
  col_names <- c("db", "id_alt", "consensus", "E.value", "p.adj", "log.p.adj",
                 "bin.location", "bin.width", "total.width", "sites.in.bin", "total.sites",
                 "p.success", "p.value", "mult.tests", "neg.sites.in.bin", "neg.sites",
                 "neg.p.adj", "neg.log.p.adj", "DROP", "fisher.p.adj", "fisher.log.p.adj")
  #col_names <- readr::read_lines(centrimo_txt_path, skip = 1, n_max = 1) %>%
  #  gsub("# ", "", .) %>%
  #  stringr::str_split(., "\t") %>%
  #  unlist
  info <- readr::read_tsv(centrimo_txt_path, comment = "#", col_names = col_names) %>%
    tidyr::separate(., id_alt, c("id", "alt"), sep = "( )+") %>%
    dplyr::select(-DROP) %>%
    data.frame
}


parseCentrimo <- function(centrimo_out_path){
  #TODO: check file exists
  centrimo_txt_info <- centrimo_txt_parse(paste0(centrimo_out_path, "centrimo.txt"))
  centrimo_counts <- centrimoSiteCounts(paste0(centrimo_out_path, "site_counts.txt"))

  centrimo_full_out <- centrimo_counts %>%
    dplyr::left_join(., centrimo_txt_info %>% dplyr::select(-alt), by = "id")

}

readTomTom <- function(txt){
  df <- readr::read_tsv(txt)

  names(df) <- names(df) %>%
    gsub("#", "", .) %>%
    gsub(" ", "_", .) %>%
    gsub("-", "_", .) %>%
    tolower

  return(df)
}


dremeStats <- function(dreme_xml, tomtom_txt){
  # Outputs info about each dreme match w/ corresponding best motif match
  dreme_stats <- dreme_motif_stats(dreme_xml)
  tomtom_stats <- readTomTom(tomtom_txt) %>%
    dplyr::group_by(query_id) %>%
    tidyr::nest(.key = "tomtom")

  dreme_stats_full <- dplyr::left_join(dreme_stats, tomtom_stats, by = c("seq" = "query_id")) %>%
  #dreme_stats_full <- dplyr::full_join(dreme_stats, tomtom_stats, by = c("seq" = "query_id")) %>%
    dplyr::mutate(bestMatch = purrr::map(tomtom, ~{
      df <- .

      if(is.null(df)) { return("None") } # in case no tomtom match

      df %>%
        dplyr::filter(dplyr::row_number(p_value) == 1) %>% # return first match w/ lowest p-value
        .$target_id

    })) %>%
    tidyr::unnest(bestMatch) %>%
    dplyr::select(id, alt, seq, bestMatch, dplyr::everything())

  return(dreme_stats_full)
}

dreme_motifCount <- function(dreme_path){
  # count number of motif entries in dreme output
  n_motifs <- system(glue::glue("grep -c MOTIF {dreme_path}"), intern = T) %>% as.integer
  return(n_motifs)
}


shuffle_fasta_seq <- function(fasta){
  # input: fasta path
  # output: shuffled input sequences
  fa <- Biostrings::readDNAStringSet(fasta)

  fa_split <- strsplit(as.character(fa), split = "", fixed = T)

  fa_shuffle <- purrr::map(fa_split, sample) %>%
    Biostrings::unstrsplit(.) %>%
    Biostrings::DNAStringSet(.)

  return(fa_shuffle)
}

readMotifAnalysisOutput <- function(dreme_files, tomtom_files){
  # Takes output of runDreme and runTomTom & creates motif_data list with stats & motifs as pwm
  # also adds column to dreme stats pointing to the meme format motif file which can be used as input to other meme-suite functions
  dreme_stats <- dremeStats(dreme_files$xml, tomtom_files$txt)

  dreme_motif_paths <- write_dreme_motifs(dreme_files$xml)

  dreme_motifs <- motifStack::importMatrix(dreme_motif_paths)

  motifPathDf <- tibble::tibble(meme_motif_path = dreme_motif_paths) %>%
    dplyr::mutate(id = basename(meme_motif_path) %>%
                    stringr::str_extract(., "m\\d+"))

  dreme_stats %<>%
    dplyr::left_join(., motifPathDf, by = "id")

  motif_data <- list(info = dreme_stats,
                     motifs = dreme_motifs)

  return(motif_data)
}

motifAnalysis <- function(input, control, database = Sys.getenv("TOMTOM_DB"), evalue = 0.05, n_motifs = 10, meme_path = Sys.getenv("MEME_DIR"),
                          tt_evalue = 10, ...){
  # Runs DREME & TOMTOM on dreme output for input vs control samples.
  # Takes path to fasta file as input.

  # TODO: Take list input for tomtom params & dreme params
  if(control == "shuffle") {
    print("Using shuffled input as background")
    shuffle_fa <- shuffle_fasta_seq(input)
    # TODO: use shuffle algorithm that preserves dinucleotide frequency:
    # https://stackoverflow.com/questions/26497583/split-a-string-every-5-characters
    # gsub("(.{2})", "\\1 ", sequence)

    shuffle_fa_path <- glue::glue("{input}.shuffle") %>% as.character
    Biostrings::writeXStringSet(shuffle_fa, filepath = shuffle_fa_path, format = "fasta")
    control <- shuffle_fa_path
    }

  dreme_files <- runDreme(input, control, evalue = evalue, n_motifs = n_motifs, meme_path = meme_path, ...)

  dreme_nMotifs <- dreme_motifCount(dreme_files$txt)

  if (dreme_nMotifs == 0){
    warning("No denovo motifs found")
    return(NULL)
  }

  tomtom_files <- runTomTom(dreme_files$txt, thresh = tt_evalue)

  motif_data <- readMotifAnalysisOutput(dreme_files, tomtom_files)

  return(motif_data)
}

find_id_match <- function(tomtom, tfid, best_id = NULL){
  # find a regex match to a motif within the $tomtom output, return the target_id of the match
  # if best_id is set, this will override the top tomtom match as the best output
  # (this is so forceBestMach can be piped together & acquire serial assignment of best match)
  # if multiple matches, will return the first entry (tomtom output is default sorted by lowest p-value, so first is most significant)
  ## TODO: return best = "first", "lowest_p", etc.

  tfid_match_row <- grep(tfid, tomtom$target_id)

  best_id <- ifelse(is.null(best_id), tomtom[1,1], best_id)

  best <- ifelse(length(tfid_match_row) != 0, tomtom[tfid_match_row,1], best_id) %>% unlist
  return(best[1])
}

forceBestMatch <- function(df, tfid){
  # take motifAnalysis$info output, search each list of tomtom hits for match to `tfid` string
  # replaces this match for 'bestMatch' column in .$info, returns .$info data frame

  df %<>%
    dplyr::mutate(bestMatch = map2_chr(tomtom, bestMatch, ~{find_id_match(.x, tfid, best_id = .y)}))
 return(df)
}


plot_denovo_matches <- function(dreme_out, database = Sys.getenv("TOMTOM_DB"), method = "bits"){
  # returns cowplot grids of de-novo matches | tomtom bestmatch motif.
  # can be input to cowplot::plot_grid(.$denovo, .$match, ncol = 2) for plotting
  info <- dreme_out$info

  denovo_motifPlots <- map(info$seq, ~{
    seq <- .
    motifPlot <- ggseqlogo::ggseqlogo(dreme_out$motifs[[seq]]@mat, method = method)
    return(motifPlot)
  })

  motif_db <- motifStack::importMatrix(database)

  match_motifPlots <- map(info$bestMatch, ~{
    match_id <- .

    if (match_id == "None") { return(NULL) } # if match is None then plot empty motif

    match_index <- grep(match_id %>% gsub("[-\\(\\)]", ".", .), names(motif_db)) # need to replace "-" and "()" for "." because importMatrix does.
    #TODO: Fix this for more sophisticated solution?
    match_index <- match_index[1] # return 1st match

    match_motif <- motif_db[[match_index]]@mat

    motifPlot <- ggseqlogo::ggseqlogo(match_motif, method = method)
    return(motifPlot)
  })

  denovo <- cowplot::plot_grid(plotlist = denovo_motifPlots, ncol = 1)

  match <- cowplot::plot_grid(plotlist = match_motifPlots, ncol = 1)

  plot_stacks <- list(denovo = denovo,
                      match = match)
  #cowplot::plot_grid(denovo, match, ncol = 2, labels = c("A", "B"))
  return(plot_stacks)
}

dreme_summary_plot <- function(dreme_output, ...){
  dreme_plots <- plot_denovo_matches(dreme_output, ...)

  info_table <- dreme_output$info %>%
    dplyr::select(seq, bestMatch, pvalue, pos_frac, neg_frac) %>%
    dplyr::mutate(ratio = pos_frac/neg_frac) %>%
    dplyr::mutate(pos_frac = signif(pos_frac * 100, 2),
                  neg_frac = signif(neg_frac * 100, 2),
                  ratio = signif(ratio, 2)) %>%
    gridExtra::tableGrob(.)

  cowplot::plot_grid(dreme_plots$denovo, dreme_plots$match, info_table, ncol = 3)
}


