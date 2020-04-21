dreme_input.character <- function(input){
  # character input should be file path unless input == "shuffle"
  if (input == "shuffle") return(input)

  dotargs::check_files_exist(input)
  return(input)
}

dreme_input.DNAStringSet <- function(input){
  # stringset is written to temporary fasta file
  write_fasta(input)
}

tomtom_input.character <- function(input){
  # TODO: Check is dreme xml format
    # parse to dreme results
    # read_xml %>% xml_name() == "dreme"
    # convert to df
    # importDremeXML()
  # check is
  dotargs::check_files_exist(input)
  return(input)
}

tomtom_input.data.frame <- function(input){
  # check data.frame is dreme_out object
  if (!is_dreme_results(input)) warn_dreme_results(input)

  path <- input$motifs %>%
    write_meme_list()

  out <- list(metadata = input,
              path = path)
  return(out)
}

tomtom_input.list <- function(input){
  # check list is universalmotif list
  if (!is_universalmotif_list(input)) error_universalmotif_list(list)

  df <- universalmotif_to_meme_df(input)

  path <- input %>%
    write_meme_list()

  out <- list(metadata = df,
              path = path)

  return(out)
}

tomtom_input.universalmotif <- function(input){

  df <- universalmotif_to_meme_df(input) %>%
    data.frame

  path <- input %>%
    write_meme_list()

  out <- list(metadata = df,
              path = path)

  return(out)
}

