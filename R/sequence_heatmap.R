# colors from:
# https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
DNA_COLORS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62839")                                                                                                                                                             
RNA_COLORS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", U = "#D62839")                                                                                                                                                             
AA_COLORS <- c(G = "#058644", S = "#058644", T = "#058644", Y = "#058644", 
                C = "#058644", Q = "#720091", N = "#720091", K = "#0046C5", 
                R = "#0046C5", H = "#0046C5", D = "#C5003E", E = "#C5003E", 
                A = "#2E2E2E", V = "#2E2E2E", L = "#2E2E2E", I = "#2E2E2E", 
                P = "#2E2E2E", W = "#2E2E2E", F = "#2E2E2E", M = "#2E2E2E")

#' Convert character vector of sequences to melted data.frame
#' 
#' @param x a character vector of sequences
#' @return a data.frame with "id", "position" and "letters" cols. where id is
#'   the sequence name, position is position within the sequence, and letter is
#'   the letter at the given position. Sequences are numbered by their order in
#'   x.
#' 
#' @noRd
sequence_to_df <- function(x){
  
  purrr::map_int(x, nchar) %>% 
    {length(unique(.)) == 1} %>% 
    if (!.){
      stop("Sequences must be equal length", call. = FALSE)
    }
  
  purrr::imap_dfr(x, ~{
    seq <- .x
    id <- .y
    letters <- strsplit(seq, "")[[1]]
    data.frame("id" = rep(id, nchar(seq)),
               "position" = seq_along(letters),
               "letters" = letters
               )
  })
}

#' Make ggplot heatmap from data.frame of sequences
#'
#' @param x result of sequence_to_df
#'
#' @param alph alphabet to use (DNA, RNA, AA), determines color scheme
#'
#' @return a ggplot object with the heatmap of sequences
#' 
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual scale_y_continuous
#'   theme_minimal
#' 
#' @noRd
sequence_df_to_heatmap <- function(x, alph = c("DNA", "RNA", "AA")){
  alph <- toupper(alph)
  alph <- match.arg(alph, c("DNA", "RNA", "AA"))
  
  colors <- switch(alph,
                   DNA = DNA_COLORS,
                   RNA = RNA_COLORS,
                   AA = RNA_COLORS)
  
  x %>% 
    dplyr::mutate(position = factor(position)) %>% 
    ggplot(aes(position, id)) +
      geom_tile(aes(fill = letters)) +
      scale_fill_manual(values = colors) +
      scale_y_continuous(expand = c(0,0)) +
      theme_minimal()
}

#' Take character vector of sequences and make ggplot2 heatmap
#'
#' @param sequences character vector of sequences 
#' @param alph alphabet to use
#' @return a ggplot2 heatmap ranked in order of sequences
#' @noRd
sequence_to_heatmap <- function(sequences, alph = c("DNA", "RNA", "AA")){
  
  sequences %>% 
    sequence_to_df %>% 
    sequence_df_to_heatmap(alph = alph)
  
}

#' Visualize a heatmap of aligned sequences
#' 
#' Sometimes it is useful to visualize individual motif matches in aggregate to
#' understand how sequence variability contributes to motif matches. This
#' function creates a heatmap where each row represents a single sequence and
#' each column represents a position. Sequences are optionally aggregated into a
#' sequence logo aligned in register with the heatmap to visualize how sequence
#' variability contributes to motif makeup.
#' 
#' @param sequence character vector of sequences, plot will be ranked in order
#'   of the sequences. Each sequence must be equal length.
#' @param title title of the plot
#' @param logo whether to include a sequence logo above the heatmap
#' @param alph alphabet colorscheme to use. One of: DNA, RNA, AA.
#' @param title_hjust value from 0 to 1 determining the horizontal justification
#'   of the title. Default: 0.
#' @param heights ratio of logo:heatmap heights. Given as: c(logo_height,
#' heatmap_height). Values are not absolute. Ignored when logo = FALSE.
#' 
#' @return a ggplot object of the sequence heatmap ranked by the order of
#'   sequences
#'   
#' @seealso runFimo
#' @export
#' 
#' @importFrom ggplot2 labs theme element_blank element_text
#' 
#' @examples
#' data(example_fimo, package = "memes")
#' genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3
#' motifs <- add_sequence(example_fimo, genome)
#' plot_sequence_heatmap(motifs$sequence)
plot_sequence_heatmap <- function(sequence, title = NULL, logo = TRUE, 
                                  alph = c("DNA", "RNA", "AA"), 
                                  title_hjust = 0, heights = c(1,5)){
  
  if (length(heights) != 2) {
    stop("heights must be a numeric vector of lenght 2", call. = FALSE)
  }
  
  pwm <- universalmotif::create_motif(sequence) %>%
    universalmotif::view_motifs() +
      labs(title = title) +
      theme(axis.text.x = element_blank(), 
            plot.title = element_text(hjust = title_hjust))
  
  heatmap <- sequence_to_heatmap(sequence, alph) + 
    labs(y = NULL) +
    theme(axis.text.y = element_blank(), 
          legend.title = element_blank(),
          legend.position = "none"
          ) 
    
  if (logo){
    return(patchwork::wrap_plots(pwm, heatmap, ncol = 1, heights = heights))
  } else {
    heatmap +
      labs(title = title) +
      theme(plot.title = element_text(hjust = title_hjust))
  }
}
