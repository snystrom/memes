% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meme.R
\name{runMeme}
\alias{runMeme}
\alias{runMeme.list}
\alias{runMeme.BStringSetList}
\alias{runMeme.default}
\title{Identify motifs with MEME}
\usage{
runMeme(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{list}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{BStringSetList}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{default}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)
}
\arguments{
\item{input}{path to fasta, Biostrings::BStringSet list, or list of
Biostrings::BStringSet (can generate using \code{get_sequence()})}

\item{control}{any data type as in \code{input}, or a character vector of
\code{names(input)} to use those regions as control sequences. Using sequences
as background requires an alternative objective function. Users must pass a non-default value of
\code{objfun} to \code{...} if using a non-NA control set (default: NA)}

\item{outdir}{(default: "auto") Directory where output data will be stored.}

\item{alph}{one of c("dna", "rna", "protein") or path to alphabet file (default: "dna").}

\item{parse_genomic_coord}{\code{logical(1)} whether to parse genomic coordinates
from fasta headers. Requires headers are in the form: "chr:start-end", or
will result in an error. Automatically set to \code{FALSE} if \code{alph = "protein"}. This setting only needs to be changed if using a custom-built
fasta file without genomic coordinates in the header.}

\item{combined_sites}{\code{logical(1)} whether to return combined sites
information (coerces output to list) (default: FALSE)}

\item{silent}{Whether to suppress printing stdout to terminal (default: TRUE)}

\item{meme_path}{path to "meme/bin/". If unset, will use default search
behavior:
\enumerate{
\item \code{meme_path} setting in \code{options()}
\item \code{MEME_PATH} setting in \code{.Renviron} or \code{.bashrc}
}}

\item{...}{additional arguments passed to MEME (see below)}
}
\value{
MEME results in universalmotif_df format (see:
\code{\link[universalmotif:tidy-motifs]{universalmotif::to_df()}}). \code{sites_hits} is a nested data.frame
column containing the position within each input sequence of matches to the
identified motif.
}
\description{
MEME performs \emph{de-novo} discovery of ungapped motifs present in the input
sequences. It can be used in both discriminative and non-discriminative
modes.
}
\details{
Note that MEME can take a long time to run. The more input sequences used,
the wider the motifs searched for, and the more motifs MEME is asked to
discover will drastically affect runtime. For this reason, MEME usually
performs best on a few (<50) short (100-200 bp) sequences, although this is
not a requirement. Additional details on how data size affects runtime can be
found \href{https://groups.google.com/g/meme-suite/c/7b7PBr6RzJk}{here}.

MEME works best when specifically tuned to the analysis question. The default
settings are unlikely to be ideal. It has several complex arguments
\href{http://meme-suite.org/doc/meme.html}{documented here}, which \code{runMeme()}
accepts as R function arguments (see details below).

If discovering motifs within ChIP-seq, ATAC-seq, or similar peaks, MEME may perform
best if using sequences flaking the summit (the site of maximum signal) of
each peak rather than the center. ChIP-seq or similar data can also benefit
from setting \verb{revcomp = TRUE, minw = 5, maxw = 20}. For more tips on using
MEME to analyze ChIP-seq data, see the following
\href{https://groups.google.com/forum/#\%21topic/meme-suite/rIbjIHbcpAE}{tips page}.
\subsection{Additional arguments}{

\code{\link[=runMeme]{runMeme()}} accepts all valid arguments to meme as arguments passed to \code{...}.
For flags without values, pass them as \code{flag = TRUE}. The \code{dna}, \code{rna}, and
\code{protein} flags should instead be passed to the \code{alph} argument of
\code{\link[=runMeme]{runMeme()}}.  The arguments passed to MEME often have many interactions
with each other, for a detailed description of each argument see
\href{meme-suite.org/doc/meme.html}{MEME Commandline Documentation}.
}
}
\section{Citation}{
If you use \code{runMeme()} in your analysis, please cite:

Timothy L. Bailey and Charles Elkan, "Fitting a mixture model by expectation
maximization to discover motifs in biopolymers", Proceedings of the Second
International Conference on Intelligent Systems for Molecular Biology, pp.
28-36, AAAI Press, Menlo Park, California, 1994.
\href{https://tlbailey.bitbucket.io/papers/ismb94.pdf}{pdf}
}

\section{Licensing}{
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}

\examples{
if (meme_is_installed()) {
seqs <- universalmotif::create_sequences("CCRAAAW", seqnum = 4)
names(seqs) <- 1:length(seqs)
runMeme(seqs, parse_genomic_coord = FALSE)

}

}
