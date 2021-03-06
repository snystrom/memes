% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_xml.R
\name{importStremeXML}
\alias{importStremeXML}
\title{Import Streme output from previous run}
\usage{
importStremeXML(streme_xml_path)
}
\arguments{
\item{streme_xml_path}{path to streme.xml file}
}
\value{
data.frame with statistics for each discovered motif. The `motifs`
  column contains a universalmotif object representation in PCM format of
  each DREME motif. If no motifs are discovered, returns NULL.
}
\description{
Import Streme output from previous run
}
\examples{
streme_xml <- system.file("extdata", "streme.xml", package = "memes")
importStremeXML(streme_xml)
}
\seealso{
[runStreme()]
}
