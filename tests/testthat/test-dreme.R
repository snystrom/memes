skip_if(T, "Only works with dreme install")

fa <- "inst/extdata/fasta_ex/fa1.fa"
fa2 <- "inst/extdata/fasta_ex/fa2.fa"

dreme_out <- runDreme(fa, "shuffle", e = 39)

###
#test tomtom

# dev tt_result_parse
# want data.frame w/ results & motif_match col w/ universalmotif object

#tomtom_xml_path <- tt_files$xml
#####
stop("dev below")
runDreme(fa, fa2)
