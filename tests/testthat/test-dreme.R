skip_if(T, "Only works with dreme install")

fa <- "inst/extdata/fasta_ex/fa1.fa"
fa2 <- "inst/extdata/fasta_ex/fa2.fa"

dreme_out <- runDreme(fa, "shuffle", e = 39)

#####
stop("dev below")
runDreme(fa, fa2)
