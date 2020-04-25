skip_if(T, "Only works with dreme install")

fa <- duplicate_file("inst/extdata/fasta_ex/fa1.fa")
fa2 <- duplicate_file("inst/extdata/fasta_ex/fa2.fa")

# won't detect @ normal threshold
expect_null(runDreme(fa, "shuffle"))
expect_true(is_dreme_results(runDreme(fa, "shuffle", e = 39)))

dreme_out <- runDreme(fa, "shuffle", e = 39)

# test dreme catch error
empty_fa <- tempfile()
file.create(empty_fa)
expect_error(suppressMessages(runDreme(empty_fa, "shuffle"), "No sequences"))

# test universalmotif validators
expect_true(is_universalmotif_list(list(universalmotif::create_motif())))
expect_false(is_universalmotif_list(list(universalmotif::create_motif(), "no")))

# Test dreme_results validators
spec <- new_dreme_results()
bad <- spec[,1:5]

expect_true(is_dreme_results(dreme_out))
expect_false(is_dreme_results(bad))

expect_null(error_dreme_results(dreme_out))
expect_error(error_dreme_results(spec), "motifs column is empty")
expect_error(error_dreme_results(bad), "Missing columns")

#####
# test multiple dispatch
peaks <- "inst/extdata/peaks/peaks.tsv" %>%
  readr::read_tsv() %>%
  GRanges

seq <- get_sequence(peaks, dm.genome)
fa <- write_fasta(seq)


## test fasta path input
runDreme(fa, "shuffle", e = 49)
runDreme(seq, "shuffle", e = 49)

runDreme(fa, seq)
runDreme(fa, fa)

## test error checking
fa <- dremeR:::duplicate_file("inst/extdata/fasta_ex/fa1.fa")
expect_error(suppressMessages(runDreme(fa, "shuffle", et = 39)),
             "\"e\" instead of", class = "error")


###
#test tomtom

# dev tt_result_parse
# want data.frame w/ results & motif_match col w/ universalmotif object

#tomtom_xml_path <- tt_files$xml
#####
stop("dev below")
runDreme(fa, fa2)
