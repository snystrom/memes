
runDreme <- function(input, control, evalue = 0.05, n_motifs = 10, sec = 18000, alph = c("dna", "rna", "protein"),
                     ngen = 100,
                     outdir = paste0(dirname(input), "/", basename(tools::file_path_sans_ext(input)),  "_vs_", basename(tools::file_path_sans_ext(control))),
                     meme_path = Sys.getenv("MEME_DIR"),
                     mink = 3,
                     maxk = 8){}


outdir_name <- function(input, control){

  paste0(dirname(input), "/",
         basename(tools::file_path_sans_ext(input)),
         "_vs_",
         basename(tools::file_path_sans_ext(control)))
}


runDreme <- function(input, control, outdir = outdir_name(input, control), ...){
  argDict <- c(n_motifs = "m",
               sec = "t",
               evalue = "e",
               seed = "s",
               input = "p",
               control = "n",
               outdir = "oc")
  print(outdir)
  # dealing with handling outdir & other args which don't get evaluated during match.call()
  args <- as.list(match.call()) %>%
    {.[-1]} %>%
    {if ("outdir" %in% names(.)) within(., rm("outdir")) else return(.)} %>%
    dotargs::dotsToArgs(argDict) %>%
    {.[!grepl("-n shuffle", .)]}

  return(args)
  #system2("dreme", c(args))

}

runDreme("input.fa", "shuffle", outdir = outdir_name("input.fa", "shuffle"))
runDreme("input.fa", "shuffle", sec = 1000)

system2(path.expand("~/meme/bin/dreme"), "-h")
processx::run(path.expand("~/meme/bin/dreme"), c(paste0("-p ", here::here("inst/extdata/fasta_ex/fa1.fa")), "-oc /inst/psx"))
