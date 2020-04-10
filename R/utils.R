
outdir_name <- function(input, control){

  paste0(dirname(input), "/",
         basename(tools::file_path_sans_ext(input)),
         "_vs_",
         basename(tools::file_path_sans_ext(control)))
}
