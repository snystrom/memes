
#' Helper for writing unique directory names
#'
#' @param input input file path
#' @param control control file path
#'
#' @return directory name in the style of: "input_vs_output".
#' @export
#'
#' @examples
#' outdir_name("condition1.fa", "backgroundSequence.fa")
outdir_name <- function(input, control){

  paste0(dirname(input), "/",
         basename(tools::file_path_sans_ext(input)),
         "_vs_",
         basename(tools::file_path_sans_ext(control)))
}
