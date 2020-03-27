# gets path to meme:
# 1. checks local env
# 2. defaults to $HOME/meme/bin/
# 3. reports path used

# Checks valid meme suite install for a tool by:
# 1. gets version of meme-suite exe
# 2. if can't get version, assume uninstalled,
# 3. throw error by reporting path used for exe

# checks/gets/sets motif db for fimo/tomtom

# startup behavior
# 1. get path to meme
# 1.1 report path used
# 2. check install of each tool (dreme-py3/ame/fimo/tomtom)
# 2.1 Report success & failures

# TODO: use {crayon} for startup messages













# uses meta functions to compile {glue} calls to exes
# takes flags as input w/ default values, returns function that takes these as input
  ## uses functions to compile {glue} {system} calls

meme_glue_constructor <- function(meme_path, ...){
  rlang::dots
}


#' Transform NULL values to empty strings with {glue}
#'
#' Transformer function for glue to drop any NULL values.
#' Useful for concatenating outputs when NULLs are involved.
#' Not ment to be called directly.
#'
#' Fix NULL addition to glue by passing to glue::glue(.transformer)
#'
#'
#' @param text
#' @param envir
#'
#' @return
#'
#' @examples
#' # This will return an empty string or glue error depending on glue version
#' myString <- ""
#' glue::glue("hello {myString}")
#'
#' # Fix by passing transformer function
#' myString <- ""
#' glue::glue("hello {myString}", .transformer = null_transformer)
null_transformer <- function(text, envir) {
  # I got this from a github issue on the glue package or a stack overflow answer but forgot to grab the hyperlink...
  out <- glue::identity_transformer(text, envir)
  if (is.null(out)) {
    return("")
  }

  out
}
