#' Handle error or printing message on successful DREME run
#'
#' @param processx_out processx ouput
#' @param silent whether or not to suppress printing all of stdout
#'
#' @return
#'
#' @examples
#' @noRd
print_process_stdout <- function(processx_out, silent = TRUE){
  process_check_error(processx_out)

  # leaving here incase needed for debugging
  #nmotifs_line <- dreme_nmotifs_line(processx_out$stdout)
  #if (!silent) message(nmotifs_line)

  if (!silent) message(processx_out$stdout)

}

#' Print stderr from processx output
#'
#' @param processx_out
#' @param silent whether to suppress stderr printing or not
#'
#' @return
#'
#' @noRd
print_process_stderr <- function(processx_out, silent = TRUE){
  process_check_error(processx_out)

  if (!silent) message(processx_out$stderr)

}

#' Check if processx process completed successfully
#'
#' @param processx_out output of processx::run(error_on_status = F)
#' @param help_fun function to produce --help lines (processx call, for example). Can be rlang formula like ~{.x}.
#'
#' @return NULL if exit status 0, otherwise print all stdout + stderr
#'
#' @examples
#'
#' @noRd
process_check_error <- function(processx_out, help_fun = NULL, user_flags = NULL, flags_fun = NULL){
  if (processx_out$status != 0 & (is.null(help_fun) | is.null(user_flags))){
    cat(processx_out$stdout)
    stop(processx_out$stderr)
  } else if (processx_out$status != 0) {
    stopifnot(is.function(help_fun) | class(help_fun) == "formula")
    stopifnot(is.function(flags_fun) | is.null(flags_fun) | class(flags_fun) == "formula")

    # Allows lambda
    help_fun <- rlang::as_function(help_fun)

    cat(processx_out$stdout)
    message(processx_out$stderr)
    #usethis::ui_warn(c("\n", processx_out$stderr))

    help_fun() %>%
      dotargs::get_help_flag_names(processx = T) %>%
      dotargs::suggest_flag_names(user_flags, .fun = flags_fun) %>%
      #suggest_flag_names(user_flags, .fun = flags_fun) -> o
      #return(o)
      dotargs::error_suggest_flag_names()

    usethis::ui_stop("Shell process had non-zero exit status.")

  }
  return(NULL)
}
