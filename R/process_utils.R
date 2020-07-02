#' Handle error or printing message on successful DREME run
#'
#' @param processx_out processx ouput
#' @param silent whether or not to suppress printing all of stdout
#'
#' @return
#'
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
#' @param processx_out output from processx::run(error_on_status = FALSE)
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
#' @param user_flags flags from cmd_args_to_flags %>% crystalizeFlags
#' @param flags_fun passed to cmdr::cmd_help_flags_similar(.fun)
#' @param default_help_fun `logical(1)` whether to use default
#'   cmd_help_parse_flags behavior, or simply evaluate help_fun fully to suggest
#'   names. Set to TRUE when help_fun returns vector of flag names to suggest against.
#'
#' @return NULL if exit status 0, otherwise print all stdout + stderr
#'
#' @examples
#'
#' @noRd
process_check_error <- function(processx_out, help_fun = NULL, user_flags = NULL, flags_fun = NULL, default_help_fun = FALSE){
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

    if (default_help_fun){
      help_fun() %>%
        cmdr::cmd_help_parse_flags(processx = T) %>%
        cmdr::cmd_help_flags_similar(user_flags, .fun = flags_fun) %>%
        cmdr::cmd_help_flags_suggest()

    } else {
      help_fun() %>%
        cmdr::cmd_help_flags_similar(user_flags, .fun = flags_fun) %>%
        cmdr::cmd_help_flags_suggest()
    }

    usethis::ui_stop("Shell process had non-zero exit status.")

  }
  return(NULL)
}
