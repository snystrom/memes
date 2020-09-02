a <- search_meme_path(util = "ame")

a_out <- list(status = 1,
              stdout = "output\n",
              stderr = "Some error")


test_that("process_check_error suggests correctly", {
  # No user flags should not check errors
  expect_error(
    suppressMessages(
      process_check_error(a_out, help_fun = ~{ame_help(a)},
      flags_fun = ~{gsub("-", "_", .x)})
      ),
    "Some error",
    class = "error"
  )

  # lambda functions should work as input
  # user flags w/typo should get suggested
  expect_error(
    process_check_error(a_out,
                        help_fun = ~{processx::run(a, "-h", error_on_status = F)$stderr},
                        user_flags = c("evalue_reorpt_threshold"),
                        flags_fun = ~{gsub("-", "_", .x)}, default_help_fun = TRUE),
    "\"evalue_report_threshold\" instead of: \"evalue_reorpt_threshold\"",
    class = "error"
  )

  # Also ensure usethis error is returned above
  expect_error(
    process_check_error(a_out,
                        help_fun = ~{processx::run(a, "-h", error_on_status = F)$stderr},
                        user_flags = c("evalue_reorpt_threshold"),
                        flags_fun = ~{gsub("-", "_", .x)}, default_help_fun = TRUE),
    "\"evalue_report_threshold\" instead of: \"evalue_reorpt_threshold\"",
    class = "usethis_error"
  )

})
