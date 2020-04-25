a <- handle_meme_path(util = "ame")

a_out <- list(status = 1,
              stdout = "output\n",
              stderr = "Some error")


# No user flags should not check errors
expect_error(
  suppressMessages(
    process_check_error(a_out, help_fun = ~{processx::run(a, "-h", error_on_status = F)$stderr},
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
                      flags_fun = ~{gsub("-", "_", .x)}),
  "\"evalue_report_threshold\" instead of: \"evalue_reorpt_threshold\"",
  class = "error"
)

# Also ensure usethis error is returned above
expect_error(
  process_check_error(a_out,
                      help_fun = ~{processx::run(a, "-h", error_on_status = F)$stderr},
                      user_flags = c("evalue_reorpt_threshold"),
                      flags_fun = ~{gsub("-", "_", .x)}),
  "\"evalue_report_threshold\" instead of: \"evalue_reorpt_threshold\"",
  class = "usethis_error"
)
