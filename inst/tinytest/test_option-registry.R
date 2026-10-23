if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
  library("tinytest")
  library("highs")
}

# highs_available_solver_options() enumerates the bundled library's own option
# registry at runtime (Highs::getOptionName/getOptionType over ALL records;
# advanced == index >= getNumOptions(), since HighsOptions registers the
# user-settable options first -- see "// Advanced options" in HighsOptions.h)
# instead of returning the static HIGHS_OPTIONS table that
# inst/parse_highs_options.R regex-scraped from HighsOptions.h.
#
# Cross-validation: the static table is retained for one release; registry and
# table must agree as sets with matching types.  Once this has held through a
# release cycle, R/constants.R and inst/parse_highs_options.R can be deleted
# and this test reduces to the shape assertions.

opts <- highs_available_solver_options()

expect_true(is.data.frame(opts))
expect_identical(colnames(opts), c("option", "type", "advanced"))
expect_true(nrow(opts) >= 100L)                       # 156 in HiGHS 1.15.1
expect_true(all(opts$type %in% c("bool", "integer", "double", "string")))
expect_false(anyDuplicated(opts$option) > 0L)
expect_true(is.logical(opts$advanced) || all(opts$advanced %in% c(0L, 1L)))

# spot-check the advanced split
expect_false(as.logical(opts$advanced[opts$option == "presolve"]))
expect_false(as.logical(opts$advanced[opts$option == "output_flag"]))
expect_true(as.logical(opts$advanced[opts$option == "log_dev_level"]))

# registry vs static table (order-independent, sets + types).
# Known phantom: `pdlp_features_off` is a HighsOptionsStruct MEMBER that HiGHS
# never registers as an OptionRecord (HighsOptions.h declares and initializes
# it but there is no `new OptionRecord...` for it), so it was never actually
# settable via Highs::setOptionValue.  The regex-scraped static table
# advertised it anyway; the runtime registry correctly drops it.
known_phantoms <- "pdlp_features_off"
tab <- getFromNamespace("HIGHS_OPTIONS", "highs")
tab$option <- as.character(tab$option)
tab$type <- as.character(tab$type)
expect_equal(setdiff(tab$option, opts$option), known_phantoms,
  info = "only the documented phantom may be table-only")
expect_equal(setdiff(opts$option, tab$option), character(0),
  info = "the registry must not contain options unknown to the static table")
common <- intersect(tab$option, opts$option)
expect_equal(opts$type[match(common, opts$option)],
             tab$type[match(common, tab$option)],
  info = "registry and HIGHS_OPTIONS must agree on every common option's type")

# the registry must drive option access end-to-end, including advanced options
solver <- suppressWarnings(hi_new_solver(example_model()))
expect_equal(hi_solver_get_option(solver, "presolve"), "choose")
expect_true(is.logical(hi_solver_get_option(solver, "output_flag")))
expect_true(is.numeric(hi_solver_get_option(solver, "log_dev_level")))
