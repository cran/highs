if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}

# The C API functions registered with R_RegisterCCallable (src/highs_ccallable.cpp
# through the generated src/highs_ccallable_list.h) and the wrappers in the
# generated inst/include/highs_ccallable.h must be the same set.

header <- system.file("include", "highs_ccallable.h", package = "highs")
expect_true(file.exists(header))

registered <- highs:::.hi_ccallable_names()
expect_true(is.character(registered))
expect_true(length(registered) >= 170L)   # 177 at HiGHS 1.15.1
expect_true(all(startsWith(registered, "Highs_")))
expect_equal(anyDuplicated(registered), 0L)

src <- readLines(header, warn = FALSE)
wrapped <- regmatches(src, regexpr("^static inline [^(]*\\b(Highs_[A-Za-z0-9_]+)\\(", src))
wrapped <- sub("^.*\\b(Highs_[A-Za-z0-9_]+)\\($", "\\1", wrapped)
resolved <- regmatches(src, regexpr("HIGHS_CCALLABLE_RESOLVE\\(Highs_[A-Za-z0-9_]+", src))
resolved <- sub("HIGHS_CCALLABLE_RESOLVE\\(", "", resolved)
expect_equal(wrapped, registered)   # same names, same (header) order
expect_equal(wrapped, resolved)     # each wrapper resolves its own name

# The callback types and the constants the C API needs are in the header.
expect_true(any(grepl("^typedef struct HighsCallbackDataOut \\{", src)))
expect_true(any(grepl("^typedef struct HighsCallbackDataIn \\{", src)))
expect_true(any(grepl("^typedef void \\(\\*HighsCCallbackType\\)", src)))
expect_true(any(grepl("^static const HighsInt kHighsModelStatusOptimal = 7;", src)))
expect_true(any(grepl("^static const HighsInt kHighsCallbackLogging = 0;", src)))

# Every registered name resolves through R_GetCCallable, as a consumer would.
ok <- highs:::.hi_ccallable_check()
expect_equal(names(ok), registered)
expect_true(all(ok))

# Drift check against the bundled HiGHS: regenerate from the headers of the
# submodule and compare with the committed files and the installed registry.
# Only possible in the source checkout (the generator and the bundled headers
# are not installed); R_HIGHS_SOURCE_DIR overrides the `inst/tinytest` guess.
src_root <- Sys.getenv("R_HIGHS_SOURCE_DIR", NA_character_)
if (is.na(src_root)) {
    guess <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(guess, "tools", "gen_ccallable.R")) &&
        file.exists(file.path(guess, "inst", "HiGHS", "highs", "interfaces", "highs_c_api.h")))
        src_root <- guess
}
if (!is.na(src_root)) {
    gen_env <- new.env()
    source(file.path(src_root, "tools", "gen_ccallable.R"), local = gen_env)
    gen <- gen_env$hi_cc_generate(src_root)
    expect_equal(gen$names, registered)
    expect_identical(readLines(file.path(src_root, "src", "highs_ccallable_list.h"), warn = FALSE),
                     gen$list_h)
    expect_identical(readLines(file.path(src_root, "inst", "include", "highs_ccallable.h"), warn = FALSE),
                     gen$header)
}

# Opt-in end-to-end test (R_HIGHS_LINK_TEST=1): compile a C consumer that
# includes only highs_ccallable.h, load it, solve a 2-variable LP through the
# registered functions with a logging callback attached.
if (identical(Sys.getenv("R_HIGHS_LINK_TEST"), "1") && nzchar(Sys.which("R"))) {
    src_c <- system.file("tinytest", "ccallable_consumer.c", package = "highs")
    expect_true(file.exists(src_c))
    wd <- tempfile("ccallable_"); dir.create(wd)
    file.copy(src_c, wd)
    owd <- setwd(wd); on.exit(setwd(owd), add = TRUE)
    rc <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "SHLIB", "ccallable_consumer.c"),
                  stdout = "shlib.log", stderr = "shlib.log",
                  env = sprintf("PKG_CPPFLAGS=-I%s", shQuote(system.file("include", package = "highs"))),
                  wait = TRUE)
    if (rc != 0L) writeLines(readLines("shlib.log"))
    expect_equal(rc, 0L)
    so <- file.path(wd, paste0("ccallable_consumer", .Platform$dynlib.ext))
    expect_true(file.exists(so))
    if (rc == 0L && file.exists(so)) {
        dll <- dyn.load(so)
        res <- .Call("ccallable_consumer", PACKAGE = dll[["name"]])
        dyn.unload(so)
        expect_equal(res$version, hi_solver_version(example_solver("LP")))
        expect_equal(res$sizeof_int, 4L)
        expect_equal(res$status_pass, 0L)
        expect_equal(res$status_run, 0L)
        expect_equal(res$model_status, 7L)      # kHighsModelStatusOptimal
        expect_equal(res$objective, 1, tolerance = 1e-9)
        expect_equal(sum(res$col_value), 1, tolerance = 1e-9)
        expect_true(res$callback_calls > 0L)    # the logging callback fired
        expect_equal(res$callback_item_ok, 1L)  # struct layout agrees with the library
    }
}
