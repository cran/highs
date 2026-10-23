#!/usr/bin/env Rscript
# audit_symbols.R -- post-build gate replicating the compiled-code check that
# R CMD check runs on the installed .so (a symbol-table scan; see Writing R
# Extensions "Checking packages" and tools:::so_symbol_names_table), so a
# violation in the bundled HiGHS static library fails the build BEFORE a CRAN
# submission.  OPT-IN: inst/build_highs.sh invokes this after `make install`
# only when R_HIGHS_SYMBOL_CHECK is set (strict|warn) -- ordinary user/CRAN
# installs never run it (CRAN's own scan is the gate there).  The Taskfile
# check targets set R_HIGHS_SYMBOL_CHECK=strict.  Scans lib*/libhighs.a under
# the given directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: audit_symbols.R <highslib-dir-or-archive>")

mode <- Sys.getenv("R_HIGHS_SYMBOL_CHECK", "strict")
if (identical(mode, "skip")) {
  message("hi-symbol-gate: skipped (R_HIGHS_SYMBOL_CHECK=skip)")
  quit(status = 0L)
}

nm <- Sys.which("nm")
if (!nzchar(nm)) {
  message("hi-symbol-gate: 'nm' not found; gate skipped (CRAN builders have it)")
  quit(status = 0L)
}

target <- args[[1L]]
archives <- if (dir.exists(target)) {
  Sys.glob(file.path(target, c("lib", "lib64"), "libhighs.a"))
} else {
  target
}
archives <- archives[file.exists(archives)]
if (!length(archives)) stop("hi-symbol-gate: no libhighs.a under ", target)

# Names flagged by CRAN's scan (tools:::so_symbol_names_table) plus platform
# aliases.  Matched against nm names after stripping ONE leading underscore
# (Mach-O).  fprintf/vfprintf/fopen/fclose/fwrite to variable FILE* are
# deliberately NOT banned -- CRAN does not flag them.
banned_exact <- c(
  "printf", "puts", "putchar", "vprintf", "sprintf", "vsprintf",
  "abort", "exit", "_exit", "_Exit",
  "stdout", "stderr",           # ELF data symbols
  "__stdoutp", "__stderrp",     # macOS (___stdoutp after one _ stripped)
  "__acrt_iob_func",            # Windows UCRT stdout/stderr accessor
  "__assert_rtn", "__assert_fail", "assert",
  "rand", "random", "rand_r", "srand", "srandom", "srand48",
  "__sprintf_chk", "__vsprintf_chk"
)
banned_regex <- "ZSt4cout|ZSt4cerr|ZSt4clog|ZNSt3__14coutE|ZNSt3__14cerrE|ZNSt3__14clogE"

report <- character()
for (a in archives) {
  out <- suppressWarnings(system2(nm, c("-Pg", shQuote(a)), stdout = TRUE,
                                  stderr = FALSE))
  # POSIX -P format: "name type [value [size]]"; undefined entries have type U.
  fields <- strsplit(out, "[[:space:]]+")
  keep <- vapply(fields, function(f) length(f) >= 2L && f[[2L]] == "U",
                 logical(1L))
  und <- vapply(fields[keep], `[[`, character(1L), 1L)
  stripped <- sub("^_", "", und)
  bad <- und[stripped %in% banned_exact | grepl(banned_regex, und)]
  if (length(bad)) {
    report <- c(report, paste0(a, ":"),
                paste0("  ", sort(unique(bad))))
  }
}

if (length(report)) {
  writeLines(c(
    "hi-symbol-gate: CRAN-banned symbols found (these WILL be flagged by",
    "R CMD check's compiled-code scan; see tools:::so_symbol_names_table):",
    report))
  if (identical(mode, "warn")) {
    message("hi-symbol-gate: R_HIGHS_SYMBOL_CHECK=warn -- continuing anyway")
    quit(status = 0L)
  }
  quit(status = 1L)
}

message("hi-symbol-gate: clean (", paste(basename(archives), collapse = ", "),
        ")")
quit(status = 0L)
