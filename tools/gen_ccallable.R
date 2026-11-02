#!/usr/bin/env Rscript
# tools/gen_ccallable.R -- generate the R_RegisterCCallable registry of the
# HiGHS C API from the headers of the bundled HiGHS.
#
# Reads   inst/HiGHS/highs/interfaces/highs_c_api.h
#         inst/HiGHS/highs/lp_data/HighsCallbackStruct.h
# Writes  src/highs_ccallable_list.h      X-macro list of every C API function;
#                                         src/highs_ccallable.cpp registers it
#         inst/include/highs_ccallable.h  the consumer header: constants and
#                                         callback structs copied from the two
#                                         files above, plus one static inline
#                                         wrapper per function that resolves
#                                         R_GetCCallable("highs", "<name>")
#
# Usage, from the package root, with the submodule checked out:
#   Rscript tools/gen_ccallable.R           # (re)write the two files
#   Rscript tools/gen_ccallable.R --check   # exit status 1 if they are stale
#
# Run it at every HiGHS bump and commit the result.  The parser accepts
# exactly the shapes highs_c_api.h uses today (one prototype per `;`, no
# function-pointer parameters other than through a typedef, no variadics)
# and stops on anything else, so a new upstream construct shows up here
# instead of as a silent omission.  inst/tinytest/test_ccallable.R compares
# the committed files against a fresh generation when it runs in the source
# checkout.

hi_cc_api_header <- "inst/HiGHS/highs/interfaces/highs_c_api.h"
hi_cc_callback_header <- "inst/HiGHS/highs/lp_data/HighsCallbackStruct.h"
hi_cc_list_out <- "src/highs_ccallable_list.h"
hi_cc_header_out <- "inst/include/highs_ccallable.h"

hi_cc_squeeze <- function(x) trimws(gsub("[[:space:]]+", " ", x))

# Collapse a header to its `;`-terminated statements: no comments, no
# preprocessor lines, no extern "C" braces.
hi_cc_statements <- function(lines) {
  lines <- lines[!grepl("^[[:space:]]*#", lines)]
  lines <- sub("//.*$", "", lines)
  txt <- paste(lines, collapse = "\n")
  txt <- gsub("(?s)/\\*.*?\\*/", " ", txt, perl = TRUE)
  txt <- gsub("extern[[:space:]]*\"C\"[[:space:]]*\\{", "", txt)
  n_brace <- lengths(regmatches(txt, gregexpr("\\}", txt)))
  if (n_brace != 1L)
    stop("gen_ccallable: expected exactly one closing brace (end of extern \"C\") in ",
         "highs_c_api.h, found ", n_brace)
  txt <- sub("\\}", "", txt)
  stmts <- hi_cc_squeeze(strsplit(txt, ";", fixed = TRUE)[[1L]])
  stmts[nzchar(stmts)]
}

hi_cc_parse_api <- function(lines) {
  stmts <- hi_cc_statements(lines)
  is_fun <- grepl("^.*[^A-Za-z0-9_]Highs_[A-Za-z0-9_]+ ?\\(.*\\)$", stmts)
  is_const <- grepl("^(static )?const [A-Za-z0-9_* ]*\\b[A-Za-z_][A-Za-z0-9_]* ?= ", stmts)
  other <- stmts[!is_fun & !is_const]
  if (length(other))
    stop("gen_ccallable: statement not understood in highs_c_api.h:\n  ",
         paste(other, collapse = "\n  "))

  # Constants: verbatim apart from whitespace and an added `static` on the
  # few upstream lines that lack it (a non-static initialized const has
  # external linkage in C and would be defined once per translation unit).
  consts <- stmts[is_const]
  consts <- sub("^const ", "static const ", consts)
  consts <- paste0(consts, ";")

  # A prototype declared `static` in the header has internal linkage: it is
  # not reachable from highs_ccallable.cpp and cannot be registered (upstream
  # has one, the deprecated Highs_compilationDate).  Skipped and reported.
  fun_stmts <- stmts[is_fun]
  is_static <- grepl("^static ", fun_stmts)
  skipped <- sub("^.*[^A-Za-z0-9_](Highs_[A-Za-z0-9_]+) ?\\(.*$", "\\1", fun_stmts[is_static])
  funs <- lapply(fun_stmts[!is_static], function(s) {
    m <- regmatches(s, regexec("^(.*?) ?\\b(Highs_[A-Za-z0-9_]+) ?\\((.*)\\)$", s, perl = TRUE))[[1L]]
    ret <- hi_cc_squeeze(m[2L])
    name <- m[3L]
    params <- hi_cc_squeeze(m[4L])
    if (grepl("[()]", params) || grepl("\\.\\.\\.", params))
      stop("gen_ccallable: parameter list of ", name, " not understood: ", params)
    if (identical(params, "void") || !nzchar(params)) {
      params <- character()
      args <- character()
    } else {
      params <- hi_cc_squeeze(strsplit(params, ",", fixed = TRUE)[[1L]])
      args <- sub("^.*[^A-Za-z0-9_]([A-Za-z0-9_]+)$", "\\1", params)
      if (any(args == params))
        stop("gen_ccallable: unnamed parameter in ", name, ": ", params[args == params][1L])
      if (anyDuplicated(args) || any(args %in% c("fn", "fn_t", "highs_ccallable_lookup")))
        stop("gen_ccallable: parameter names of ", name, " clash with the wrapper: ",
             paste(args, collapse = ", "))
    }
    if (!nzchar(ret)) stop("gen_ccallable: no return type for ", name)
    list(ret = ret, name = name, params = params, args = args)
  })
  names(funs) <- vapply(funs, `[[`, "", "name")
  if (anyDuplicated(names(funs)))
    stop("gen_ccallable: duplicate prototype: ", names(funs)[duplicated(names(funs))][1L])
  list(consts = consts, funs = funs, skipped = skipped)
}

# The typedefs of HighsCallbackStruct.h, verbatim: everything between the
# extern "C" line and the #ifdef __cplusplus that closes it.
hi_cc_callback_block <- function(lines) {
  open <- grep("^extern \"C\" \\{", lines)
  close <- grep("^#ifdef __cplusplus", lines)
  if (length(open) != 1L || length(close) != 2L || close[2L] <= open)
    stop("gen_ccallable: HighsCallbackStruct.h does not have the expected layout")
  block <- lines[(open + 1L):(close[2L] - 1L)]
  if (!grepl("^#endif", block[1L]))   # closes the `#ifdef __cplusplus` of extern "C"
    stop("gen_ccallable: HighsCallbackStruct.h does not have the expected layout")
  block <- block[-1L]
  while (length(block) && !nzchar(trimws(block[1L]))) block <- block[-1L]
  while (length(block) && !nzchar(trimws(block[length(block)]))) block <- block[-length(block)]
  if (!any(grepl("typedef struct HighsCallbackDataOut", block)) ||
      !any(grepl("HighsCCallbackType", block)))
    stop("gen_ccallable: HighsCallbackStruct.h block lacks the expected typedefs")
  block
}

# prefix item, item, ... itemsuffix  filled greedily to `width` columns,
# continuation lines indented by `indent`; if the first item does not fit
# after the prefix, break right after the prefix (clang-format style).
hi_cc_fill <- function(prefix, items, suffix, indent, width = 80L) {
  if (!length(items)) return(paste0(prefix, suffix))
  pieces <- paste0(items, c(rep(",", length(items) - 1L), suffix))
  out <- character()
  if (nchar(prefix) + nchar(pieces[1L]) > width) {
    out <- prefix
    line <- paste0(indent, pieces[1L])
  } else {
    line <- paste0(prefix, pieces[1L])
  }
  for (p in pieces[-1L]) {
    if (nchar(line) + 1L + nchar(p) > width) {
      out <- c(out, line)
      line <- paste0(indent, p)
    } else {
      line <- paste0(line, " ", p)
    }
  }
  c(out, line)
}

hi_cc_wrapper <- function(f) {
  c(hi_cc_fill(paste0("static inline ", f$ret, " ", f$name, "("), f$params, ") {", "    "),
    hi_cc_fill(paste0("  typedef ", f$ret, " (*fn_t)("), f$params, ");", "      "),
    paste0("  HIGHS_CCALLABLE_RESOLVE(", f$name, ", fn_t);"),
    hi_cc_fill(if (identical(f$ret, "void")) "  fn(" else "  return fn(", f$args, ");", "      "),
    "}")
}

hi_cc_const_lines <- function(consts) {
  unlist(lapply(consts, function(s) {
    if (nchar(s) <= 80L) return(s)
    # long string constants: break after `=`, as upstream does
    strsplit(sub(" = ", " =\n    ", s, fixed = TRUE), "\n", fixed = TRUE)[[1L]]
  }))
}

hi_cc_generate <- function(root = ".") {
  api_lines <- readLines(file.path(root, hi_cc_api_header), warn = FALSE)
  cb_lines <- readLines(file.path(root, hi_cc_callback_header), warn = FALSE)
  api <- hi_cc_parse_api(api_lines)
  cb_block <- hi_cc_callback_block(cb_lines)
  nms <- names(api$funs)
  n <- length(nms)

  # -- src/highs_ccallable_list.h -----------------------------------------
  pad <- max(nchar(nms)) + 4L
  entries <- sprintf("  X(%s)%s\\", nms, strrep(" ", pad - nchar(nms)))
  entries[n] <- sub("[[:space:]]*\\\\$", "", entries[n])
  list_h <- c(
    "// src/highs_ccallable_list.h -- GENERATED by tools/gen_ccallable.R from",
    "// inst/HiGHS/highs/interfaces/highs_c_api.h; do not edit, regenerate.",
    sprintf("// Every function of the HiGHS C API, in header order (%d).", n),
    "// src/highs_ccallable.cpp registers each one with R_RegisterCCallable.",
    if (length(api$skipped))
      sprintf("// Not registrable, declared static in the header: %s.",
              paste(api$skipped, collapse = ", ")),
    "#ifndef HIGHS_CCALLABLE_LIST_H",
    "#define HIGHS_CCALLABLE_LIST_H",
    "",
    "#define HIGHS_CCALLABLE_LIST(X) \\",
    entries,
    "",
    "#endif /* HIGHS_CCALLABLE_LIST_H */")

  # -- inst/include/highs_ccallable.h -------------------------------------
  header <- c(
    "/* highs_ccallable.h -- call the HiGHS C API bundled in the R package 'highs'",
    " * from another R package's compiled code, without linking to HiGHS.",
    " *",
    " * GENERATED by tools/gen_ccallable.R from highs_c_api.h and",
    " * HighsCallbackStruct.h of the bundled HiGHS (MIT license, Copyright (c)",
    " * HiGHS); do not edit, regenerate.",
    " *",
    " * Usage in the consumer package:",
    " *",
    " *   DESCRIPTION:  Imports: highs      (so the namespace, hence highs.so, loads)",
    " *                 LinkingTo: highs    (so this header is on the include path)",
    " *   NAMESPACE:    import at least one symbol from highs, or call",
    " *                 loadNamespace(\"highs\") -- this header also does the latter",
    " *                 the first time any function below is used.",
    " *   C or C++:     #include <highs_ccallable.h>",
    " *",
    sprintf(" * Every function of the C API (%d in this HiGHS) is here under its own name", n),
    " * and signature, as a static inline wrapper that looks the real function up",
    " * once through R_GetCCallable(\"highs\", \"<name>\") and caches the pointer (per",
    " * translation unit, which is harmless).  Constants and the callback structs",
    " * are copied from the HiGHS headers, so code written against this header",
    " * also compiles against a regular HiGHS installation by including",
    " * highs_c_api.h instead; this header defines the include guards of the",
    " * HiGHS headers it replaces, so including both in that order is harmless.",
    " * HighsInt is int: the R package builds HiGHS without HIGHSINT64",
    " * (Highs_getSizeofHighsInt() reports 4 at run time).",
    if (length(api$skipped)) c(
    " *",
    " * Not here, declared static in highs_c_api.h (internal linkage):",
    sprintf(" *   %s.", paste(api$skipped, collapse = ", "))),
    " */",
    "",
    "#ifndef HIGHS_CCALLABLE_H",
    "#define HIGHS_CCALLABLE_H",
    "",
    "#include <stdint.h>",
    "",
    "#ifndef R_NO_REMAP",
    "#define R_NO_REMAP",
    "#endif",
    "#include <Rinternals.h>",
    "#include <R_ext/Rdynload.h>",
    "",
    "/* ---- util/HighsInt.h (32-bit build) ---------------------------------- */",
    "",
    "#ifndef UTIL_HIGHS_INT_H_",
    "#define UTIL_HIGHS_INT_H_",
    "typedef int HighsInt;",
    "typedef unsigned int HighsUInt;",
    "#define HIGHSINT_FORMAT \"d\"",
    "#endif",
    "",
    "#ifdef __cplusplus",
    "extern \"C\" {",
    "#endif",
    "",
    "/* ---- lp_data/HighsCallbackStruct.h, verbatim -------------------------- */",
    "",
    "#ifndef LP_DATA_HIGHSCALLBACKSTRUCT_H_",
    "#define LP_DATA_HIGHSCALLBACKSTRUCT_H_",
    cb_block,
    "#endif /* LP_DATA_HIGHSCALLBACKSTRUCT_H_ */",
    "",
    "/* ---- interfaces/highs_c_api.h: constants ------------------------------ */",
    "",
    "#ifndef HIGHS_C_API",
    "#define HIGHS_C_API",
    "",
    hi_cc_const_lines(api$consts),
    "",
    "/* ---- resolution through R's registry ---------------------------------- */",
    "",
    "static DL_FUNC highs_ccallable_lookup(const char* name) {",
    "  static int namespace_loaded = 0;",
    "  if (!namespace_loaded) {",
    "    SEXP pkg = PROTECT(Rf_mkString(\"highs\"));",
    "    SEXP call = PROTECT(Rf_lang2(Rf_install(\"loadNamespace\"), pkg));",
    "    Rf_eval(call, R_BaseEnv);",
    "    UNPROTECT(2);",
    "    namespace_loaded = 1;",
    "  }",
    "  return R_GetCCallable(\"highs\", name);",
    "}",
    "",
    "#define HIGHS_CCALLABLE_RESOLVE(name, fn_t) \\",
    "  static fn_t fn = NULL;                    \\",
    "  if (fn == NULL) fn = (fn_t) highs_ccallable_lookup(#name)",
    "",
    sprintf("/* ---- interfaces/highs_c_api.h: the %d functions ---------------------- */", n),
    "",
    unlist(lapply(api$funs, function(f) c(hi_cc_wrapper(f), "")), use.names = FALSE),
    "#undef HIGHS_CCALLABLE_RESOLVE",
    "",
    "#endif /* HIGHS_C_API */",
    "",
    "#ifdef __cplusplus",
    "}",
    "#endif",
    "",
    "#endif /* HIGHS_CCALLABLE_H */")

  list(names = nms, list_h = list_h, header = header)
}

hi_cc_write <- function(root = ".") {
  gen <- hi_cc_generate(root)
  writeLines(gen$list_h, file.path(root, hi_cc_list_out))
  writeLines(gen$header, file.path(root, hi_cc_header_out))
  cat(sprintf("gen_ccallable: %d functions -> %s, %s\n",
              length(gen$names), hi_cc_list_out, hi_cc_header_out))
  invisible(gen)
}

hi_cc_check <- function(root = ".") {
  gen <- hi_cc_generate(root)
  same_list <- identical(readLines(file.path(root, hi_cc_list_out), warn = FALSE), gen$list_h)
  same_header <- identical(readLines(file.path(root, hi_cc_header_out), warn = FALSE), gen$header)
  cat(sprintf("gen_ccallable --check: %d functions; %s %s; %s %s\n", length(gen$names),
              hi_cc_list_out, if (same_list) "up to date" else "STALE",
              hi_cc_header_out, if (same_header) "up to date" else "STALE"))
  same_list && same_header
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--check")) {
    quit(status = if (hi_cc_check()) 0L else 1L)
  } else if (length(args)) {
    stop("usage: Rscript tools/gen_ccallable.R [--check]")
  } else {
    hi_cc_write()
  }
}
