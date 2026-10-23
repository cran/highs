# **R** HIGHS Interface

## 1.15.1 (development)
- `highs_available_solver_options()` now enumerates the bundled library's own
  option registry at runtime (`Highs::getOptionName()`/`getOptionType()` over
  all option records via new internal `solver_option_registry()`) instead of a
  static table regex-scraped from `HighsOptions.h` -- it can no longer drift
  from the compiled solver.  New `advanced` column (options registered past
  HiGHS's user-settable set, e.g. `log_dev_level`).  Drops `pdlp_features_off`,
  a struct member the scrape advertised but HiGHS never registers as an
  option -- it was never actually settable.  `R/constants.R` and
  `inst/parse_highs_options.R` are retained for one release as a cross-check
  (`inst/tinytest/test_option-registry.R`), then will be removed
- compliance shim gains a thread-local console mute (`hi_shim_set_quiet`):
  the `Highs` constructor banner is silenced for the enumerating instance
  without `capture.output()`; worker-thread messages are unaffected
- update bundled HiGHS to v1.15.1
- regenerate `R/constants.R` from v1.15.1's `HighsOptions.h` (156 options; new
  ones include `presolve_rule_test`, `qp_allow_hot_start`,
  `mip_search_simulate_concurrency`, `simplex_initial_condition_check`/
  `_tolerance`, `write_matrix_image`, `write_hessian_image`)
- `inst/patch_min.R`: drop the new unconditional `add_subdirectory(extern)`
  (HiGHS >= 1.15 builds a shared `highs_extras` library by default; the R
  package needs nothing from `extern/`)
- `.gitmodules`: record the actual bundled tag (was stale at v1.9.0)
- `hi_new_solver()` banner filter extended for v1.15's second banner line
  (THIRD_PARTY_NOTICES.md)

## 1.14.0-3 (development)
- CRAN output/exit compliance reworked from release-time source rewriting to a
  compile-time shim: `src/hi_shim.h` is force-included into every bundled
  HiGHS translation unit (C and C++), remapping `printf`/`stdout`/`stderr`/
  `abort`/`exit`/... tokens to package-local functions implemented over
  `Rvprintf`/`REvprintf` (`src/hi_shim.cpp`), with worker-thread messages
  queued and drained on the R main thread
- new post-build symbol gate `tools/audit_symbols.R` replicates CRAN's
  compiled-code scan on `libhighs.a` and fails the build on any banned
  symbol; maintainer-side and OPT-IN (`R_HIGHS_SYMBOL_CHECK=strict`, set by
  the Taskfile check targets) -- ordinary user/CRAN installs skip it
- `inst/patch_highs.R` retired to `inst/archive/`; successor
  `inst/patch_min.R` keeps only the structural patches (std::cout rewrite,
  pdlp/check/examples stubs, cmake surgery) with mandatory match-count
  assertions -- pattern drift now stops the script instead of silently
  doing nothing
- solver log routing simplified: the deprecated `setLogCallback` handler is
  removed; HiGHS's own console path now lands in `Rprintf` (sink()- and
  capture.output()-capturable), still gated by `output_flag`/`log_to_console`;
  the 1024-byte callback truncation is gone
- installed-package size handled declaratively via top-level `.Rinstignore`
  (`/HiGHS`) instead of deleting `inst/HiGHS` during configure
- `-DNDEBUG` passed explicitly to the cmake build (previously inherited from
  HiGHS's Release default); bundled CMakeLists' `CMAKE_CXX_STANDARD` fixed to
  17 via patch (HiGHS 1.14 requires C++17)
- `inst/tinytest/test_stdout-stderr.R` activated: as shipped in 1.14.0-2 its
  three test functions were defined but never called (0 expectations ran);
  also fixed its `presolve` option value ("off", not FALSE)
- `hi_new_solver()` strips the HiGHS version banner from captured
  construction output so warnings carry only their real text

## 1.14.0-2 2026-06-08 (reconstructed from the CRAN release; fix was not previously in git)
- solver output made controllable/capturable per CRAN request:
  `R_message_handler` now receives the solver as callback data and prints only
  when `output_flag` is true; `hi_solver_run(solver, verbose = FALSE)` captures
  solver log output and returns an integer status code (-1 error / 0 success /
  1 warning); `highs_solve()` gains `control$log_to_console`
- `hi_solver_set_option()` accepts non-character scalars (new `as_bool()`
  coercion); `structure(mat, class = "highs_sparse_matrix")` class-attribute fix
- `static_cast` signed/unsigned cleanups in `src/highs_interface.cpp`
- new regression test `inst/tinytest/test_stdout-stderr.R`

## 1.14 2026-05-18
- update bundled HiGHS to v1.14.0
- regenerate `R/constants.R` from `HighsOptions.h` for v1.14.0
- strip `CITATION.cff` from CMakeLists install target (new in v1.13.1) via `inst/patch_highs.R`
- HiGHS v1.14.0 requires C++17; bump `inst/build_highs.sh` accordingly

## 1.12 2025-11-23
- Thank you to Benjamin Schwendinger for fixing the
  `HighsHashTree.h:855:32: runtime error: index 1 out of bounds for type 'NodePtr [1]'`
  in **HiGHS**
  https://github.com/ERGO-Code/HiGHS/pull/2608/commits/c9cbe0b2a283ff0f2aa70a769609ca816189ed0b

## 1.10 2025-04-20
- update to HiGHS v.10
- improve the low-level API
- improve README

## 0.1-11 2024-08-03
- update to HiGHS v1.7.2
- improve API
