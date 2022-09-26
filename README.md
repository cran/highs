**R** HIGHS Interface
================
Florian Schwendinger</br>
Updated: 2022-09-09

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/highs)](https://CRAN.R-project.org/package=highs)
[![Licence](https://img.shields.io/cran/l/highs)](https://www.gnu.org/licenses/gpl-2.0.en.html)
<!-- badges: end -->

This repository contains an **R** interface to the
[**HiGHS**](https://github.com/ERGO-Code/HiGHS) solver. The
[**HiGHS**](https://github.com/ERGO-Code/HiGHS) solver, is a
**high**-performance open-source **solver** for solving linear
programming (LP), mixed-integer programming (MIP) and quadratic
programming (QP) optimization problems.

# 1 Installation

The package can be installed from
[**CRAN**](https://CRAN.R-project.org/package=highs)

``` r
install.packages("highs")
```

or [**GitLab**](https://gitlab.com/roigrp/solver/highs).

``` r
remotes::install_gitlab("roigrp/solver/highs")
```

### 1.0.1 Using a preinstalled HiGHS library

It is possible to use a precompile HiGHS library by providing the system
variable `R_HIGHS_LIB_DIR`. For example I used

``` sh
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/Z/bin/highslib -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON -DSHARED:bool=OFF -DBUILD_TESTING:bool=OFF
make install
```

to install the **HiGHS** library to `/Z/bin/highslib`

``` r
Sys.setenv(R_HIGHS_LIB_DIR = "/Z/bin/highslib")
install.packages("highs")
# or 
# remotes::install_gitlab("roigrp/solver/highs")
```

# 2 Basic usage

``` r
library("highs")

args(highs_solve)
#> function (Q = NULL, L, lower, upper, A, lhs, rhs, types, maximum = FALSE, 
#>     offset = 0, control = list(), dry_run = FALSE) 
#> NULL
```

## 2.1 LP

``` r
# Minimize
#  x_0 +  x_1 + 3
# Subject to
#                 x_1 <= 7
#  5 <=   x_0 + 2 x_1 <= 15
#  6 <= 3 x_0 + 2 x_1
#  0 <=   x_0         <= 4
#  1 <=           x_1
A <- rbind(c(0, 1), c(1, 2), c(3, 2))
s <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                 A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                 offset = 3)
str(s)
#> List of 6
#>  $ primal_solution: num [1:2] 0.5 2.25
#>  $ objective_value: num 5.75
#>  $ status         : int 7
#>  $ status_message : chr "Optimal"
#>  $ solver_msg     :List of 6
#>   ..$ value_valid: logi TRUE
#>   ..$ dual_valid : logi TRUE
#>   ..$ col_value  : num [1:2] 0.5 2.25
#>   ..$ col_dual   : num [1:2] 0 0
#>   ..$ row_value  : num [1:3] 2.25 5 6
#>   ..$ row_dual   : num [1:3] 0 0.25 0.25
#>  $ info           :List of 18
#>   ..$ valid                     : logi TRUE
#>   ..$ mip_node_count            : num -1
#>   ..$ simplex_iteration_count   : int 2
#>   ..$ ipm_iteration_count       : int 0
#>   ..$ qp_iteration_count        : int 0
#>   ..$ crossover_iteration_count : int 0
#>   ..$ primal_solution_status    : chr "Feasible"
#>   ..$ dual_solution_status      : chr "Feasible"
#>   ..$ basis_validity            : int 1
#>   ..$ objective_function_value  : num 5.75
#>   ..$ mip_dual_bound            : num 0
#>   ..$ mip_gap                   : num Inf
#>   ..$ num_primal_infeasibilities: int 0
#>   ..$ max_primal_infeasibility  : num 0
#>   ..$ sum_primal_infeasibilities: num 0
#>   ..$ num_dual_infeasibilities  : int 0
#>   ..$ max_dual_infeasibility    : num 0
#>   ..$ sum_dual_infeasibilities  : num 0
```

## 2.2 QP

``` r
# Minimize
#  0.5 x^2 - 2 x + y
# Subject to
#  x <= 3
zero <- .Machine$double.eps * 100
Q <- rbind(c(1, 0), c(0, zero))
L <- c(-2, 1)
A <- t(c(1, 0))

cntrl <- list(log_dev_level = 0L)
s <- highs_solve(Q = Q, L = L, A = A, lhs = 0, rhs = 3, control = cntrl)
str(s)
#> List of 6
#>  $ primal_solution: num [1:2] 2e+00 -1e+07
#>  $ objective_value: num -1e+07
#>  $ status         : int 7
#>  $ status_message : chr "Optimal"
#>  $ solver_msg     :List of 6
#>   ..$ value_valid: logi TRUE
#>   ..$ dual_valid : logi TRUE
#>   ..$ col_value  : num [1:2] 2e+00 -1e+07
#>   ..$ col_dual   : num [1:2] 0 0
#>   ..$ row_value  : num 2
#>   ..$ row_dual   : num 0
#>  $ info           :List of 18
#>   ..$ valid                     : logi TRUE
#>   ..$ mip_node_count            : num -1
#>   ..$ simplex_iteration_count   : int 0
#>   ..$ ipm_iteration_count       : int 0
#>   ..$ qp_iteration_count        : int 5
#>   ..$ crossover_iteration_count : int 0
#>   ..$ primal_solution_status    : chr "Feasible"
#>   ..$ dual_solution_status      : chr "Feasible"
#>   ..$ basis_validity            : int 0
#>   ..$ objective_function_value  : num -1e+07
#>   ..$ mip_dual_bound            : num 0
#>   ..$ mip_gap                   : num Inf
#>   ..$ num_primal_infeasibilities: int 0
#>   ..$ max_primal_infeasibility  : num 0
#>   ..$ sum_primal_infeasibilities: num 0
#>   ..$ num_dual_infeasibilities  : int 0
#>   ..$ max_dual_infeasibility    : num 0
#>   ..$ sum_dual_infeasibilities  : num 0
```

# 3 Options

The function `highs_available_solver_options` lists the available solver
options

``` r
d <- highs_available_solver_options()
d[["option"]] <- sprintf("`%s`", d[["option"]])
knitr::kable(d, row.names = FALSE)
```

| option                                          | type    | category |
|:------------------------------------------------|:--------|:---------|
| `allow_unbounded_or_infeasible`                 | bool    | advanced |
| `allowed_cost_scale_factor`                     | integer | advanced |
| `allowed_matrix_scale_factor`                   | integer | advanced |
| `cost_scale_factor`                             | integer | advanced |
| `dual_simplex_cost_perturbation_multiplier`     | double  | advanced |
| `dual_simplex_pivot_growth_tolerance`           | double  | advanced |
| `dual_steepest_edge_weight_error_tolerance`     | double  | advanced |
| `dual_steepest_edge_weight_log_error_threshold` | double  | advanced |
| `factor_pivot_threshold`                        | double  | advanced |
| `factor_pivot_tolerance`                        | double  | advanced |
| `keep_n_rows`                                   | integer | advanced |
| `less_infeasible_DSE_check`                     | bool    | advanced |
| `less_infeasible_DSE_choose_row`                | bool    | advanced |
| `log_dev_level`                                 | integer | advanced |
| `lp_presolve_requires_basis_postsolve`          | bool    | advanced |
| `max_dual_simplex_cleanup_level`                | integer | advanced |
| `max_dual_simplex_phase1_cleanup_level`         | integer | advanced |
| `mps_parser_type_free`                          | bool    | advanced |
| `no_unnecessary_rebuild_refactor`               | bool    | advanced |
| `presolve_pivot_threshold`                      | double  | advanced |
| `presolve_rule_logging`                         | bool    | advanced |
| `presolve_rule_off`                             | integer | advanced |
| `presolve_substitution_maxfillin`               | integer | advanced |
| `primal_simplex_bound_perturbation_multiplier`  | double  | advanced |
| `rebuild_refactor_solution_error_tolerance`     | double  | advanced |
| `run_crossover`                                 | bool    | advanced |
| `simplex_dualise_strategy`                      | integer | advanced |
| `simplex_initial_condition_check`               | bool    | advanced |
| `simplex_initial_condition_tolerance`           | double  | advanced |
| `simplex_permute_strategy`                      | integer | advanced |
| `simplex_price_strategy`                        | integer | advanced |
| `simplex_unscaled_solution_strategy`            | integer | advanced |
| `start_crossover_tolerance`                     | double  | advanced |
| `use_implied_bounds_from_presolve`              | bool    | advanced |
| `use_original_HFactor_logic`                    | bool    | advanced |
| `dual_feasibility_tolerance`                    | double  | file     |
| `glpsol_cost_row_location`                      | integer | file     |
| `highs_analysis_level`                          | integer | file     |
| `highs_debug_level`                             | integer | file     |
| `infinite_bound`                                | double  | file     |
| `infinite_cost`                                 | double  | file     |
| `ipm_iteration_limit`                           | integer | file     |
| `ipm_optimality_tolerance`                      | double  | file     |
| `large_matrix_value`                            | double  | file     |
| `log_file`                                      | string  | file     |
| `objective_bound`                               | double  | file     |
| `objective_target`                              | double  | file     |
| `primal_feasibility_tolerance`                  | double  | file     |
| `random_seed`                                   | integer | file     |
| `simplex_crash_strategy`                        | integer | file     |
| `simplex_dual_edge_weight_strategy`             | integer | file     |
| `simplex_iteration_limit`                       | integer | file     |
| `simplex_max_concurrency`                       | integer | file     |
| `simplex_min_concurrency`                       | integer | file     |
| `simplex_primal_edge_weight_strategy`           | integer | file     |
| `simplex_scale_strategy`                        | integer | file     |
| `simplex_strategy`                              | integer | file     |
| `simplex_update_limit`                          | integer | file     |
| `small_matrix_value`                            | double  | file     |
| `solution_file`                                 | string  | file     |
| `threads`                                       | integer | file     |
| `write_model_file`                              | string  | file     |
| `write_model_to_file`                           | bool    | file     |
| `write_solution_style`                          | integer | file     |
| `write_solution_to_file`                        | bool    | file     |
| `icrash`                                        | bool    | icrash   |
| `icrash_approx_iter`                            | integer | icrash   |
| `icrash_breakpoints`                            | bool    | icrash   |
| `icrash_dualize`                                | bool    | icrash   |
| `icrash_exact`                                  | bool    | icrash   |
| `icrash_iterations`                             | integer | icrash   |
| `icrash_starting_weight`                        | double  | icrash   |
| `icrash_strategy`                               | string  | icrash   |
| `log_to_console`                                | bool    | logging  |
| `output_flag`                                   | bool    | logging  |
| `mip_abs_gap`                                   | double  | mip      |
| `mip_detect_symmetry`                           | bool    | mip      |
| `mip_feasibility_tolerance`                     | double  | mip      |
| `mip_heuristic_effort`                          | double  | mip      |
| `mip_lp_age_limit`                              | integer | mip      |
| `mip_max_improving_sols`                        | integer | mip      |
| `mip_max_leaves`                                | integer | mip      |
| `mip_max_nodes`                                 | integer | mip      |
| `mip_max_stall_nodes`                           | integer | mip      |
| `mip_min_cliquetable_entries_for_parallelism`   | integer | mip      |
| `mip_pool_age_limit`                            | integer | mip      |
| `mip_pool_soft_limit`                           | integer | mip      |
| `mip_pscost_minreliable`                        | integer | mip      |
| `mip_rel_gap`                                   | double  | mip      |
| `mip_report_level`                              | integer | mip      |
| `parallel`                                      | string  | run-time |
| `presolve`                                      | string  | run-time |
| `ranging`                                       | string  | run-time |
| `solver`                                        | string  | run-time |
| `time_limit`                                    | double  | run-time |

for additional information see the [HiGHS homepage](https://highs.dev/).

# 4 Status codes

HiGHS currently has the following status codes defined in `HConst.h"`.

| enumerator               | status | message                            |
|--------------------------|-------:|------------------------------------|
| `kNotset`                |      0 | `"Not Set"`                        |
| `kLoadError`             |      1 | `"Load error"`                     |
| `kModelError`            |      2 | `"Model error"`                    |
| `kPresolveError`         |      3 | `"Presolve error"`                 |
| `kSolveError`            |      4 | `"Solve error"`                    |
| `kPostsolveError`        |      5 | `"Postsolve error"`                |
| `kModelEmpty`            |      6 | `"Empty"`                          |
| `kOptimal`               |      7 | `"Optimal"`                        |
| `kInfeasible`            |      8 | `"Infeasible"`                     |
| `kUnboundedOrInfeasible` |      9 | `"Primal infeasible or unbounded"` |
| `kUnbounded`             |     10 | `"Unbounded"`                      |
| `kObjectiveBound`        |     11 | `"Bound on objective reached"`     |
| `kObjectiveTarget`       |     12 | `"Target for objective reached"`   |
| `kTimeLimit`             |     13 | `"Time limit reached"`             |
| `kIterationLimit`        |     14 | `"Iteration limit reached"`        |
| `kUnknown`               |     15 | `"Unknown"`                        |
| `kMin`                   |      0 | `"Not Set"`                        |
| `kMax`                   |     15 | `"Unknown"`                        |
