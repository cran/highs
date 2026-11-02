**R** HIGHS Interface
================
Florian Schwendinger</br>
Updated: 2026-10-30

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
```

# 2 Package

The **highs** package provides an API similar to **Rglpk** and a low
level API to the **HiGHS** solver. For most users using `highs_solve` as
shown below should be the best choice.

> **Note:** The low-level API is only meant for developers. Most users
> should avoid it and just use the high-level API (`highs_solve`).

The package functions can be grouped into the following categories:

1.  The main function **`highs_solve`**.
2.  Object style wrappers for the model and the solver **`highs_model`**
    and **`highs_solver`**.
3.  Function **`highs_control`** to construct the control object for
    **`highs_solve`**and **`highs_solver`**.
4.  Low-level API wrapper functions to create and modify models
    **`hi_new_model`** and other functions starting with
    **`hi_model_`**.
5.  Low-level API wrapper functions to create and modify solvers
    **`hi_new_solver`** and other functions starting with
    **`hi_solver_`**.
6.  Functions to create example models **`example_model`** and solvers
    **`example_solver`** for the documentation examples.
7.  Function **`highs_available_solver_options`** to get the available
    solver options.
8.  Function **`highs_write_model`** to write the model to a file.

``` r
library("highs")
```

## 2.1 Examples

The the example models and solvers are included to have small examples
available for the manual.

``` r
writeLines(ls("package:highs", pattern = "^example"))
#> example_model
#> example_solver
```

## 2.2 Low level model functions

The low-level model functions allow to create and modify models. More
details and examples can be found in the manual.

``` r
writeLines(ls("package:highs", pattern = "^hi(|_new)_model"))
#> hi_model_get_ncons
#> hi_model_get_nvars
#> hi_model_set_constraint_matrix
#> hi_model_set_hessian
#> hi_model_set_lhs
#> hi_model_set_lower
#> hi_model_set_ncol
#> hi_model_set_nrow
#> hi_model_set_objective
#> hi_model_set_offset
#> hi_model_set_rhs
#> hi_model_set_sense
#> hi_model_set_upper
#> hi_model_set_vartype
#> hi_new_model
```

## 2.3 Low level solver functions

The low-level solver functions allow to create and modify solvers. More
details and examples can be found in the manual.

``` r
writeLines(ls("package:highs", pattern = "^hi(|_new)_solver"))
#> hi_new_solver
#> hi_solver_add_cols
#> hi_solver_add_rows
#> hi_solver_add_vars
#> hi_solver_change_constraint_bounds
#> hi_solver_change_variable_bounds
#> hi_solver_clear
#> hi_solver_clear_basis
#> hi_solver_clear_model
#> hi_solver_clear_solver
#> hi_solver_get_basis
#> hi_solver_get_bool_option
#> hi_solver_get_constraint_bounds
#> hi_solver_get_constraint_matrix
#> hi_solver_get_dbl_option
#> hi_solver_get_dual_ray
#> hi_solver_get_int_option
#> hi_solver_get_lp_costs
#> hi_solver_get_num_col
#> hi_solver_get_num_row
#> hi_solver_get_option
#> hi_solver_get_options
#> hi_solver_get_primal_ray
#> hi_solver_get_ranging
#> hi_solver_get_run_time
#> hi_solver_get_sense
#> hi_solver_get_solution
#> hi_solver_get_str_option
#> hi_solver_get_variable_bounds
#> hi_solver_get_vartype
#> hi_solver_infinity
#> hi_solver_info
#> hi_solver_postsolve
#> hi_solver_presolve
#> hi_solver_read_model
#> hi_solver_run
#> hi_solver_set_basis
#> hi_solver_set_coeff
#> hi_solver_set_constraint_bounds
#> hi_solver_set_integrality
#> hi_solver_set_objective
#> hi_solver_set_offset
#> hi_solver_set_option
#> hi_solver_set_options
#> hi_solver_set_sense
#> hi_solver_set_solution
#> hi_solver_set_sparse_solution
#> hi_solver_set_start
#> hi_solver_set_variable_bounds
#> hi_solver_status
#> hi_solver_status_message
#> hi_solver_version
#> hi_solver_write_basis
#> hi_solver_write_model
```

## 2.4 High level functions

The high level functions allow to work with models and solvers. More
details and examples can be found in the manual.

``` r
args(highs_model)
#> function (Q = NULL, L, lower, upper, A = NULL, lhs = NULL, rhs = NULL, 
#>     types = rep.int(1L, length(L)), maximum = FALSE, offset = 0) 
#> NULL
args(highs_solver)
#> function (model, control = highs_control()) 
#> NULL
args(highs_control)
#> function (threads = 1L, time_limit = Inf, log_to_console = FALSE, 
#>     ...) 
#> NULL
args(highs_write_model)
#> function (model, file) 
#> NULL
```

## 2.5 Main function

The main function `highs_solve`.

``` r
library("highs")

args(highs_solve)
#> function (Q = NULL, L, lower, upper, A = NULL, lhs = NULL, rhs = NULL, 
#>     types = rep.int(1L, length(L)), maximum = FALSE, offset = 0, 
#>     start = NULL, control = highs_control()) 
#> NULL
```

## 2.6 LP

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

## 2.7 QP

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
#>   ..$ qp_iteration_count        : int 3
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

# 3 Additional information

## 3.1 Sparse matrices

The **HiGHs** **C++** library internally supports the matrix formats csc
(compressed sparse column matrix) and csr (compressed Sparse Row array).
The **highs** package currently supports the following matrix classes:

1.  `"matrix"` dense matrices,  
2.  `"dgCMatrix"` compressed sparse column matrix from the **Matrix**
    package,  
3.  `"dgRMatrix"` compressed sparse row matrix from the **Matrix**
    package,  
4.  `"matrix.csc"` compressed sparse column matrix from the **SparseM**
    package,  
5.  `"matrix.csr"` compressed sparse row matrix from the **SparseM**
    package,  
6.  `"simple_triplet_matrix"` coordinate format from the **slam**
    package.

If the constraint matrix `A` is provided as `dgCMatrix`, `dgRMatrix`,
`matrix.csc` or `matrix.csr` the underlying data is directly passed to
**HiGHs** otherwise it is first transformed into the csc format an
afterwards passed to **HiGHs**

``` r
library("Matrix")

A <- rbind(c(0, 1), c(1, 2), c(3, 2))
csc <- as(A, "CsparseMatrix")  # dgCMatrix
s0 <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                  A = csc, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                  offset = 3)

csr <- as(A, "RsparseMatrix")  # dgRMatrix
s1 <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                  A = csr, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                  offset = 3)

library("SparseM")

csc <- as.matrix.csc(A)
s2 <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                  A = csc, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                  offset = 3)

csr <- as.matrix.csr(A)
s3 <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                  A = csr, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                  offset = 3)

library("slam")
stm <- as.simple_triplet_matrix(A)
s4 <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                  A = stm, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                  offset = 3)
```

# 4 Options

The function `highs_available_solver_options` lists the available solver
options

``` r
d <- highs_available_solver_options()
d[["option"]] <- sprintf("`%s`", d[["option"]])
knitr::kable(d, row.names = FALSE)
```

| option                                          | type    | advanced |
|:------------------------------------------------|:--------|---------:|
| `presolve`                                      | string  |        0 |
| `solver`                                        | string  |        0 |
| `parallel`                                      | string  |        0 |
| `threads`                                       | integer |        0 |
| `run_crossover`                                 | string  |        0 |
| `time_limit`                                    | double  |        0 |
| `ranging`                                       | string  |        0 |
| `infinite_cost`                                 | double  |        0 |
| `infinite_bound`                                | double  |        0 |
| `small_matrix_value`                            | double  |        0 |
| `large_matrix_value`                            | double  |        0 |
| `kkt_tolerance`                                 | double  |        0 |
| `primal_feasibility_tolerance`                  | double  |        0 |
| `dual_feasibility_tolerance`                    | double  |        0 |
| `primal_residual_tolerance`                     | double  |        0 |
| `dual_residual_tolerance`                       | double  |        0 |
| `optimality_tolerance`                          | double  |        0 |
| `objective_bound`                               | double  |        0 |
| `objective_target`                              | double  |        0 |
| `random_seed`                                   | integer |        0 |
| `user_objective_scale`                          | integer |        0 |
| `user_bound_scale`                              | integer |        0 |
| `highs_debug_level`                             | integer |        0 |
| `highs_analysis_level`                          | integer |        0 |
| `simplex_strategy`                              | integer |        0 |
| `simplex_scale_strategy`                        | integer |        0 |
| `simplex_crash_strategy`                        | integer |        0 |
| `simplex_dual_edge_weight_strategy`             | integer |        0 |
| `simplex_primal_edge_weight_strategy`           | integer |        0 |
| `simplex_iteration_limit`                       | integer |        0 |
| `simplex_update_limit`                          | integer |        0 |
| `simplex_min_concurrency`                       | integer |        0 |
| `simplex_max_concurrency`                       | integer |        0 |
| `output_flag`                                   | bool    |        0 |
| `log_to_console`                                | bool    |        0 |
| `timeless_log`                                  | bool    |        0 |
| `log_file`                                      | string  |        0 |
| `write_model_to_file`                           | bool    |        0 |
| `write_presolved_model_to_file`                 | bool    |        0 |
| `write_solution_to_file`                        | bool    |        0 |
| `write_solution_style`                          | integer |        0 |
| `glpsol_cost_row_location`                      | integer |        0 |
| `icrash`                                        | bool    |        0 |
| `icrash_dualize`                                | bool    |        0 |
| `icrash_strategy`                               | string  |        0 |
| `icrash_starting_weight`                        | double  |        0 |
| `icrash_iterations`                             | integer |        0 |
| `icrash_approx_iter`                            | integer |        0 |
| `icrash_exact`                                  | bool    |        0 |
| `icrash_breakpoints`                            | bool    |        0 |
| `read_solution_file`                            | string  |        0 |
| `read_basis_file`                               | string  |        0 |
| `write_model_file`                              | string  |        0 |
| `solution_file`                                 | string  |        0 |
| `write_basis_file`                              | string  |        0 |
| `write_presolved_model_file`                    | string  |        0 |
| `write_iis_model_file`                          | string  |        0 |
| `mip_detect_symmetry`                           | bool    |        0 |
| `mip_allow_restart`                             | bool    |        0 |
| `mip_max_nodes`                                 | integer |        0 |
| `mip_max_stall_nodes`                           | integer |        0 |
| `mip_max_start_nodes`                           | integer |        0 |
| `mip_improving_solution_save`                   | bool    |        0 |
| `mip_improving_solution_report_sparse`          | bool    |        0 |
| `mip_improving_solution_file`                   | string  |        0 |
| `mip_root_presolve_only`                        | bool    |        0 |
| `mip_lifting_for_probing`                       | integer |        0 |
| `mip_max_leaves`                                | integer |        0 |
| `mip_max_improving_sols`                        | integer |        0 |
| `mip_lp_age_limit`                              | integer |        0 |
| `mip_pool_age_limit`                            | integer |        0 |
| `mip_pool_soft_limit`                           | integer |        0 |
| `mip_pscost_minreliable`                        | integer |        0 |
| `mip_min_cliquetable_entries_for_parallelism`   | integer |        0 |
| `mip_report_level`                              | integer |        0 |
| `mip_feasibility_tolerance`                     | double  |        0 |
| `mip_heuristic_effort`                          | double  |        0 |
| `mip_heuristic_run_feasibility_jump`            | bool    |        0 |
| `mip_heuristic_run_rins`                        | bool    |        0 |
| `mip_heuristic_run_rens`                        | bool    |        0 |
| `mip_heuristic_run_root_reduced_cost`           | bool    |        0 |
| `mip_heuristic_run_zi_round`                    | bool    |        0 |
| `mip_heuristic_run_shifting`                    | bool    |        0 |
| `mip_allow_cut_separation_at_nodes`             | bool    |        0 |
| `mip_rel_gap`                                   | double  |        0 |
| `mip_abs_gap`                                   | double  |        0 |
| `mip_min_logging_interval`                      | double  |        0 |
| `mip_lp_solver`                                 | string  |        0 |
| `mip_ipm_solver`                                | string  |        0 |
| `ipm_optimality_tolerance`                      | double  |        0 |
| `mip_search_simulate_concurrency`               | bool    |        0 |
| `ipm_iteration_limit`                           | integer |        0 |
| `hipo_system`                                   | string  |        0 |
| `hipo_parallel_type`                            | string  |        0 |
| `hipo_ordering`                                 | string  |        0 |
| `hipo_block_size`                               | integer |        0 |
| `pdlp_iteration_limit`                          | integer |        0 |
| `pdlp_scaling_mode`                             | integer |        0 |
| `pdlp_ruiz_iterations`                          | integer |        0 |
| `pdlp_restart_strategy`                         | integer |        0 |
| `pdlp_cupdlpc_restart_method`                   | integer |        0 |
| `pdlp_step_size_strategy`                       | integer |        0 |
| `pdlp_optimality_tolerance`                     | double  |        0 |
| `qp_allow_hot_start`                            | bool    |        0 |
| `qp_iteration_limit`                            | integer |        0 |
| `qp_nullspace_limit`                            | integer |        0 |
| `qp_regularization_value`                       | double  |        0 |
| `iis_strategy`                                  | integer |        0 |
| `iis_time_limit`                                | double  |        0 |
| `blend_multi_objectives`                        | bool    |        0 |
| `log_dev_level`                                 | integer |        1 |
| `log_githash`                                   | bool    |        1 |
| `solve_relaxation`                              | bool    |        1 |
| `allow_unbounded_or_infeasible`                 | bool    |        1 |
| `use_implied_bounds_from_presolve`              | bool    |        1 |
| `lp_presolve_requires_basis_postsolve`          | bool    |        1 |
| `mps_parser_type_free`                          | bool    |        1 |
| `use_warm_start`                                | bool    |        1 |
| `write_matrix_image`                            | bool    |        1 |
| `write_hessian_image`                           | bool    |        1 |
| `keep_n_rows`                                   | integer |        1 |
| `cost_scale_factor`                             | integer |        1 |
| `allowed_matrix_scale_factor`                   | integer |        1 |
| `allowed_cost_scale_factor`                     | integer |        1 |
| `ipx_dualize_strategy`                          | integer |        1 |
| `simplex_dualize_strategy`                      | integer |        1 |
| `simplex_permute_strategy`                      | integer |        1 |
| `max_dual_simplex_cleanup_level`                | integer |        1 |
| `max_dual_simplex_phase1_cleanup_level`         | integer |        1 |
| `simplex_price_strategy`                        | integer |        1 |
| `simplex_unscaled_solution_strategy`            | integer |        1 |
| `no_unnecessary_rebuild_refactor`               | bool    |        1 |
| `rebuild_refactor_solution_error_tolerance`     | double  |        1 |
| `dual_steepest_edge_weight_error_tolerance`     | double  |        1 |
| `dual_steepest_edge_weight_log_error_threshold` | double  |        1 |
| `dual_simplex_cost_perturbation_multiplier`     | double  |        1 |
| `primal_simplex_bound_perturbation_multiplier`  | double  |        1 |
| `dual_simplex_pivot_growth_tolerance`           | double  |        1 |
| `presolve_pivot_threshold`                      | double  |        1 |
| `presolve_reduction_limit`                      | integer |        1 |
| `restart_presolve_reduction_limit`              | integer |        1 |
| `presolve_rule_off`                             | integer |        1 |
| `presolve_rule_test`                            | integer |        1 |
| `presolve_rule_logging`                         | bool    |        1 |
| `presolve_remove_slacks`                        | bool    |        1 |
| `presolve_substitution_maxfillin`               | integer |        1 |
| `factor_pivot_threshold`                        | double  |        1 |
| `factor_pivot_tolerance`                        | double  |        1 |
| `start_crossover_tolerance`                     | double  |        1 |
| `use_original_HFactor_logic`                    | bool    |        1 |
| `less_infeasible_DSE_check`                     | bool    |        1 |
| `less_infeasible_DSE_choose_row`                | bool    |        1 |
| `run_centring`                                  | bool    |        1 |
| `max_centring_steps`                            | integer |        1 |
| `centring_ratio_tolerance`                      | double  |        1 |

for additional information see the [HiGHS homepage](https://highs.dev/).

# 5 Status codes

HiGHS currently has the following status codes defined in `HConst.h"`.

| enumerator               | status | message                            |
|:-------------------------|-------:|:-----------------------------------|
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
| `kSolutionLimit`         |     16 | `"Solution limit reached"`         |
| `kInterrupt`             |     17 | `"Interrupted by user"`            |
| `kMemoryLimit`           |     18 | `"Memory limit reached"`           |
| `kHighsInterrupt`        |     19 | `"Interrupted by HiGHS"`           |
| `kMin`                   |      0 | `"Not Set"`                        |
| `kMax`                   |     19 | `"Interrupted by HiGHS"`           |
