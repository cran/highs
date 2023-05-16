
#' Available Solver Options
#'
#' Reference for the available solver options.
#'
#' @returns A \code{data.frame} containing the available solver options.
#' @examples
#' highs_available_solver_options()
#' @export
highs_available_solver_options <- function() {
    option_names <- rbind(
        c("presolve", "string", "run-time"),
        c("solver", "string", "run-time"),
        c("parallel", "string", "run-time"),
        c("ranging", "string", "run-time"),
        c("time_limit", "double", "run-time"),
        c("infinite_cost", "double", "file"),
        c("infinite_bound", "double", "file"),
        c("small_matrix_value", "double", "file"),
        c("large_matrix_value", "double", "file"),
        c("primal_feasibility_tolerance", "double", "file"),
        c("dual_feasibility_tolerance", "double", "file"),
        c("ipm_optimality_tolerance", "double", "file"),
        c("objective_bound", "double", "file"),
        c("objective_target", "double", "file"),
        c("random_seed", "integer", "file"),
        c("threads", "integer", "file"),
        c("highs_debug_level", "integer", "file"),
        c("highs_analysis_level", "integer", "file"),
        c("simplex_strategy", "integer", "file"),
        c("simplex_scale_strategy", "integer", "file"),
        c("simplex_crash_strategy", "integer", "file"),
        c("simplex_dual_edge_weight_strategy", "integer", "file"),
        c("simplex_primal_edge_weight_strategy", "integer", "file"),
        c("simplex_iteration_limit", "integer", "file"),
        c("simplex_update_limit", "integer", "file"),
        c("simplex_min_concurrency", "integer", "file"),
        c("simplex_max_concurrency", "integer", "file"),
        c("ipm_iteration_limit", "integer", "file"),
        c("write_model_file", "string", "file"),
        c("solution_file", "string", "file"),
        c("log_file", "string", "file"),
        c("write_model_to_file", "bool", "file"),
        c("write_solution_to_file", "bool", "file"),
        c("write_solution_style", "integer", "file"),
        c("glpsol_cost_row_location", "integer", "file"),
        c("output_flag", "bool", "logging"),
        c("log_to_console", "bool", "logging"),
        c("log_dev_level", "integer", "advanced"),
        c("run_crossover", "bool", "advanced"),
        c("allow_unbounded_or_infeasible", "bool", "advanced"),
        c("use_implied_bounds_from_presolve", "bool", "advanced"),
        c("lp_presolve_requires_basis_postsolve", "bool", "advanced"),
        c("mps_parser_type_free", "bool", "advanced"),
        c("keep_n_rows", "integer", "advanced"),
        c("cost_scale_factor", "integer", "advanced"),
        c("allowed_matrix_scale_factor", "integer", "advanced"),
        c("allowed_cost_scale_factor", "integer", "advanced"),
        c("simplex_dualise_strategy", "integer", "advanced"),
        c("simplex_permute_strategy", "integer", "advanced"),
        c("max_dual_simplex_cleanup_level", "integer", "advanced"),
        c("max_dual_simplex_phase1_cleanup_level", "integer", "advanced"),
        c("simplex_price_strategy", "integer", "advanced"),
        c("simplex_unscaled_solution_strategy", "integer", "advanced"),
        c("presolve_substitution_maxfillin", "integer", "advanced"),
        c("presolve_rule_off", "integer", "advanced"),
        c("presolve_rule_logging", "bool", "advanced"),
        c("simplex_initial_condition_check", "bool", "advanced"),
        c("no_unnecessary_rebuild_refactor", "bool", "advanced"),
        c("simplex_initial_condition_tolerance", "double", "advanced"),
        c("rebuild_refactor_solution_error_tolerance", "double", "advanced"),
        c("dual_steepest_edge_weight_error_tolerance", "double", "advanced"),
        c("dual_steepest_edge_weight_log_error_threshold", "double", "advanced"),
        c("dual_simplex_cost_perturbation_multiplier", "double", "advanced"),
        c("primal_simplex_bound_perturbation_multiplier", "double", "advanced"),
        c("dual_simplex_pivot_growth_tolerance", "double", "advanced"),
        c("presolve_pivot_threshold", "double", "advanced"),
        c("factor_pivot_threshold", "double", "advanced"),
        c("factor_pivot_tolerance", "double", "advanced"),
        c("start_crossover_tolerance", "double", "advanced"),
        c("less_infeasible_DSE_check", "bool", "advanced"),
        c("less_infeasible_DSE_choose_row", "bool", "advanced"),
        c("use_original_HFactor_logic", "bool", "advanced"),
        c("icrash", "bool", "icrash"),
        c("icrash_dualize", "bool", "icrash"),
        c("icrash_strategy", "string", "icrash"),
        c("icrash_starting_weight", "double", "icrash"),
        c("icrash_iterations", "integer", "icrash"),
        c("icrash_approx_iter", "integer", "icrash"),
        c("icrash_exact", "bool", "icrash"),
        c("icrash_breakpoints", "bool", "icrash"),
        c("mip_detect_symmetry", "bool", "mip"),
        c("mip_max_nodes", "integer", "mip"),
        c("mip_max_stall_nodes", "integer", "mip"),
        c("mip_max_leaves", "integer", "mip"),
        c("mip_max_improving_sols", "integer", "mip"),
        c("mip_lp_age_limit", "integer", "mip"),
        c("mip_pool_age_limit", "integer", "mip"),
        c("mip_pool_soft_limit", "integer", "mip"),
        c("mip_pscost_minreliable", "integer", "mip"),
        c("mip_min_cliquetable_entries_for_parallelism", "integer", "mip"),
        c("mip_report_level", "integer", "mip"),
        c("mip_feasibility_tolerance", "double", "mip"),
        c("mip_rel_gap", "double", "mip"),
        c("mip_abs_gap", "double", "mip"),
        c("mip_heuristic_effort", "double", "mip")
    )
    colnames(option_names) <- c("option", "type", "category")
    option_names <- as.data.frame(option_names)
    option_names[order(option_names[["category"]], option_names[["option"]], option_names[["type"]]), ]
}

.available__solver__options_ <- function() {
    option_names <- list(
        bool = c("write_solution_to_file", "output_flag", "log_to_console",
                 "run_crossover", "allow_unbounded_or_infeasible",
                 "use_implied_bounds_from_presolve", "lp_presolve_requires_basis_postsolve",
                 "mps_parser_type_free", "simplex_initial_condition_check",
                 "no_unnecessary_rebuild_refactor", "less_infeasible_DSE_check",
                 "less_infeasible_DSE_choose_row", "use_original_HFactor_logic",
                 "mip_detect_symmetry"),
        integer = c("log_dev_level", "random_seed", "threads", "highs_debug_level",
                    "highs_analysis_level", "simplex_strategy", "simplex_scale_strategy",
                    "simplex_crash_strategy", "simplex_dual_edge_weight_strategy",
                    "simplex_primal_edge_weight_strategy", "simplex_iteration_limit",
                    "simplex_update_limit", "simplex_min_concurrency", "simplex_max_concurrency",
                    "ipm_iteration_limit", "write_solution_style", "keep_n_rows",
                    "cost_scale_factor", "allowed_matrix_scale_factor", "allowed_cost_scale_factor",
                    "simplex_dualise_strategy", "simplex_permute_strategy",
                    "max_dual_simplex_cleanup_level", "max_dual_simplex_phase1_cleanup_level",
                    "simplex_price_strategy", "simplex_unscaled_solution_strategy",
                    "presolve_substitution_maxfillin", "mip_max_nodes", "mip_max_stall_nodes",
                    "mip_max_leaves", "mip_lp_age_limit", "mip_pool_age_limit",
                    "mip_pool_soft_limit", "mip_pscost_minreliable", "mip_report_level"),
        double = c("time_limit", "infinite_cost", "infinite_bound", "small_matrix_value",
                   "large_matrix_value", "primal_feasibility_tolerance",
                   "dual_feasibility_tolerance", "ipm_optimality_tolerance",
                   "objective_bound", "objective_target", "simplex_initial_condition_tolerance",
                   "rebuild_refactor_solution_error_tolerance",
                   "dual_steepest_edge_weight_log_error_threshold",
                   "dual_simplex_cost_perturbation_multiplier",
                   "primal_simplex_bound_perturbation_multiplier",
                   "dual_simplex_pivot_growth_tolerance", "presolve_pivot_threshold",
                   "factor_pivot_threshold", "factor_pivot_tolerance", "start_crossover_tolerance",
                   "mip_feasibility_tolerance", "mip_heuristic_effort"),
        string = c("presolve", "solver", "parallel", "ranging", "solution_file", "log_file")
    )
    option_names
}


solver_get_options <- function(solver, keys = NULL) {
    assert_character(keys, null.ok = TRUE)
    option_names <- highs_available_solver_options()
    getters <- list(bool = solver_get_bool_option, integer = solver_get_int_option,
                    double = solver_get_dbl_option, string = solver_get_str_option)
    opts <- vector("list", NROW(option_names))
    names(opts) <- option_names[["option"]]
    for (i in seq_len(NROW(option_names))) {
        row <- option_names[i,]
        solver_get_option <- getters[[row[["type"]]]]
        key <- row[["option"]]
        opts[[key]] <- solver_get_option(solver, key)
    }
    if (length(keys) > 0L) {
        opts[which(names(opts) %in% keys)]
    } else {
        opts
    }
}


force_type <- function(obj, type) {
    force <- list(bool = as.logical, integer = as.integer,
                  double = as.double, string = as.character)[[type]]
    force(obj)
}


test_type <- function(obj, type) {
    if (type == "bool") {
        checkmate::test_logical(obj, len = 1L, any.missing = FALSE)
    } else if (type == "integer") {
        checkmate::test_integer(obj, len = 1L, any.missing = FALSE)
    } else if (type == "double") {
        checkmate::test_double(obj, len = 1L, any.missing = FALSE)
    } else if (type == "string") {
        checkmate::test_string(obj)
    } else {
        stop(sprintf("type '%s' is not allowed", type))
    }
}


solver_set_options <- function(solver, kwargs = list()) {
    if (length(kwargs) == 0) return(invisible(NULL))
    option_names <- highs_available_solver_options()
    assert_character(names(kwargs), min.len = 1L, any.missing = FALSE)
    m <- match(names(kwargs), option_names[["option"]])
    if (anyNA(m)) {
        stop("unknown option")
    }
    option_names <- option_names[m, ]
    option_types <- setNames(option_names[["type"]], option_names[["option"]])
    for (i in seq_along(kwargs)) {
        key <- names(kwargs)[i]
        value <- kwargs[[i]]
        dtype <- option_types[key]
        val <- force_type(value, dtype)
        if (!test_type(val, dtype)) {
            msg <- sprintf("option '%s' expects object of class '%s' got '%s'",
                           key, dtype, paste(class(val), collapse = "."))
            stop(msg)
        }
        solver_set_option(solver, key, val)
    }
    return(invisible(NULL))
}

