#' Create a new solver instance.
#'
#' This function initializes a new Highs solver instance using the provided model pointer.
#'
#' @param model an object of class \code{"highs_model"}
#' 
#' @return A new solver instance.
#' 
#' @examples
#' model <- example_model()
#' solver <- hi_new_solver(model)
#' @export 
hi_new_solver <- function(model) {
    checkmate::assert_class(model, classes = "highs_model")
    solver <- new_solver(model)
    class(solver) <- c("highs_solver", class(solver))
    solver
}


#' Get the optimization sense of the solver instance.
#'
#' This function returns the optimization sense (e.g., minimization or maximization) of the provided solver instance.
#'
#' @param solver An object of class "highs_solver" representing the solver instance.
#' 
#' @return The optimization sense of the solver instance.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_get_sense(solver)
#' 
#' @export 
hi_solver_get_sense <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_sense(solver)
}

#' Set the optimization sense of the solver instance.
#'
#' This function updates the optimization sense for the given solver instance. Use TRUE for maximization and FALSE for minimization.
#'
#' @param solver An object of class "highs_solver".
#' @param maximum A boolean indicating whether to set maximization (TRUE) or minimization (FALSE).
#' 
#' @return The updated solver instance with the new optimization sense.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_sense(solver, TRUE)
#' 
#' @export
hi_solver_set_sense <- function(solver, maximum) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_logical(maximum, len = 1, any.missing = FALSE)
    solver_set_sense(solver, maximum)
}


#' Set the objective offset for the solver.
#'
#' This function sets the objective offset in the solver instance.
#'
#' @param solver An object of class "highs_solver".
#' @param ext_offset A numeric value representing the offset.
#' 
#' @return The solver instance with the updated offset.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_offset(solver, 5.0)
#' 
#' @export
hi_solver_set_offset <- function(solver, ext_offset) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_numeric(ext_offset, len = 1, any.missing = FALSE)
    solver_set_offset(solver, as.double(ext_offset))
}


#' Set integrality for a set of variables in the solver.
#'
#' This function defines whether a variable is categorized as integral or continuous.
#'
#' @param solver An object of class "highs_solver".
#' @param index An integer vector specifying the variable index.
#' @param type An integer vector representing the integrality type.
#' 
#' @return The solver instance with updated integrality settings.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_integrality(solver, 1, 1)
#' 
#' @export
hi_solver_set_integrality <- function(solver, index, type) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_integerish(type, lower = 0, any.missing = FALSE, len = length(index))
    solver_set_integrality(solver, as.integer(index), as.integer(type))
}


#' Set the objective coefficient for a set of variables.
#'
#' This function assigns a coefficient to a variable in the objective function.
#'
#' @param solver An object of class "highs_solver".
#' @param index The variable index.
#' @param coeff A numeric value representing the objective coefficient.
#' 
#' @return The solver instance with the updated objective.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_objective(solver, 2, 3.5)
#' 
#' @export
hi_solver_set_objective <- function(solver, index, coeff) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(coeff, len = length(index), any.missing = FALSE)
    solver_set_objective(solver, as.integer(index), as.double(coeff))
}


#' Set variable bounds for a set of variables.
#'
#' This function sets the lower and upper bounds for a set of variables.
#'
#' @param solver An object of class "highs_solver".
#' @param index The variable index.
#' @param lower The lower bound.
#' @param upper The upper bound.
#' 
#' @return The solver instance with updated variable bounds.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_variable_bounds(solver, 2, 0, 10)
#' 
#' @export
hi_solver_set_variable_bounds <- function(solver, index, lower, upper) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(lower, len = length(index), any.missing = FALSE)
    checkmate::assert_numeric(upper, len = length(index), any.missing = FALSE)
    solver_set_variable_bounds(solver, as.integer(index), as.double(lower), as.double(upper))
}


#' Set constraint bounds for a given constraint.
#'
#' This function sets the lower and upper bounds for a specific constraint.
#'
#' @param solver An object of class "highs_solver".
#' @param index The constraint index.
#' @param lower The lower bound.
#' @param upper The upper bound.
#' 
#' @return The solver instance with updated constraint bounds.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_constraint_bounds(solver, 1, -Inf, 100)
#' 
#' @export
hi_solver_set_constraint_bounds <- function(solver, index, lower, upper) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, len = 1)
    checkmate::assert_numeric(lower, len = 1, any.missing = FALSE)
    checkmate::assert_numeric(upper, len = 1, any.missing = FALSE)
    solver_set_constraint_bounds(solver, as.integer(index), as.double(lower), as.double(upper))
}


#' Set a coefficient in the constraint matrix.
#'
#' This function assigns a coefficient value to a specific entry in the constraint matrix.
#'
#' @param solver An object of class "highs_solver".
#' @param row The row index.
#' @param col The column index.
#' @param val The coefficient value.
#' 
#' @return The solver instance with the updated coefficient.
#' 
#' @examples
#' solver <- example_solver()
#' hi_solver_set_coeff(solver, 1, 1, 4.2)
#' 
#' @export
hi_solver_set_coeff <- function(solver, row, col, val) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(row, lower = 0, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(col, lower = 0, any.missing = FALSE, len = 1)
    checkmate::assert_numeric(val, len = 1, any.missing = FALSE)
    solver_set_coeff(solver, as.integer(row), as.integer(col), as.double(val))
}



#' Add Variables to the Solver
#'
#' This function adds new variables to the solver with specified bounds.
#'
#' @param solver An object of class "highs_solver".
#' @param lower A numeric vector of lower bounds for the new variables.
#' @param upper A numeric vector of upper bounds for the new variables.
#'
#' @return The solver instance with the new variables added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_vars(solver, lower = c(0, 0, 0), upper = c(10, 10, 10))
#'
#' @export
hi_solver_add_vars <- function(solver, lower, upper) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_numeric(lower, any.missing = FALSE)
    checkmate::assert_numeric(upper, any.missing = FALSE)
    checkmate::assert_true(length(lower) == length(upper))
    solver_add_vars(solver, as.double(lower), as.double(upper))
}


#' Clear All Solver Data
#'
#' This function clears all data from the solver instance, including the model and solution.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The cleared solver instance.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear(solver)
#'
#' @export
hi_solver_clear <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_clear(solver)
}


#' Clear Model Data
#'
#' This function clears only the model data from the solver instance.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The solver instance with cleared model data.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear_model(solver)
#'
#' @export
hi_solver_clear_model <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_clear_model(solver)
}


#' Clear Solver State
#'
#' This function clears the internal solver state while preserving the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The solver instance with cleared solver state.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear_solver(solver)
#'
#' @export
hi_solver_clear_solver <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_clear_solver(solver)
}


#' Run the Solver
#'
#' This function executes the optimization solver on the current model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The solver instance after optimization.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#'
#' @export
hi_solver_run <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_run(solver)
}


# Get Model from Solver
#
# This function retrieves the current optimization model from the solver.
#
# @param solver An object of class "highs_solver".
#
# @return The current optimization model.
#
# @examples
# solver <- example_solver()
# model <- hi_solver_get_model(solver)
#
hi_solver_get_model <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_model(solver)
}


#' Get Number of Variables
#'
#' This function returns the number of variables (columns) in the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer representing the number of variables.
#'
#' @examples
#' solver <- example_solver()
#' n_vars <- hi_solver_get_num_col(solver)
#'
#' @export
hi_solver_get_num_col <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_num_col(solver)
}


#' Get Number of Constraints
#'
#' This function returns the number of constraints (rows) in the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer representing the number of constraints.
#'
#' @examples
#' solver <- example_solver()
#' n_constraints <- hi_solver_get_num_row(solver)
#'
#' @export
hi_solver_get_num_row <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_num_row(solver)
}


#' Write Model to File
#'
#' This function writes the current optimization model to a file.
#'
#' @param solver An object of class "highs_solver".
#' @param filename A character string specifying the output file path.
#'
#' @return Invisible NULL.
#'
#' @examples
#' solver <- example_solver()
#' model_file <- tempfile(fileext = ".mps")
#' hi_solver_write_model(solver, model_file)
#'
#' @export
hi_solver_write_model <- function(solver, filename) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_write_model(solver, filename)
}


#' Write Basis to File
#'
#' This function writes the current basis information to a file.
#'
#' @param solver An object of class "highs_solver".
#' @param filename A character string specifying the output file path.
#'
#' @return Invisible NULL.
#'
#' @examples
#' solver <- example_solver()
#' basis_file <- tempfile(fileext = ".txt")
#' hi_solver_write_basis(solver, basis_file)
#'
#' @export
hi_solver_write_basis <- function(solver, filename) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_write_basis(solver, filename)
}


#' Get Solver Status Message
#'
#' This function returns a human-readable message describing the current solver status.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A character string containing the status message.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' message <- hi_solver_status_message(solver)
#'
#' @export
hi_solver_status_message <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_status_message(solver)
}


#' Get Solver Status
#'
#' This function returns the current status of the solver.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A status code indicating the solver state.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' status <- hi_solver_status(solver)
#'
#' @export
hi_solver_status <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_status(solver)
}


#' Get Solver Infinity Value
#'
#' This function returns the value that the solver uses to represent infinity.
#'
#' @return A numeric value representing infinity in the solver.
#'
#' @examples
#' inf <- hi_solver_infinity()
#'
#' @export
hi_solver_infinity <- function() {
    solver_infinity()
}


#' Reset Global Scheduler
#'
#' This function resets the global scheduler used by the solver.
#'
#' @param blocking A logical value indicating whether to wait for completion.
#'
#' @return Invisible NULL.
#'
#' @examples
#' hi_reset_global_scheduler(TRUE)
#'
#' @export
hi_reset_global_scheduler <- function(blocking) {
    checkmate::assert_logical(blocking, len = 1, any.missing = FALSE)
    reset_global_scheduler(blocking)
}


#' Get Solver Information
#'
#' This function retrieves detailed information about the solver's state and performance.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing solver information.
#'
#' @examples
#' solver <- example_solver()
#' info <- hi_solver_info(solver)
#'
#' @export
hi_solver_info <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_info(solver)
}


#' Get Solution
#'
#' This function retrieves the solution from the solver after optimization.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing the solution information.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' solution <- hi_solver_get_solution(solver)
#'
#' @export
hi_solver_get_solution <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_solution(solver)
}


#' Get Boolean Option Value
#'
#' This function retrieves the value of a boolean solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A logical value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_bool_option(solver, "mip_detect_symmetry")
#'
#' @export
hi_solver_get_bool_option <- function(solver, key) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_character(key, len = 1, any.missing = FALSE)
    solver_get_bool_option(solver, key)
}


#' Get Integer Option Value
#'
#' This function retrieves the value of an integer solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return An integer value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_int_option(solver, "log_dev_level")
#'
#' @export
hi_solver_get_int_option <- function(solver, key) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_character(key, len = 1, any.missing = FALSE)
    solver_get_int_option(solver, key)
}


#' Get Double Option Value
#'
#' This function retrieves the value of a double precision solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A numeric value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_dbl_option(solver, "time_limit")
#'
#' @export
hi_solver_get_dbl_option <- function(solver, key) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_character(key, len = 1, any.missing = FALSE)
    solver_get_dbl_option(solver, key)
}


#' Get String Option Value
#'
#' This function retrieves the value of a string solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A character string.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_str_option(solver, "solver")
#'
#' @export
hi_solver_get_str_option <- function(solver, key) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_character(key, len = 1, any.missing = FALSE)
    solver_get_str_option(solver, key)
}


#' Change Variable Bounds
#'
#' This function updates the bounds of an existing variable in the model.
#'
#' @param solver An object of class "highs_solver".
#' @param idx An integer specifying the variable index.
#' @param lower The new lower bound.
#' @param upper The new upper bound.
#'
#' @return The solver instance with updated bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_change_variable_bounds(solver, 1, 0, 10)
#'
#' @export
hi_solver_change_variable_bounds <- function(solver, idx, lower, upper) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(idx, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(lower, len = length(idx), any.missing = FALSE)
    checkmate::assert_numeric(upper, len = length(idx), any.missing = FALSE)
    solver_change_variable_bounds(solver, as.integer(idx), as.double(lower), as.double(upper))
}


#' Change Constraint Bounds
#'
#' This function updates the bounds of an existing constraint in the model.
#'
#' @param solver An object of class "highs_solver".
#' @param idx An integer vector specifying the constraint indices.
#' @param lhs The new left-hand side bound.
#' @param rhs The new right-hand side bound.
#'
#' @return The solver instance with updated constraint bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_change_constraint_bounds(solver, 1, -Inf, 100)
#'
#' @export
hi_solver_change_constraint_bounds <- function(solver, idx, lhs, rhs) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_integerish(idx, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(lhs, len = length(idx), any.missing = FALSE)
    checkmate::assert_numeric(rhs, len = length(idx), any.missing = FALSE)
    solver_change_constraint_bounds(solver, as.integer(idx), as.double(lhs), as.double(rhs))
}


#' Add Constraints to Model
#'
#' This function adds new constraints (rows) to the optimization model.
#'
#' @param solver An object of class "highs_solver".
#' @param lhs A numeric vector of left-hand side bounds.
#' @param rhs A numeric vector of right-hand side bounds.
#' @param start An integer vector of starting positions in the sparse matrix.
#' @param index An integer vector of column indices.
#' @param value A numeric vector of coefficient values.
#'
#' @return The solver instance with new constraints added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_rows(solver, c(-Inf), c(10), c(0, 2), c(0, 1), c(1, 2))
#'
#' @export
hi_solver_add_rows <- function(solver, lhs, rhs, start, index, value) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_numeric(lhs, min.len = 1, any.missing = FALSE)
    checkmate::assert_numeric(rhs, len = length(lhs), any.missing = FALSE)
    checkmate::assert_integerish(start, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(value, len = length(index), any.missing = FALSE, min.len = 1)
    solver_add_rows(solver, as.double(lhs), as.double(rhs), as.integer(start), as.integer(index), as.double(value))
}


#' Add Variables to Model
#'
#' This function adds new variables (columns) to the optimization model.
#'
#' @param solver An object of class "highs_solver".
#' @param costs A numeric vector of objective coefficients.
#' @param lower A numeric vector giving the lower bounds of the new variables.
#' @param upper A numeric vector giving the upper bounds of the new variables.
#' @param start An integer vector of starting positions in the sparse matrix.
#' @param index An integer vector of row indices.
#' @param value A numeric vector of coefficient values.
#'
#' @return The solver instance with new variables added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_cols(solver, c(1), c(0), c(10), c(0, 1), c(0), c(2))
#'
#' @export
hi_solver_add_cols <- function(solver, costs, lower, upper, start, index, value) {
    checkmate::assert_class(solver, classes = "highs_solver")
    checkmate::assert_numeric(costs, min.len = 1, any.missing = FALSE)
    checkmate::assert_numeric(lower, len = length(costs), any.missing = FALSE)
    checkmate::assert_numeric(upper, len = length(costs), any.missing = FALSE)
    checkmate::assert_integerish(start, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, min.len = 1)
    checkmate::assert_numeric(value, len = length(index), any.missing = FALSE, min.len = 1)
    solver_add_cols(solver, as.double(costs), as.double(lower), as.double(upper), as.integer(start), as.integer(index), as.double(value))
}


#' Get Objective Coefficients
#'
#' This function retrieves the objective coefficients of the linear program.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A numeric vector of objective coefficients.
#'
#' @examples
#' solver <- example_solver()
#' costs <- hi_solver_get_lp_costs(solver)
#'
#' @export
hi_solver_get_lp_costs <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_lp_costs(solver)
}


#' Get Variable Bounds
#'
#' This function retrieves the bounds for all variables in the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing lower and upper bounds for all variables.
#'
#' @examples
#' solver <- example_solver()
#' bounds <- hi_solver_get_variable_bounds(solver)
#'
#' @export
hi_solver_get_variable_bounds <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_variable_bounds(solver)
}


#' Get Constraint Bounds
#'
#' This function retrieves the bounds for all constraints in the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing lower and upper bounds for all constraints.
#'
#' @examples
#' solver <- example_solver()
#' bounds <- hi_solver_get_constraint_bounds(solver)
#'
#' @export
hi_solver_get_constraint_bounds <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_constraint_bounds(solver)
}


#' Get Constraint Matrix
#'
#' This function retrieves the constraint matrix of the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A sparse matrix representing the constraints.
#'
#' @examples
#' solver <- example_solver()
#' matrix <- hi_solver_get_constraint_matrix(solver)
#'
#' @export
hi_solver_get_constraint_matrix <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_constraint_matrix(solver)
}


#' Get Variable Types
#'
#' This function retrieves the type (continuous, integer, etc.) of all variables.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A character vector of variable types.
#'
#' @examples
#' solver <- example_solver()
#' types <- hi_solver_get_vartype(solver)
#'
#' @export
hi_solver_get_vartype <- function(solver) {
    checkmate::assert_class(solver, classes = "highs_solver")
    solver_get_vartype(solver)
}
