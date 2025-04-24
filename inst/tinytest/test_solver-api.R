library(tinytest)
library(checkmate)

# Helper function to create a simple test model
create_test_model <- function() {
    model <- new_model()
    model_add_variables(model, 2, c(0, 0), c(10, 10))
    model_add_constraints(model, 1, c(-Inf), c(5))
    model_add_coefficients(model, 1, c(1, 2), c(1, 1))
    model
}

# Test new_solver creation
test_new_solver <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    expect_class(solver, "externalptr")
    expect_true(!is.null(solver))
}

# Test optimization sense functions
test_solver_sense <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test default sense (should be minimization - FALSE)
    expect_false(solver_get_sense(solver))
    
    # Test setting to maximization
    solver_set_sense(solver, TRUE)
    expect_true(solver_get_sense(solver))
    
    # Test setting back to minimization
    solver_set_sense(solver, FALSE)
    expect_false(solver_get_sense(solver))
}

# Test objective offset
test_solver_offset <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting offset
    offset_value <- 5.5
    solver_set_offset(solver, offset_value)
    # Note: We can't test getting the offset as there's no getter function
    # This test just ensures the function call doesn't error
    expect_true(TRUE)
}

# Test integrality settings
test_solver_integrality <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting integrality for first variable
    expect_silent(solver_set_integrality(solver, 0, 1))  # 1 for integer
    expect_silent(solver_set_integrality(solver, 1, 0))  # 0 for continuous
}

# Test objective coefficient setting
test_solver_objective <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting objective coefficients
    expect_silent(solver_set_objective(solver, 0, 2.5))
    expect_silent(solver_set_objective(solver, 1, -1.0))
}

# Test variable bounds
test_solver_variable_bounds <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting variable bounds
    expect_silent(solver_set_variable_bounds(solver, 0, 1.0, 5.0))
    expect_silent(solver_set_variable_bounds(solver, 1, 0.0, 10.0))
    
    # Test invalid bounds should error (lower > upper)
    expect_error(solver_set_variable_bounds(solver, 0, 5.0, 1.0))
}

# Test constraint bounds
test_solver_constraint_bounds <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting constraint bounds
    expect_silent(solver_set_constraint_bounds(solver, 0, -Inf, 10.0))
    expect_silent(solver_set_constraint_bounds(solver, 0, 0.0, Inf))
}

# Test coefficient setting
test_solver_coefficients <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test setting coefficients
    expect_silent(solver_set_coeff(solver, 0, 0, 2.5))
    expect_silent(solver_set_coeff(solver, 0, 1, 1.5))
}

# Test model writing
test_solver_write_model <- function() {
    model <- create_test_model()
    solver <- hi_new_solver(model)
    
    # Test writing model to file
    temp_file <- tempfile(fileext = ".mps")
    expect_silent(solver_write_model(solver, temp_file))
    expect_true(file.exists(temp_file))
    
    # Cleanup
    unlink(temp_file)
}


# Run all tests
test_new_solver()
test_solver_sense()
test_solver_offset()
test_solver_integrality()
test_solver_objective()
test_solver_variable_bounds()
test_solver_constraint_bounds()
test_solver_coefficients()
test_solver_write_model()
