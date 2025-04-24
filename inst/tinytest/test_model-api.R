library(tinytest)
library(checkmate)

# Test model creation
test_hi_new_model <- function() {
    model <- hi_new_model()
    expect_class(model, "highs_model")
    expect_class(model, "externalptr")
}

# Test setting number of columns
test_hi_model_set_ncol <- function() {
    model <- hi_new_model()
    
    # Test valid input
    expect_silent(hi_model_set_ncol(model, 5L))
    
    # Test invalid inputs
    expect_error(hi_model_set_ncol(model, -1L))  # negative number
    expect_error(hi_model_set_ncol(model, 1.5))  # non-integer
    expect_error(hi_model_set_ncol(model, c(1L, 2L)))  # wrong length
    expect_error(hi_model_set_ncol("not_a_model", 5L))  # wrong class
}

# Test setting number of rows
test_hi_model_set_nrow <- function() {
    model <- hi_new_model()
    
    # Test valid input
    expect_silent(hi_model_set_nrow(model, 3L))
    
    # Test invalid inputs
    expect_error(hi_model_set_nrow(model, -1L))  # negative number
    expect_error(hi_model_set_nrow(model, 1.5))  # non-integer
    expect_error(hi_model_set_nrow(model, c(1L, 2L)))  # wrong length
    expect_error(hi_model_set_nrow("not_a_model", 3L))  # wrong class
}

# Test setting optimization sense
test_hi_model_set_sense <- function() {
    model <- hi_new_model()
    
    # Test valid inputs
    expect_silent(hi_model_set_sense(model, TRUE))   # maximization
    expect_silent(hi_model_set_sense(model, FALSE))  # minimization
    
    # Test invalid inputs
    expect_error(hi_model_set_sense(model, NA))  # NA not allowed
    expect_error(hi_model_set_sense(model, c(TRUE, FALSE)))  # wrong length
    expect_error(hi_model_set_sense("not_a_model", TRUE))  # wrong class
}

# Test setting offset
test_hi_model_set_offset <- function() {
    model <- hi_new_model()
    
    # Test valid inputs
    expect_silent(hi_model_set_offset(model, 0))
    expect_silent(hi_model_set_offset(model, -1.5))
    expect_silent(hi_model_set_offset(model, 100))
    
    # Test invalid inputs
    expect_error(hi_model_set_offset(model, c(1, 2)))  # wrong length
    expect_error(hi_model_set_offset(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_offset("not_a_model", 0))  # wrong class
}

# Test setting objective
test_hi_model_set_objective <- function() {
    model <- hi_new_model()
    hi_model_set_ncol(model, 3L)
    
    # Test valid inputs
    expect_silent(hi_model_set_objective(model, c(1, 2, 3)))
    expect_silent(hi_model_set_objective(model, c(-1.5, 0, 1.5)))
    
    # Test invalid inputs
    expect_error(hi_model_set_objective(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_objective("not_a_model", c(1, 2, 3)))  # wrong class
}

# Test setting bounds
test_hi_model_set_bounds <- function() {
    model <- hi_new_model()
    hi_model_set_ncol(model, 3L)
    
    # Test valid inputs for lower bounds
    expect_silent(hi_model_set_lower(model, c(0, 0, 0)))
    expect_silent(hi_model_set_lower(model, c(-Inf, -1, 0)))
    
    # Test valid inputs for upper bounds
    expect_silent(hi_model_set_upper(model, c(1, 1, 1)))
    expect_silent(hi_model_set_upper(model, c(Inf, 10, 100)))
    
    # Test invalid inputs
    expect_error(hi_model_set_lower(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_upper(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_lower("not_a_model", c(0, 0, 0)))  # wrong class
    expect_error(hi_model_set_upper("not_a_model", c(1, 1, 1)))  # wrong class
}

# Test setting constraint matrix
test_hi_model_set_constraint_matrix <- function() {
    model <- hi_new_model()
    hi_model_set_ncol(model, 2L)
    hi_model_set_nrow(model, 2L)
    
    # Test valid inputs
    matrix1 <- matrix(c(1, 0, 0, 1), nrow = 2)
    expect_silent(hi_model_set_constraint_matrix(model, matrix1))
    
    matrix2 <- matrix(c(-1, 2, 3, -4), nrow = 2)
    expect_silent(hi_model_set_constraint_matrix(model, matrix2))
    
    # Test invalid inputs
    expect_error(hi_model_set_constraint_matrix(model, matrix(NA_real_)))  # NA not allowed
    expect_error(hi_model_set_constraint_matrix("not_a_model", matrix1))  # wrong class
}

# Test setting LHS/RHS
test_hi_model_set_sides <- function() {
    model <- hi_new_model()
    hi_model_set_nrow(model, 2L)
    
    # Test valid inputs
    expect_silent(hi_model_set_lhs(model, c(-Inf, -Inf)))
    expect_silent(hi_model_set_rhs(model, c(10, 20)))
    
    # Test invalid inputs
    expect_error(hi_model_set_lhs(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_rhs(model, NA_real_))  # NA not allowed
    expect_error(hi_model_set_lhs("not_a_model", c(-Inf, -Inf)))  # wrong class
    expect_error(hi_model_set_rhs("not_a_model", c(10, 20)))  # wrong class
}

# Run all tests
test_hi_new_model()
test_hi_model_set_ncol()
test_hi_model_set_nrow()
test_hi_model_set_sense()
test_hi_model_set_offset()
test_hi_model_set_objective()
test_hi_model_set_bounds()
test_hi_model_set_constraint_matrix()
test_hi_model_set_sides()
