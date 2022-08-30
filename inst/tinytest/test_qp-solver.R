if (FALSE) {
    library("highs")
}


# Test unconstrained QP
# minimize -x_2 - 3x_3 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
L <- c(0, -1, -3)
Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
str(s <- highs_solve(Q = Q, L = L))
tinytest::expect_equal(s[["objective_value"]], -5.5)

# Test constrained QP
# minimize -x_2 - 3x_3 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
# subject to x_1 + x_3 <= 2; x>=0
A <- rbind(c(1, 0, 1))
str(s <- highs_solve(Q = Q, L = L, lower = 0, A = A, rhs = 2))
tinytest::expect_equal(s[["objective_value"]], -5.25)