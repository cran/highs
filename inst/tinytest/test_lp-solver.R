if (FALSE) {
    library("tinytest")
    library("highs")
}


L <- c(2, 4, 3)
A <- matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE)
rhs <- c(60, 40, 80)
s <- highs_solve(L = L, lower = 0, A = A, rhs = rhs, maximum = TRUE)
tinytest::expect_equal(s[["objective_value"]], 230 / 3)


s <- highs_solve(L = L, lower = 0, maximum = TRUE)
tinytest::expect_equal(s[["primal_solution"]], c(0, 0, 0))


L <- c(3, 1, 3)
A <- rbind(c(-1,  2,  1), c( 0,  4, -3), c( 1, -3,  2))
rhs <- c(4, 2, 3)
lower <- c(-Inf, 0, 2)
upper <- c(4, 100, Inf)
types <- c("I", "C", "I")
s <- highs_solve(L = L, lower = lower, upper = upper, A = A, rhs = rhs,
                 types = types, maximum = TRUE)
tinytest::expect_equal(s[["objective_value"]], 23.5)


L <- c(3, 1, 3)
A <- rbind(c(-1, 2, 1), c(0, 4, -3), c(1, -3, 2))
rhs <- c(4, 2, 3)
types <- c("I", "C", "I")
s <- highs_solve(L = L, lower = 0, A = A, rhs = rhs, types = types, maximum = TRUE)
tinytest::expect_equal(s[["objective_value"]], 26.75)
