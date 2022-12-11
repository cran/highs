#' @import checkmate
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames
#' @importFrom utils modifyList capture.output tail
#' @useDynLib highs, .registration=TRUE
NULL


highs_globals <- local({
    globals <- list(threads = NA_integer_)
    function(key, value) {
        if (missing(key)) return(globals)
        if (missing(value))
            globals[[key]]
        else
            globals[[key]] <<- value
    }
})


highs_infinity <- function() {
    if (is.null(highs_globals("Inf"))) {
        highs_globals("Inf", solver_infinity())
    }
    highs_globals("Inf")
}


get_number_of_threads <- function() {
    highs_globals("threads")
}


set_number_of_threads <- function(threads) {
    if (!isTRUE(highs_globals("threads") == threads)) {
        reset_global_scheduler(TRUE)
        highs_globals("threads", threads)
    }
}


csc_to_matrix <- function(start, index, value, nrow = max(index + 1L), ncol = length(start) - 1L) {
    stopifnot(length(index) == length(value))
    ind <- index + 1L
    M <- matrix(0, nrow, ncol)
    for (i in seq_along(index)) {
        row_id <- ind[i]
        col_id <- min(which(start >= i) - 1L)
        M[row_id, col_id] <- value[i]
    }
    M
}


#' Solve an Optimization Problems
#'
#' Solve linear and quadratic mixed integer optimization problems.
#'
#' @param Q a numeric symmetric matrix giving the quadratic part of the objective.
#' @param L a numeric vector giving the linear part of the objective function.
#' @param lower a numeric vector giving the lower bounds of the variables.
#' @param upper a numeric vector giving the upper bounds of the variables.
#' @param A a numeric matrix giving the linear part of the constraints. Rows are
#'   constraints, and columns are decision variables.
#' @param lhs a numeric vector giving the left hand-side of the linear constraints.
#' @param rhs a numeric vector giving the right hand-side of the linear constraints.
#' @param types a integer vector or character vector giving the variable types.
#'      \code{'C'} or \code{'1'} for continuous,
#'      \code{'I'} or \code{'2'} for integer,
#'      \code{'SC'} or \code{'3'} for semi continuous,
#'      \code{'SI'} or \code{'4'} for semi integer and
#'      \code{'II'} or \code{'5'} for implicit integer.
#' @param maximum a logical if \code{TRUE} the solver searches for a maximum,
#'                if \code{FALSE} the solver searches for a minimum.
#' @param offset a numeric value giving the offset (default is \code{0}).
#' @param control a list giving additional options for the solver,
#'                see \link{highs_available_solver_options} or the \code{README} file
#'                for a list of all available options.
#' @param dry_run a logical if true only the model is returned.
#'
#' @return A \code{list} containing the result provided by the solver,
#'  containing the following named objects:
#' \item{\code{primal_solution}}{a numeric vector giving the primal solution.}
#' \item{\code{objective_value}}{a numeric giving the objective value.}
#' \item{\code{status}}{an integer giving the status code}
#' \item{\code{status_message}}{a character string giving the status message (explanation of the \code{status_code}).}
#' \item{\code{solver_msg}}{a list giving the original (not canonicalized) solver message.}
#' \item{\code{info}}{a list giving additional information provided by the solver.}
#'
#' Additional information on can be found in the \code{README} file.
#'
#' @examples
#' library("highs")
#' # Minimize:
#' #  x_0 +  x_1 + 3
#' # Subject to:
#' #               x_1 <=  7
#' #  5 <=  x_0 + 2x_1 <= 15
#' #  6 <= 3x_0 + 2x_1
#' #  0 <= x_0 <= 4
#' #  1 <= x_1
#' A <- rbind(c(0, 1), c(1, 2), c(3, 2))
#' s <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
#'                  A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
#'                  offset = 3)
#' s[["objective_value"]]
#' s[["primal_solution"]]
#'
#' # Minimize:
#' #  -x_2 - 3x_3 + (1/2) * (2 x_1^2 - 2 x_1x_3 + 0.2 x_2^2 + 2 x_3^2)
#' # Subject to:
#' #  x_1 + x_3 <= 2
#' #  0 <= x
#' L <- c(0, -1, -3)
#' Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
#' A <- cbind(1, 0, 1)
#' s <- highs_solve(Q = Q, L = L, lower = 0, A = A, rhs = 2)
#' s[["objective_value"]]
#' s[["primal_solution"]]
#' @export
highs_solve <- function(Q = NULL, L, lower, upper, A, lhs, rhs, types, maximum = FALSE,
                        offset = 0, control = list(), dry_run = FALSE) {
    assert_numeric(L, any.missing = FALSE)
    default_control <- list(log_to_console = FALSE, threads = 1L, parallel = "off")
    control <- modifyList(default_control, control)
    checkmate::assert_integerish(control$threads, len = 1L, any.missing = FALSE)
    set_number_of_threads(control$threads)
    control$parallel <- if (control$threads > 1) "on" else "off"

    nvars <- length(L)
    if (missing(A) || NROW(A) == 0L) {
        A <- lhs <- rhs <- NULL
        ncons <- 0L
    } else {
        stopifnot(is.vector(L), !(missing(lhs) & missing(rhs)))
        cscA <- as_csc(A)
        ncons <- cscA[["nrow"]]
    }
    model <- new_model()
    INF <- highs_infinity()
    model_set_ncol(model, nvars)
    model_set_nrow(model, ncons)
    model_set_sense(model, maximum)
    model_set_objective(model, L)
    if (!is.null(Q)) {
        cscQ <- as_csc(Q)
        model_set_hessian(model, format = "square", dim = nvars,
            start = cscQ[["col_ptr"]], index = cscQ[["row_id"]], value = cscQ[["value"]])
    }
    if (missing(types) || length(types) == 0L) {
        types <- rep.int(0L, nvars)
    } else {
        if (is.character(types)) {
            types <- match(types, c("C", "I", "SC", "SI", "II")) - 1L
        } else {
            types <- types - 1L
        }
        assert_integerish(types, lower = 0, upper = 4L, any.missing = FALSE)
        model_set_vartype(model, as.integer(types))
    }
    if (missing(lower) || length(lower) == 0L) {
        lower <- rep.int(-INF, nvars)
    } else if (length(lower) == 1L) {
        lower <- rep.int(lower, nvars)
    }
    if (missing(upper) || length(upper) == 0L) {
        upper <- rep.int(INF, nvars)
    } else if (length(upper) == 1L) {
        upper <- rep.int(upper, nvars)
    }

    lower <- replace(lower, lower == -Inf, -INF)
    upper <- replace(upper, upper ==  Inf, INF)
    model_set_lower(model, lower)
    model_set_upper(model, upper)
    if (ncons > 0L) {
        model_set_constraint_matrix(model, "colwise",
            start = cscA[["col_ptr"]], index = cscA[["row_id"]], value = cscA[["value"]])
        if (missing(lhs) || length(lhs) == 0L) {
            lhs <- rep.int(-INF, ncons)
        }
        if (missing(rhs) || length(rhs) == 0L) {
            rhs <- rep.int(INF, ncons)
        }
        lhs <- replace(lhs, lhs == -Inf, -INF)
        rhs <- replace(rhs, rhs ==  Inf, INF)
        model_set_lhs(model, lhs)
        model_set_rhs(model, rhs)
    }
    if (offset != 0) {
        model_set_offset(model, offset)
    }
    if (dry_run) return(model)
    init_msg <- capture.output(solver <- new_solver(model))
    if (is.null(solver)) {
        stop(paste(tail(init_msg, -3), collapse = "\n"))
    }
    solver_set_options(solver, control)

    run_status <- solver_run(solver)
    status <- solver_status(solver)
    status_message <- solver_status_message(solver)

    solution <- solver_solution(solver)
    info <- solver_info(solver)
    list(primal_solution = solution[["col_value"]], objective_value = info[["objective_function_value"]],
         status = status, status_message = status_message, solver_msg = solution, info = info)
}
