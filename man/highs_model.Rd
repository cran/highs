% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highs.R
\name{highs_model}
\alias{highs_model}
\title{Create a Highs Model}
\usage{
highs_model(
  Q = NULL,
  L,
  lower,
  upper,
  A = NULL,
  lhs = NULL,
  rhs = NULL,
  types = rep.int(1L, length(L)),
  maximum = FALSE,
  offset = 0
)
}
\arguments{
\item{Q}{a numeric symmetric matrix giving the quadratic part of the objective.}

\item{L}{a numeric vector giving the linear part of the objective function.}

\item{lower}{a numeric vector giving the lower bounds of the variables.}

\item{upper}{a numeric vector giving the upper bounds of the variables.}

\item{A}{a numeric matrix giving the linear part of the constraints. Rows are
constraints, and columns are decision variables.}

\item{lhs}{a numeric vector giving the left hand-side of the linear constraints.}

\item{rhs}{a numeric vector giving the right hand-side of the linear constraints.}

\item{types}{a integer vector or character vector giving the variable types.
\code{'C'} or \code{'1'} for continuous,
\code{'I'} or \code{'2'} for integer,
\code{'SC'} or \code{'3'} for semi continuous,
\code{'SI'} or \code{'4'} for semi integer and
\code{'II'} or \code{'5'} for implicit integer.}

\item{maximum}{a logical if \code{TRUE} the solver searches for a maximum,
if \code{FALSE} the solver searches for a minimum.}

\item{offset}{a numeric value giving the offset (default is \code{0}).}
}
\value{
A an object of class \code{highs_model}.
}
\description{
Solve linear and quadratic mixed integer optimization problems.
}
\examples{
library("highs")
# Minimize:
#  x_0 +  x_1 + 3
# Subject to:
#               x_1 <=  7
#  5 <=  x_0 + 2x_1 <= 15
#  6 <= 3x_0 + 2x_1
#  0 <= x_0 <= 4
#  1 <= x_1
A <- rbind(c(0, 1), c(1, 2), c(3, 2))
m <- highs_model(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
                 A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
                 offset = 3)
m

# Minimize:
#  -x_2 - 3x_3 + (1/2) * (2 x_1^2 - 2 x_1x_3 + 0.2 x_2^2 + 2 x_3^2)
# Subject to:
#  x_1 + x_3 <= 2
#  0 <= x
L <- c(0, -1, -3)
Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
A <- cbind(1, 0, 1)
m <- highs_model(Q = Q, L = L, lower = 0, A = A, rhs = 2)
m
}
