% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver-api.R
\name{hi_solver_set_coeff}
\alias{hi_solver_set_coeff}
\title{Set a coefficient in the constraint matrix.}
\usage{
hi_solver_set_coeff(solver, row, col, val)
}
\arguments{
\item{solver}{An object of class "highs_solver".}

\item{row}{The row index.}

\item{col}{The column index.}

\item{val}{The coefficient value.}
}
\value{
The solver instance with the updated coefficient.
}
\description{
This function assigns a coefficient value to a specific entry in the constraint matrix.
}
\examples{
solver <- example_solver()
hi_solver_set_coeff(solver, 1, 1, 4.2)

}
