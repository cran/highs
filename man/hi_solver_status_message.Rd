% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver-api.R
\name{hi_solver_status_message}
\alias{hi_solver_status_message}
\title{Get Solver Status Message}
\usage{
hi_solver_status_message(solver)
}
\arguments{
\item{solver}{An object of class "highs_solver".}
}
\value{
A character string containing the status message.
}
\description{
This function returns a human-readable message describing the current solver status.
}
\examples{
solver <- example_solver()
hi_solver_run(solver)
message <- hi_solver_status_message(solver)

}
