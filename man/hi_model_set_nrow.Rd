% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-api.R
\name{hi_model_set_nrow}
\alias{hi_model_set_nrow}
\title{Set the number of rows in the model}
\usage{
hi_model_set_nrow(model, nrow)
}
\arguments{
\item{model}{an object of class \code{"highs_model"}.}

\item{nrow}{an integer giving the number of rows (variables) to set in the model}
}
\value{
\code{NULL}
}
\description{
This function sets the number of rows in the given model.
}
\examples{
model <- hi_new_model()
hi_model_set_nrow(model, 5L) # Sets the model to have 5 rows
}
