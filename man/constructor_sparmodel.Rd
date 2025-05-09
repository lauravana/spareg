% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparmodel_classes.R
\name{constructor_sparmodel}
\alias{constructor_sparmodel}
\title{Constructor function for building \code{'sparmodel'} objects}
\usage{
constructor_sparmodel(name, model_fun, update_fun = NULL)
}
\arguments{
\item{name}{character}

\item{model_fun}{function for estimating the marginal models which returns the
function should have arguments  and   \code{y} (vector of responses -- standardized
for Gaussian family), \code{z} (the matrix of projected predictors) and a
\code{'sparmodel'} \code{object}.}

\item{update_fun}{optional function for updating the \code{'sparmodel'} object
before the
start of the algorithm.}
}
\value{
a function which in turn creates an
object of class \code{'sparmodel'}.
}
\description{
The created function will return a object of class \code{'sparmodel'} which
constitutes of a list.
}
\details{
Creates an object of class \code{'sparmodel'} using arguments passed by user.
}
