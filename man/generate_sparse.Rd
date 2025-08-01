% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rp_classes.R
\name{generate_sparse}
\alias{generate_sparse}
\title{Sparse Random Projection Matrix}
\usage{
generate_sparse(rp, m, included_vector, x = NULL, y = NULL)
}
\arguments{
\item{rp}{object of class  \code{'randomprojection'}}

\item{m}{goal dimension, which will be randomly sampled in the SPAR algorithm}

\item{included_vector}{integer vector of column indices for the variables to be
included in the random projection. These indices are produced in the
screening step of the SPAR algorithm.}

\item{x}{matrix of predictors}

\item{y}{vector of response variable}
}
\value{
(possibly sparse) matrix with m rows and
\code{length(included_vector)} columns.
}
\description{
Sparse Random Projection Matrix
}
\keyword{internal}
