% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparegcv.R
\name{plot.spar.cv}
\alias{plot.spar.cv}
\title{plot.spar.cv}
\usage{
\method{plot}{spar.cv}(
  x,
  plot_type = c("Val_Measure", "Val_numAct", "res-vs-fitted", "coefs"),
  plot_along = c("nu", "nummod"),
  nummod = NULL,
  nu = NULL,
  xfit = NULL,
  yfit = NULL,
  opt_par = c("best", "1se"),
  prange = NULL,
  coef_order = NULL,
  digits = 2,
  ...
)
}
\arguments{
\item{x}{result of \link{spar.cv} function of class  \code{'spar.cv'}.}

\item{plot_type}{one of  \code{c("Val_Measure","Val_numAct","res-vs-fitted","coefs")}.}

\item{plot_along}{one of  \code{c("nu","nummod")}; ignored when  \code{plot_type="res-vs-fitted"}.}

\item{nummod}{fixed value for  \code{nummod} when  \code{plot_along="nu"} for
\code{plot_type="Val_Measure"} or  \code{"Val_numAct"};
same as for \code{\link{predict.spar.cv}} when plot_type="res-vs-fitted".}

\item{nu}{fixed value for \eqn{\nu} when  \code{plot_along="nummod"}
for  \code{plot_type="Val_Measure"} or  \code{"Val_numAct"}; same as for \code{\link{predict.spar.cv}} when  \code{plot_type="res-vs-fitted"}.}

\item{xfit}{data used for predictions in  \code{"res-vs-fitted"}.}

\item{yfit}{data used for predictions in  \code{"res-vs-fitted"}.}

\item{opt_par}{one of  \code{c("best","1se")}, only needed for
\code{plot_type="res-vs-fitted"} to set type of predictions, see \code{\link{predict.spar.cv}}.}

\item{prange}{optional vector of length 2 for  \code{"coefs"}-plot to give the limits of the predictors' plot range; defaults to  \code{c(1, p)}.}

\item{coef_order}{optional index vector of length p for \code{"coefs"}-plot to give the order of the predictors; defaults to  \code{1 : p}.}

\item{digits}{number of significant digits to be displayed in the axis; defaults to 2L.}

\item{...}{further arguments passed to or from other methods}
}
\value{
\code{'\link[ggplot2:ggplot]{ggplot2::ggplot}'}  object
}
\description{
Plot cross-validation measure or number of active variables over different thresholds or number
of models of \code{'spar.cv'} object, produce a residuals vs fitted plot,
or a plot of the estimated coefficients in each marginal model, sorted by their absolute value.
}
