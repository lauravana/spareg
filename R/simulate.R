#' Simulate Sparse Regression Data
#'
#' Generates synthetic data for sparse linear regression problems.
#' Returns training and test sets along with model parameters.
#'
#' @param n Integer. Number of training samples.
#' @param p Integer. Number of predictors (features).
#' @param ntest Integer. Number of test samples.
#' @param a Integer. Number of non-zero coefficients in the true beta vector. Default is min(100, p/4).
#' @param snr Numeric. Signal-to-noise ratio. Default is 10.
#' @param rho Numeric between 0 and 1. Correlation coefficient among predictors. Default is 0.5.
#' @param mu Numeric. Intercept term (mean of response). Default is 1.
#' @param seed Integer. Random seed for reproducibility. Default is NULL.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{x}{Training design matrix (n x p).}
#'   \item{y}{Training response vector (length n).}
#'   \item{xtest}{Test design matrix (ntest x p).}
#'   \item{ytest}{Test response vector (length ntest).}
#'   \item{mu}{Intercept used in data generation.}
#'   \item{beta}{True coefficient vector (length p).}
#'   \item{sigma2}{Noise variance used in data generation. Equals }
#' }
#'
#' @examples
#' data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' str(data)
#'
#' @export
simulate_spareg_data <- function(n, p, ntest, a = min(100, p/4),
                                 snr = 10, rho = 0.5, mu = 1,
                                 seed = NULL) {
  beta <- numeric(p)
  if (!is.null(seed)) set.seed(seed)
  beta[1:a] <- sample(c(-3:3)[-4],a,replace = TRUE)
  x <- sqrt(rho)* matrix(rep(rnorm((n+ntest),0,1),p),
                        n+ntest,p) +
    sqrt(1-rho)*matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
  bSb <- rho*sum(beta)^2 + (1-rho)*sum(beta^2)
  sigma2 <- bSb/snr
  y <- mu + x %*% beta + rnorm(n + ntest, 0, sqrt(sigma2))

  xtest <- x[-(1:n),]
  x <- x[1:n,]
  ytest <- y[-(1:n)]
  y <- y[1:n]
  example_data <- list(x = x, y = y,
                       xtest = xtest,
                       ytest = ytest,
                       mu = mu,beta = beta,
                       sigma2 = sigma2)
  example_data
}
