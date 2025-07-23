#' @keywords internal
update_rp <- function(...) {
  args <- list2(...)
  if (is.null(attr(args$rp, "family"))) {
    family_string <- paste0(args$family$family, "(", args$family$link, ")")
    attr(args$rp, "family_string") <- family_string
  }
  args$rp
}


#' Constructor Function for Building \code{'randomprojection'} Object
#'
#' Creates an object class \code{'randomprojection'} using arguments passed by user.
#' @param name character
#' @param generate_fun function for generating the random projection matrix. This
#' function should have with arguments \code{rp}, which is a \code{'randomprojection'}
#' object, \code{m}, the target dimension and a vector of indexes
#' \code{included_vector}, \code{x} matrix of predictors and \code{y} matrix of predictors.
#' Vector \code{included_vector} shows the column index of the original variables in the
#' \code{x} matrix to be projected using the random projection. This is needed
#' due to the fact that screening is employed pre-projection.
#' @param update_fun function for updating the \code{'randomprojection'} object with
#' information from the data. This
#' function should have arguments \code{rp}, which is a \code{'randomprojection'}
#' object and `x` (the matrix of predictors)
#' and `y` (the vector of responses).
#' @param update_rpm_w_data function for updating the random projection matrix with data.
#' This can be used for the case where a list of random projection matrices is
#' provided by argument \code{RPMs}. In this case, the random structure is kept
#' fixed, but the data-dependent part gets updated with the provided data. Defaults
#' to NULL. If not provided, the values of the provided RPMs do not change.
#' @param control list of controls for random projection. Can include minimum and
#' maximum dimension for the projection defaults to
#' \code{list(mslow = NULL, msup = NULL)}
#' @return a function which in turn creates an object of class \code{'randomprojection'}
#'
#' @export
constructor_randomprojection <- function(name,
                                         generate_fun,
                                         update_fun = NULL,
                                         update_rpm_w_data = NULL,
                                         control = list()) {
  ## Checks
  stopifnot("Function generate_fun needs arguments rp, m, included_vector, x, y."=
              names(formals(generate_fun)) %in% c("rp", "m", "included_vector", "x", "y"))
  if (!is.null(update_fun)) {
    stopifnot("Function update_fun should have as argument .... All arguments of spar are passed through ..."=names(formals(update_fun)) %in% c("..."))
  } else {
    update_fun <- update_rp
  }
  if (!is.null(update_rpm_w_data)) {
    stopifnot(
      "Function update_rpm_w_data should have arguments rpm,  rp, included_vector."=names(formals(update_rpm_w_data)) %in% c("rpm", "rp", "included_vector"))
  }
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                generate_fun = generate_fun,
                update_fun = update_fun,
                update_rpm_w_data = update_rpm_w_data,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    class(out) <- c("randomprojection")
    return(out)
  }
}


#'
#' Gaussian Random Projection Matrix
#'
#' @param rp object of class  \code{'randomprojection'}
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @param x matrix of predictors
#' @param y vector of response variable
#' @return matrix with m rows and
#'  \code{length(included_vector)} columns sampled from the normal distribution.
#' @keywords internal
generate_gaussian <- function(rp, m, included_vector, x = NULL, y = NULL) {
  p <- length(included_vector)
  control_rnorm <- c(
    rp$control[names(rp$control) %in% names(formals(rnorm))],
    attributes(rp)[names(attributes(rp)) %in% names(formals(rnorm))])
  # remove duplicates
  control_rnorm <-  control_rnorm[!duplicated(names(control_rnorm))]

  vals <- do.call(function(...)
    rnorm(m * p, ...), control_rnorm)
  RM <- matrix(vals, nrow = m, ncol = p)
  return(RM)
}
#'
#' Gaussian Random Projection Matrix
#'
#' @description
#' Creates an object class \code{'randomprojection'} using arguments passed by
#' user which in turn can be employed to generate a random matrix with normally
#' distributed entries (mean 0 and standard deviation 1 by default).
#'
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_fun}, \code{update_rpm_w_data}
#'
#' @return object of class \code{'randomprojection'} which is a list with
#' elements name,
#' \code{generate_fun},  \code{update_fun},  \code{control}
#'
#' @details
#' Arguments related to the random projection procedure can
#' be passed to the \code{rp_gaussian()} function through \code{...}, and
#' will be saved as attributes of the \code{'randomprojection'} object.
#' The following attributes are relevant for [spar] and [spar.cv]:
#'  \itemize{
#'  \item \code{mslow}: integer giving the minimum dimension to which the predictors
#'  should be projected; defaults to \eqn{\log(p)}.
#'  \item \code{msup}: integer giving the maximum dimension to which the predictors
#'  should be projected; defaults to \eqn{n/2}.
#'  }
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   rp = rp_gaussian(control = list(sd = 1/sqrt(ncol(example_data$x)))))
#'
#' @export
#'
rp_gaussian <- constructor_randomprojection(
  "rp_gaussian",
  generate_fun = generate_gaussian
)

#'
#' Sparse Random Projection Matrix
#'
#' @param rp object of class  \code{'randomprojection'}
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @param x matrix of predictors
#' @param y vector of response variable
#' @return (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_sparse <- function(rp, m, included_vector, x = NULL, y = NULL) {
  p <- length(included_vector)
  psi <- rp$control$psi
  if (is.null(psi)) psi <- attr(rp, "psi")
  if (is.null(psi)) psi <- 1
  if (psi > 1 | psi <= 0) stop("For a sparse rpm, psi should lie in interval (0,1].")
  v <- sample(c(-1, 0, 1), size = m * p,
              prob = c(psi/2, 1 - psi, psi/2), replace=TRUE)
  RM <- matrix(v/sqrt(psi), nrow = m, ncol = p)
  RM <- RM[rowSums(abs(RM)) > 0, ]
  RM <- Matrix(RM, sparse = TRUE)
  return(RM)
}

#'
#' Sparse Random Projection Matrix
#'
#' @description
#' Creates an object class \code{'randomprojection'} using arguments passed by
#'  user which in turn can be employed to generate a sparse embedding matrix as
#'  in  \insertCite{ACHLIOPTAS2003JL}{spareg}.
#'
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix.
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_fun}, \code{update_rpm_w_data}
#'
#' @return object of class \code{'randomprojection'} which is a list with
#' elements name,
#' \code{generate_fun},  \code{update_fun},  \code{control}
#'
#' @details
#' The sparse matrix used in \insertCite{ACHLIOPTAS2003JL}{spareg} with entries equal to
#' \eqn{\Psi_{ij} = \pm 1/\sqrt{\psi}} with probability \eqn{\psi/2} and zero otherwise
#' for \eqn{\psi\in (0,1]}. Default is \code{psi = 1}.
#'
#' Arguments related to the random projection procedure can
#' be passed to the \code{rp_gaussian()} function through \code{...}, and
#' will be saved as attributes of the \code{'randomprojection'} object.
#' The following attributes are relevant for [spar] and [spar.cv]:
#'  \itemize{
#'  \item \code{mslow}: integer giving the minimum dimension to which the predictors
#'  should be projected; defaults to \eqn{\log(p)}.
#'  \item \code{msup}: integer giving the maximum dimension to which the predictors
#'  should be projected; defaults to \eqn{n/2}.
#' }
#'
#' @references{
#'   \insertRef{ACHLIOPTAS2003JL}{spareg}
#' }
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   rp = rp_sparse(control = list(psi = 1/3)))
#'
#' @export
#'
rp_sparse <- constructor_randomprojection(
  "rp_sparse",
  generate_fun = generate_sparse
)


#'
#' Sparse Embedding Matrix
#'
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @param x matrix of predictors
#' @param y vector of response variable
#' @return (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_cw <- function(rp, m, included_vector, x = NULL, y = NULL) {
  p <- length(included_vector)
  use_data <- attr(rp, "data")
  if (is.null(use_data)) {
    diagvals <- sample(c(-1, 1), p, replace = TRUE)
  } else {
    if (use_data) {
      if (is.null(attr(rp, "diagvals")))
        stop("Must provide vector of coefficients for data-driven RP.")
      diagvals <- attr(rp, "diagvals")[included_vector]
    } else {
      diagvals <- sample(c(-1, 1), p, replace = TRUE)
    }
  }

  goal_dims <- sample(m, p, replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in seq_len(m)) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims > goal_dim - counter] <-
        goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }
  RM <- Matrix(0, nrow = m - counter, ncol = p,sparse = TRUE)
  RM@i <- as.integer(goal_dims - 1)
  RM@p <- 0:p
  RM@x <- diagvals
  return(RM)
}

update_rp_cw <- function(...) {
  args <- list2(...)
  n <- NROW(args$x)
  p <- NCOL(args$x)
  if (attr(args$screencoef, "reuse_in_rp") &&
      !is.null(attr(args$screencoef, "importance"))) {
    scr_coef <- attr(args$screencoef, "importance")
    inc_probs <- attr(args$screencoef, "inc_prob")
    attr(args$rp, "diagvals") <- scr_coef/max(inc_probs)
  } else {
    if (is.null(args$rp$control$family)) {
      family_string <- paste0(args$family$family, "(", args$family$link, ")")
      args$rp$control$family_string <- family_string
    }
    family <- args$family
    if (family$family=="gaussian" & family$link=="identity") {
      fit_family <- "gaussian"
    } else {
      if (family$family=="binomial" & family$link=="logit") {
        fit_family <- "binomial"
      } else if (family$family=="poisson" & family$link=="log") {
        fit_family <- "poisson"
      } else {
        fit_family <- family
      }
    }
    if (family$family=="gaussian") {
      dev.ratio_cutoff <- 0.999
    } else {
      dev.ratio_cutoff <- 0.8
    }

    if (is.null(args$rp$control$alpha)) args$rp$control$alpha <-  0
    if (is.null(args$rp$control$lambda.min.ratio)) {
      tmp_sc <- apply(args$x, 2, function(col) sqrt(var(col)*(n-1)/n))
      x2 <- scale(args$x, center = colMeans(args$x), scale = tmp_sc)
      ytX <- crossprod(args$y, x2[,tmp_sc > 0])
      lam_max <- 1000 * max(abs(ytX))/n *
        family$mu.eta(family$linkfun(mean(args$y)))/
        family$variance(mean(args$y))
      args$rp$control$lambda.min.ratio <- min(0.01, 1e-4 / lam_max)
    }

    control_glmnet <- args$rp$control[names(args$rp$control)  %in% names(formals(glmnet))]
    glmnet_res <- do.call(function(...)
      glmnet(x = args$x, y = args$y, family = fit_family, ...), control_glmnet)

    lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
    scr_coef <- coef(glmnet_res,s=lam)[-1]
    inc_probs <- abs(scr_coef)
    max_inc_probs <- max(inc_probs)
    attr(args$rp, "diagvals") <- scr_coef/max_inc_probs
  }
  return(args$rp)
}

update_rpm_w_data_cw <- function(rpm, rp, included_vector) {
  rpm@x <-  attr(rp, "diagvals")[included_vector]
  return(rpm)
}



#'
#' Sparse Embedding Matrix
#'
#' @description
#' Creates an object class \code{'randomprojection'} using arguments passed by
#'  user which in turn can be employed to generate a sparse embedding matrix as
#'  in  \insertCite{Clarkson2013LowRankApprox}{spareg}.
#'
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_fun}, \code{update_rpm_w_data}
#'
#' @return object of class \code{'randomprojection'} which is a list with
#' elements name,
#' \code{generate_fun},  \code{update_fun},  \code{control}
#'
#' @details
#' The entries of the matrix are generated based on \insertCite{Clarkson2013LowRankApprox}{spareg}.
#' This matrix is constructed as \eqn{\Phi=BD\in \mathbb{R}^{m\times p}}, where
#' \eqn{B} is a \eqn{(p\times p)} binary matrix, where for each column \eqn{j}
#' an index is uniformly sampled from \eqn{\{1,\ldots,m\}} and the corresponding
#' entry is set to one, and \eqn{D} is a \eqn{(p\times p)} diagonal matrix,
#' with entries \eqn{d_j \sim \text{Unif}(\{-1, 1\})}.
#' If specified as \code{rp_cw(data = TRUE)}, the random elements on the diagonal
#' are replaced by the ridge coefficients with a small penalty, as introduced in
#' \insertCite{parzer2024glms}{spareg}.
#'
#' @references{
#'   \insertRef{Clarkson2013LowRankApprox}{spareg}
#'
#'   \insertRef{parzer2024glms}{spareg}.
#' }
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   rp = rp_cw(data = TRUE))
#'
#' @export
rp_cw <- constructor_randomprojection(
  "rp_cw",
  generate_fun = generate_cw,
  update_fun = update_rp_cw,
  update_rpm_w_data = update_rpm_w_data_cw
)

#' Print Method for a \code{'randomprojection'} Object
#'
#' Print method for a \code{'randomprojection'} object
#' @param x object of class \code{'randomprojection'}
#' @param ... further arguments passed to or from other methods
#' @return text summary
#'
#' @export
print.randomprojection <- function(x, ...) {
  cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  # cat("* Data-dependent:", attr(x,"data"), "\n")
  cat("* Lower bound on goal dimension m:",
      ifelse(is.null(attr(x, "mslow")),
             "not provided, will default to log(p).",
             attr(x, "mslow")), "\n")
  cat("* Upper bound on goal dimension m:",
      ifelse(is.null(attr(x, "msup")),
             "not provided, will default to n/2.",
             attr(x, "mslow")), "\n")
}
