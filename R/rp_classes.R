#' Random Projection Object Class
#'
#' @description
#' The `'randomprojection'` class' represents a configuration for generating and managing random projection matrices in the **SPAR (Sparse Projected Averaged Regression)** framework. Objects of this class encapsulate:
#'   - Functions for generating, updating, and modifying random projection matrices.
#'   - Control parameters for the random projection process.
#'   - Metadata (e.g., name, attributes) for customization.
#'
#' @details
#' Objects of class `'randomprojection'` are created using the `constructor_randomprojection` function. They are used in the `spar` and `spar.cv` functions to define how predictors are projected into a lower-dimensional space.
#'
#' The class includes the following components:
#'   - **`name`**: A character string describing the random projection method (e.g., `"rp_gaussian"`, `"rp_sparse"`).
#'   - **`generate_fun`**: A function to generate the random projection matrix. This function must accept arguments like `object`, `x`, `y`, `m`, `included_vector`, and `...`.
#'   - **`update_fun`**: A function to update the `randomprojection` object with data-specific information. This function must accept arguments like `object`, `x`, `y`, `family`, and `...`.
#'   - **`update_rpm_w_data`**: A function to update an already-generated random projection matrix with data-dependent information. This function must accept arguments like `rpm`, `object`, `included_vector`, `x`, `y`, `family`, and `...`.
#'   - **`control`**: A list of control parameters for the random projection process (e.g., `mslow`, `msup`, `psi`).
#'
#' @section Attributes:
#' The following attributes are commonly used in `randomprojection` objects:
#'   - **`mslow`**: Integer. The minimum dimension for projection. Default: `\eqn{\log(p)}`.
#'   - **`msup`**: Integer. The maximum dimension for projection. Default: `\eqn{n/2}`.
#'   - **`data`**: Logical. If `TRUE`, the projection matrix is updated with data-dependent coefficients (e.g., for `rp_cw`).
#'
#' @section Usage:
#' `'randomprojection'` objects are typically created using predefined constructors like:
#'   - `rp_gaussian()`: Gaussian random projection.
#'   - `rp_sparse()`: Sparse random projection (Achlioptas, 2003).
#'   - `rp_cw()`: Sparse Embedding (Clarkson-Woodruff) random projection.
#'
#' Users can also define custom random projection methods by providing their own `generate_fun`, `update_fun`, and `update_rpm_w_data` functions to `constructor_randomprojection`.
#'
#' @examples
#' rp_gauss <- rp_gaussian(mslow = 5, msup = 10)
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(
#'   example_data$x,
#'   example_data$y,
#'   xval = example_data$xtest,
#'   yval = example_data$ytest,
#'   rp = rp_gauss
#' )
#'
#' @seealso
#' [`constructor_randomprojection`], [`rp_gaussian`], [`rp_sparse`], [`rp_cw`], [`spar`]
#'
#' @references{
#'  \insertRef{ACHLIOPTAS2003JL}{spareg}
#'
#'   \insertRef{Clarkson2013LowRankApprox}{spareg}
#'
#'   \insertRef{parzer2024glms}{spareg}.
#' }
#'
#' @name randomprojection-class
NULL

#' Default Update Function for Random Projection Objects
#'
#' @description
#' A default function to update a [`randomprojection-class`] object with the family information.
#' This function is used internally by `constructor_randomprojection` if no custom `update_fun` is provided.
#'
#' @param object An object of class [`randomprojection-class`].
#' @param x A matrix of standardized predictors (unused in this default function).
#' @param y A vector of standardized responses (unused in this default function).
#' @param family A [`stats::family`] object specifying the GLM family and link function.
#' @param ... Additional arguments (ignored).
#'
#' @return
#' The updated \code{'randomprojection'} object with the
#' `family_string` attribute set.
#'
#' @keywords internal
update_rp_default <- function(object, x, y, family, ...) {
  family_string <- paste0(family$family, "(", family$link, ")")
  attr(object, "family_string") <- family_string
  object
}

#' Default Function to Update Random Projection Matrices
#'
#' @description
#' A default function that returns the input random projection matrix unchanged.
#' This function is used internally by `constructor_randomprojection` if no custom `update_rpm_w_data` is provided.
#'
#' @param rpm A random projection matrix.
#' @param object An object of class [`randomprojection-class`] (unused in this default function).
#' @param included_vector A vector of column indices for variables included in the projection (unused in this default function).
#' @param x A matrix of standardized predictors (unused in this default function).
#' @param y A vector of standardized responses (unused in this default function).
#' @param family A [`stats::family`] object (unused in this default function).
#' @param ... Additional arguments (ignored).
#'
#' @return
#' The input `rpm` matrix, unchanged.
#'
#' @keywords internal
update_rpm_identity <- function(rpm, object, included_vector, x, y, family, ...) {
  return(rpm)
}

#' Constructor Function for Building \code{'randomprojection'} Object
#'
#' Creates an object class \code{'randomprojection'} using arguments passed by user.
#' @param name optional string describing the random projection method. This is used for printing.
#' @param generate_fun A function for generating the random projection matrix. This function must accept the following arguments:
#'   - `object`: An object of class [`randomprojection-class`].
#'   - `x`: A matrix of standardized predictors.
#'   - `y`: A vector of standardized responses.
#'   - `m`: The target dimension for the projection.
#'   - `included_vector`: A vector of column indices for the variables to be included in the projection.
#'   - `...`: Additional arguments passed from `spar()` or other functions.' due to the fact that screening is employed pre-projection.
#' @param update_fun A function for updating the [`randomprojection-class`] object with data-specific information.
#'   This function must accept:
#'   - `object`: An object of class [`randomprojection-class`].
#'   - `x`: A matrix of standardized predictors.
#'   - `y`: A vector of standardized responses.
#'   - `family`: A [`stats::family`] object.
#'   - `...`: Additional arguments passed from `spar()`.
#'   If not provided, the default `update_rp_default` is used.
#' @param update_rpm_w_data A function for updating an already-generated random projection matrix with data-dependent information.
#'   This function must accept:
#'   - `rpm`: The random projection matrix to update.
#'   - `object`: An object of class [`randomprojection-class`].
#'   - `included_vector`: A vector of column indices for the variables included in the projection.
#'   - `x`: A matrix of standardized predictors.
#'   - `y`: A vector of standardized responses.
#'   - `family`: A [`stats::family`] object.
#'   - `...`: Additional arguments passed from `spar()`.
#'   If not provided, the default `update_rpm_identity` is used, which leaves the matrix unchanged.
#' @param control A list of control parameters for the random projection. Default: `list()`.
#'   These parameters are passed to `generate_fun`, `update_fun`, and `update_rpm_w_data`.
#' @return Returns a function that, when called, creates and returns an object of class [`randomprojection-class`].
#' @details
#' The `update_rpm_w_data` function is particularly relevant for cross-validation procedures where
#' random projection matrices are precomputed (e.g., `precompute_mode = "precompute_all"` or `precompute_mode = "precompute_rpm"`).
#' In such cases, the cross-validation procedure uses the precomputed matrices, but you may want to update
#' data-dependent entries (e.g., diagonal elements) with the training data in each fold.
#' For example, in `rp_cw(data = TRUE)`, the diagonal elements of the projection matrices are updated
#' to reflect screening coefficients computed on the training data for each fold, while the random elements remain unchanged.
#'
#' @examples
#' generate_cauchy <- function(object, x, y, m, included_vector, ...) {
#'   p <- length(included_vector)
#'   control_rcauchy <- c(object$control[names(object$control) %in% names(formals(rcauchy))],
#'     attributes(object)[names(attributes(object)) %in% names(formals(rcauchy))])
#'   control_rcauchy <-  control_rcauchy[!duplicated(names(control_rcauchy))]
#'   vals <- do.call(function(...)
#'     rcauchy(m * p, ...), control_rcauchy)
#'   RM <- matrix(vals, nrow = m, ncol = p)
#'   return(RM)
#' }
#' rp_cauchy <- constructor_randomprojection(
#'   generate_fun = generate_cauchy, name = "rp_cauchy")
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, rp = rp_cauchy(scale = 1/400))
#' spar_res
#' @export
constructor_randomprojection <- function(name = NULL,
                                         generate_fun,
                                         update_fun = update_rp_default,
                                         update_rpm_w_data = update_rpm_identity,
                                         control = list()) {
  ## Checks
  stopifnot("Function generate_fun needs arguments object,x,y,m,included_vector and ...(ellipsis)."=
              names(formals(generate_fun)) %in% c("object", "m", "included_vector", "x", "y", "..."))
  stopifnot("Function update_fun should have as argument object, x, y, family, ..." = names(formals(update_fun)) %in% c("object","x", "y", "family", "..."))
  stopifnot("Function update_rpm_w_data should have arguments rpm, object, included_vector, x, y, family, ... ."=names(formals(update_rpm_w_data)) %in% c("rpm", "object", "included_vector", "x", "y", "family", "..."))
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
#' @param object object of class  \code{'randomprojection'}.
#' @param x matrix of standardized predictors.
#' @param y vector of standardized response variable.
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return Returns a matrix with \eqn{m} rows and
#'  \code{length(included_vector)} columns sampled from the normal distribution.
#' @keywords internal
generate_gaussian <- function(object,x = NULL, y = NULL, m, included_vector, ...) {
  p <- length(included_vector)
  control_rnorm <- c(
    object$control[names(object$control) %in% names(formals(rnorm))],
    attributes(object)[names(attributes(object)) %in% names(formals(rnorm))])
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
#' @return Returns an object of class \code{'randomprojection'} which is a list with
#' elements \code{name},
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
  name = "rp_gaussian",
  generate_fun = generate_gaussian
)

#'
#' Sparse Random Projection Matrix
#'
#' @param object object of class  \code{'randomprojection'}
#' @param x matrix of standardized predictors
#' @param y vector of standardized response variable
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return Returns a (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_sparse <-  function(object,x = NULL, y = NULL, m, included_vector, ...)  {
  p <- length(included_vector)
  psi <- object$control$psi
  if (is.null(psi)) psi <- attr(object, "psi")
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
#' elements \code{name},
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
  name = "rp_sparse",
  generate_fun = generate_sparse
)


#'
#' Sparse Embedding Matrix
#'
#' @param object object of class  \code{'randomprojection'}
#' @param x matrix of standardized predictors
#' @param y vector of standardized response variable
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_cw <-  function(object, x, y, m, included_vector, ...)  {
  p <- length(included_vector)
  use_data <- attr(object, "data")
  if (is.null(use_data)) {
    diagvals <- sample(c(-1, 1), p, replace = TRUE)
  } else {
    if (use_data) {
      if (is.null(attr(object, "diagvals")))
        stop("Must provide vector of coefficients for data-driven RP.")
      diagvals <- attr(object, "diagvals")[included_vector]
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

update_rp_cw <- function(object, x, y, family, ...) {
  args <- list2(...)
  n <- NROW(x)
  p <- NCOL(x)
  if (attr(args$screencoef, "reuse_in_rp") &&
      !is.null(attr(args$screencoef, "importance"))) {
    scr_coef <- attr(args$screencoef, "importance")
    inc_probs <- attr(args$screencoef, "inc_prob")
    attr(object, "diagvals") <- scr_coef/max(inc_probs)
  } else {
    if (!is.null(object$control$family)) {
      stopifnot("Family provided in control should be of class family." = class(object$control$family) == "family")
      fam <- object$control$family
      object$control$family <- NULL
    } else {
      fam <- family
    }
    object$control$family_string <- paste0(fam$family, "(", fam$link, ")")
    if (fam$family=="gaussian" & fam$link=="identity") {
      fit_family <- "gaussian"
    } else {
      if (fam$family=="binomial" & fam$link=="logit") {
        fit_family <- "binomial"
      } else if (fam$family=="poisson" & fam$link=="log") {
        fit_family <- "poisson"
      } else {
        fit_family <- fam
      }
    }
    if (fam$family=="gaussian") {
      dev.ratio_cutoff <- 0.999
    } else {
      dev.ratio_cutoff <- 0.8
    }

    if (is.null(object$control$alpha)) object$control$alpha <-  0
    if (is.null(object$control$lambda.min.ratio)) {
      object$control$lambda.min.ratio <-
        compute_default_lambda_min_ratio(
          x, y, family = fam)
    }

    control_glmnet <- object$control[names(object$control)  %in% names(formals(glmnet))]
    glmnet_res <- do.call(function(...)
      glmnet(x = x, y = y, family = fit_family, ...), control_glmnet)

    lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
    scr_coef <- coef(glmnet_res,s=lam)[-1]
    inc_probs <- abs(scr_coef)
    max_inc_probs <- max(inc_probs)
    attr(object, "diagvals") <- scr_coef/max_inc_probs
  }
  return(object)
}


update_rpm_w_data_cw <- function(rpm, object, included_vector, x, y, family, ...) {
  rpm@x <-  attr(object, "diagvals")[included_vector]
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
#' elements \code{name},
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
#' Arguments related to the random projection procedure can
#' be passed to the \code{rp_cw()} function through \code{...}, and
#' will be saved as attributes of the \code{'randomprojection'} object.
#' The following attributes are relevant for [spar] and [spar.cv]:
#'  \itemize{
#'  \item \code{mslow}: integer giving the minimum dimension to which the predictors
#'  should be projected; defaults to \eqn{\log(p)}.
#'  \item \code{msup}: integer giving the maximum dimension to which the predictors
#'  should be projected; defaults to \eqn{n/2}.
#'  }
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
  name = "rp_cw",
  generate_fun = generate_cw,
  update_fun = update_rp_cw,
  update_rpm_w_data = update_rpm_w_data_cw
)

#' Print Method for a \code{'randomprojection'} Object
#'
#' @description
#' Prints a summary of a [`randomprojection-class`] object, including its name, and key attributes
#' such as `mslow` and `msup`.
#'
#' @param x An object of class [`randomprojection-class`].
#' @param ... Additional arguments (ignored).
#'
#' @return
#' Invisibly returns the input object `x`.
#'
#' @examples
#' rp <- rp_gaussian(mslow = 5, msup = 10)
#' print(rp)
#'
#' @export
#' @method print randomprojection
print.randomprojection <- function(x, ...) {
  if (!is.null(x$name)) cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  # cat("* Data-dependent:", attr(x,"data"), "\n")
  cat("* Lower bound on goal dimension m:",
      ifelse(is.null(attr(x, "mslow")),
             "not provided, will default to log(p).",
             attr(x, "mslow")), "\n")
  cat("* Upper bound on goal dimension m:",
      ifelse(is.null(attr(x, "msup")),
             "not provided, will default to n/2.",
             attr(x, "msup")), "\n")
}
