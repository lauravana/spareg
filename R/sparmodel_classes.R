#' SPAR Model Object Class
#'
#' @description
#' The `sparmodel` class represents a configuration for **fitting marginal models** in the
#' **SPAR (Sparse Projected Averaged Regression)** framework. Objects of this class encapsulate:
#'   - Functions for estimating marginal models (e.g., GLM, GLMNET, or custom models).
#'   - Control parameters for the model fitting process.
#'   - Metadata (e.g., name, attributes) for customization.
#'
#' Marginal models are fitted to **projected predictors** (after screening and random projection) to compute
#' coefficients and intercepts for each model in the ensemble.
#'
#' @details
#' Objects of class `sparmodel` are created using the `constructor_sparmodel` function or predefined constructors
#' like `spar_glm` and `spar_glmnet`. They are used in the `spar` and `spar.cv` functions to define how marginal models
#' are fitted to the projected data.
#'
#' The class includes the following components:
#'   - **`name`**: A character string describing the model type (e.g., `"glm"`, `"glmnet"`, `"glmrob"`).
#'   - **`generate_fun`**: A function to estimate the marginal model coefficients and intercept. This function must accept:
#'     - `object`: An object of class `sparmodel`.
#'     - `x`: A matrix of standardized predictors (unused in most cases, as models are fitted on projected data).
#'     - `y`: A vector of standardized responses.
#'     - `z`: A matrix of projected predictors (the design matrix for the marginal model).
#'     - `...`: Additional arguments passed from `spar()`.
#'     The function must return a list with two elements:
#'       - `gammas`: A vector of regression coefficients for the projected predictors.
#'       - `intercept`: The intercept of the model.
#'   - **`update_fun`**: A function to update the `sparmodel` object with data-specific information before fitting.
#'     This function must accept:
#'     - `object`: An object of class `sparmodel`.
#'     - `x`: A matrix of standardized predictors.
#'     - `y`: A vector of standardized responses.
#'     - `family`: A [`stats::family`] object.
#'     - `...`: Additional arguments passed from `spar()`.
#'     If not provided, the default `update_sparmodel_default` is used, which only updates the `family` if missing.
#'   - **`control`**: A list of control parameters for the model fitting process (e.g., `family`, `alpha` for `glmnet`).
#'

#' @section Model Types:
#' The following predefined model types are available:
#'   - **`spar_glm`**: Fits marginal models using [`stats::glm`]. Supports all standard GLM families.
#'     For Gaussian models with identity link, it uses a fast OLS solver.
#'   - **`spar_glmnet`**: Fits marginal models using [`glmnet::glmnet`]. Supports penalized regression (Ridge, Lasso, or Elastic Net).
#'     Automatically converts `family` objects to strings for compatibility with `glmnet`.
#'   - **Custom Models**: Users can define their own model types by providing custom `generate_fun` and `update_fun`
#'     functions to `constructor_sparmodel`. For example, `spar_glmrob` (not exported by default) uses
#'     [`robustbase::glmrob`] for robust regression.
#'
#' @examples
#' model_glmnet <- spar_glmnet(alpha = 0.5)
#'
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(
#'   example_data$x,
#'   example_data$y,
#'   xval = example_data$xtest,
#'   yval = example_data$ytest,
#'   model = model_glmnet
#' )
#'
#' @seealso
#' [`constructor_sparmodel`], [`spar_glm`], [`spar_glmnet`], [`spar`]
#'
#'
#' @name sparmodel-class
NULL


#' @keywords internal
update_sparmodel_default <- function(object, x, y, family, ...) {
  if (is.null(object$control$family)) object$control$family <- family
  object
}

#' Constructor Function for Building [`sparmodel-class`] Object
#'
#' Creates an object of class [`sparmodel-class`] using arguments passed by user.
#' @param generate_fun function for estimating the marginal models which returns the
#     intercept and the vector of coefficients. This
#'    function should have arguments  a
#'    [`sparmodel-class`] \code{object}, \code{x} and \code{y} (matrix of predictors and
#'    vector of responses supplied by `spar()`, which have already been standardized),
#'    \code{z} (the matrix of projected predictors) and `...`,
#'    whereas all other potentially relevant arguments of `spar()`
#'   are passed internally to this function through the ellipsis.
#'
#' @param name optional string describing the model employed. This is used for printing.
#' @param update_fun optional function for updating the [`sparmodel-class`] object. This
#' function should have arguments \code{object}, which is a \code{'screencoef'}
#' object, `x`, `y` (the predictor matrix and response vector supplied by `spar()`
#' where both have already been standardized by `spar()`),
#' `family` and `...`, whereas all other potentially relevant arguments of `spar()`
#' are passed internally to this function through the ellipsis.
#' @return Returns a function that, when called, creates and returns an object of class [`sparmodel-class`].
#' @description
#' The created function will return a object of class [`sparmodel-class`] which
#' constitutes of a list.
#' @examples
#' model_glmrob <- function(object, x, y, z, ...) {
#'   requireNamespace("robustbase")
#'   glmrob_res <- do.call(function(...)
#'     robustbase::glmrob(y ~ as.matrix(z), ...),
#'     object$control)
#'   intercept <- coef(glmrob_res)[1]
#'   gammas <- coef(glmrob_res)[-1]
#'   list(gammas = gammas, intercept = intercept)
#' }
#' spar_glmrob <- constructor_sparmodel(
#'   generate_fun = model_glmrob,name = "glmrob")
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest,
#'   model = spar_glmrob())
#' spar_res
#' @export
constructor_sparmodel <- function(name = NULL, generate_fun,
                                  update_fun = update_sparmodel_default) {
  ## Checks ----
  args_generate_fun <- formals(generate_fun)
  stopifnot("Function generate_fun should contain five arguments: object, x, y, z and ... (ellipsis)." =
              length(args_generate_fun) == 5)
  stopifnot("Function generate_fun should contain argument 'x', the matrix of responses." =
              "x" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'y', the vector of responses." =
              "y" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'z', the matrix of reduced predictors." =
              "z" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'object', an object
            of class \"sparmodel\"" =
              "object" %in% names(args_generate_fun))
  stopifnot("Function update_fun should have as argument object, x, y, family, ..." = names(formals(update_fun)) %in% c("object","x", "y", "family", "..."))

  ## Function to return  ----
  function(..., control = list()) {
    out <- list(name = name,
                generate_fun = generate_fun,
                update_fun = update_fun,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    class(out) <- c("sparmodel")
    return(out)
  }
}
#'
#' Penalized  GLM Marginal  [`sparmodel-class`]
#'
#' @description
#' Creates an object class [`sparmodel-class`] using arguments passed by user, where
#' the generating function computes the coefficients of the marginal models
#' based on a penalized GLM. Computation relies on [`glmnet::glmnet`].
#' By default, the models assume \eqn{\alpha=0} and return the coefficients
#' obtained with the penalty \eqn{\lambda_\text{min}}.
#'
#' @param ... includes arguments which can be passed as attributes to the
#' [`sparmodel-class`] object
#' @param control list of controls to be passed to the model function
#' @return Returns an object of class [`sparmodel-class`].
#' @details
#' Relies on \link[glmnet]{glmnet}.
#' @seealso [`sparmodel-class`], [`glmnet::glmnet`]
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y,
#'   xval = example_data$xtest, yval = example_data$ytest,
#'   model = spar_glmnet(alpha = 0.1))
#'
#' @export
#'
spar_glmnet <- function(..., control = list()) {
  out <-  list(name = "glmnet",
               generate_fun = model_glmnet,
               update_fun = update_sparmodel_glmnet,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

ols_fun <- function(y, z) {
  solve(crossprod(z), crossprod(z,y))
}

ols_fun_corrected <- function(y, z) {
  ZtZ <- crossprod(z) + 0.01 * diag(ncol = ncol(z), nrow = ncol(z))
  solve(ZtZ, crossprod(z,y))
}



update_sparmodel_glmnet <- function(object, z, yz, family, ...) {
  # Family compatibility with glmnet
  if (!is.null(object$control$family)) {
    fam <- object$control$family
  } else {
    fam <- family
  }
  if (fam$family == "gaussian" & fam$link=="identity") {
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
  object$control$family <- fit_family
  ## Set defaults
  if (is.null(object$control$alpha)) {
    object$control$alpha <- 0
  }
  ## Return
  object
}

model_glmnet <- function(object, x = NULL, y, z, ...) {
  ## y - vector of n responses
  ## z - matrix of reduced predictors with n rows
  glmnet_res <- do.call(function(...) glmnet(x = z, y = y, ...),
                        object$control)
  mar_coef <- coef(glmnet_res, s = min(glmnet_res$lambda))
  intercept <- mar_coef[1]
  gammas <- mar_coef[-1]
  list(gammas = gammas, intercept = intercept)
}

#' GLM Marginal [`sparmodel-class`]
#'
#' @description
#' Creates an object class [`sparmodel-class`] using arguments passed by user.
#' The generating function computes the coefficients of the marginal models
#' based on a GLM. Computation relies on [`stats::glm`].
#'
#' @param ... includes arguments which can be passed as attributes to the
#' [`sparmodel-class`] object
#' @param control list of controls to be passed to the model function
#' @return Returns an object of class [`sparmodel-class`].
#' @details
#' Relies on \link[stats]{glm}.
#' @seealso [`sparmodel-class`], [`stats::glm`]
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y,
#'   xval = example_data$xtest, yval = example_data$ytest,
#'   model = spar_glm())
#' @export
#'
spar_glm <- function(..., control = list()) {
  out <-  list(name = "glm",
               generate_fun = model_glm,
               update_fun = update_sparmodel_default,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

model_glm <- function(object, x = NULL, y, z, ...) {
  ## y - vector of n responses
  ## z - matrix with n rows
  family <- object$control$family
  if (family$family=="gaussian" & family$link=="identity") {
    intercept <- 0
    gammas <- tryCatch(ols_fun(y, z),
                       error = function(error_message) {
                         return(ols_fun_corrected(y, z))
                       })
  } else {
    glm_res <- do.call(function(...) glm.fit(y = y, x = cbind(1, z), ...),
                       object$control)
    intercept <- glm_res$coefficients[1]
    gammas <- glm_res$coefficients[-1]
  }

  list(gammas = gammas, intercept = intercept)
}

