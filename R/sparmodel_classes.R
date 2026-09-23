#' @keywords internal
update_sparmodel_default <- function(object, x, y, family, ...) {
  if (is.null(object$control$family)) object$control$family <- family
  object
}

#' Constructor Function for Building \code{'sparmodel'} Object
#'
#' Creates an object of class \code{'sparmodel'} using arguments passed by user.
#' @param generate_fun function for estimating the marginal models which returns the
#     intercept and the vector of coefficients. This
#'    function should have arguments  a
#'    \code{'sparmodel'} \code{object}, \code{x} and \code{y} (matrix of predictors and
#'    vector of responses supplied by `spar()`, which have already been standardized),
#'    \code{z} (the matrix of projected predictors) and `...`.
#' @param name optional string describing the model employed. This is used for printing.
#' @param update_fun optional function for updating the \code{'sparmodel'} object
#'  before the
#' start of the algorithm.
#' @return a function which in turn creates an
#'    object of class \code{'sparmodel'}.
#' @description
#' The created function will return a object of class \code{'sparmodel'} which
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
#' Penalized  GLM Marginal  \code{'sparmodel'}
#'
#' @description
#' Creates an object class \code{'sparmodel'} using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' \code{'sparmodel'} object
#' @param control list of controls to be passed to the model function
#' @return object of class \code{'sparmodel'} which is a list with elements
#' \itemize{
#'  \item \code{name} (character, optional, used for printing)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  function for estimating the marginal models which returns the
#     intercept and the vector of coefficients. This
#'    function should have arguments  a
#'    \code{'sparmodel'} \code{object}, \code{x} and \code{y} (matrix of predictors and
#'    vector of responses supplied by `spar()`, which have already been standardized),
#'    \code{z} (the matrix of projected predictors) and `...`.
#'   Returns a list with
#'   two elements: \code{gammas} which is the vector of regression coefficients
#'    for the projected predictors and \code{intercept} which is the intercept
#'    of the model.
#'  \item \code{update_fun}  optional function for updating the \code{'sparmodel'}
#'   object before the start of the algorithm.This
#'   function should have arguments \code{object}, which is a \code{'sparmodel'}
#'   object, `x` (the matrix of standardized predictors), `y` (the vector of standardized responses),
#'   `family` and `...`, whereas all other potentially relevant arguments of `spar()`
#'    are passed internally to this function through `...`. For
#'    \code{spar_glmnet()} this function manipulates the
#'     \code{'family'} object for compatibility with
#'    \link[glmnet]{glmnet}. I.e., In the case of families Gaussian,
#'     binomial and Poisson with  canonical link, the family object is
#'     replaced by a string containing the name of the family.
#'     This leads to  \link[glmnet]{glmnet} using the faster
#'     specialized algorithms rather than the general algorithm.
#' }
#' @details
#' Relies on \link[glmnet]{glmnet}.
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

#' GLM Marginal \code{'sparmodel'}
#'
#' @description
#' Creates an object class \code{'sparmodel'} using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' \code{'sparmodel'} object
#' @param control list of controls to be passed to the model function
#' @return object of class \code{'sparmodel'} which is a list with elements
#' \itemize{
#'  \item \code{name} (character, optional, used for printing)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun} function for estimating the model coefficients and the intercept.
#'  This
#'    function should have arguments  a
#'    \code{'sparmodel'} \code{object}, \code{x} and \code{y} (matrix of predictors and
#'    vector of responses supplied by `spar()`, which have already been standardized),
#'    \code{z} (the matrix of projected predictors) and `...`. Returns a list with
#'    two elements: \code{gammas} which is the vector of regression coefficients
#'    for the projected predictors and \code{intercept} which is the intercept of the model
#' }
#' @details
#' Relies on \link[stats]{glm}.
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

