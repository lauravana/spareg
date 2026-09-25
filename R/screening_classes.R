#' Screening Coefficient Object Class
#'
#' @description
#' The `'screencoef'` class represents a configuration for computing and managing **screening coefficients**
#' in the **SPAR (Sparse Projected Averaged Regression)** framework. Objects of this class encapsulate:
#'   - Functions for generating and updating screening coefficients.
#'   - Control parameters for the screening process.
#'   - Metadata (e.g., name, attributes) for customization.
#'
#' Screening coefficients are used to reduce the dimensionality of the predictor space by selecting
#' the most relevant variables before applying random projections.
#'
#' @details
#' Objects of class `screencoef` are created using the `constructor_screencoef` function. They are used in the
#' `spar` and `spar.cv` functions to define how predictors are screened.
#'
#' The class includes the following components:
#'   - **`name`**: A character string describing the screening method (e.g., `"screen_marglik"`, `"screen_cor"`, `"screen_glmnet"`).
#'   - **`generate_fun`**: A function to compute the screening coefficients. This function must accept:
#'     - `x`: A matrix of standardized predictors.
#'     - `y`: A vector of standardized responses.
#'     - `object`: An object of class `screencoef`.
#'     - `...`: Additional arguments passed from `spar()`.
#'   - **`update_fun`**: A function to update the `screencoef` object with data-specific information. This function must accept:
#'     - `object`: An object of class `screencoef`.
#'     - `x`: A matrix of standardized predictors.
#'     - `y`: A vector of standardized responses.
#'     - `family`: A [`stats::family`] object.
#'     - `...`: Additional arguments passed from `spar()`.
#'     If not provided, the default `update_screen_default` is used, which leaves the object unchanged.
#'   - **`control`**: A list of control parameters for the screening process (e.g., `nscreen`, `split_data_prop`).
#'
#' @section Attributes:
#' The following attributes are commonly used in `screencoef` objects:
#'   - **`type`**: Character. The type of screening to employ:
#'     - `"prob"`: Probabilistic screening (variables are selected probabilistically based on their screening coefficients).
#'     - `"fixed"`: Fixed screening (the top `nscreen` variables are selected).
#'     Default: `"prob"`.
#'   - **`nscreen`**: Integer. The number of variables to retain after screening. Default: `2n` (twice the number of observations).
#'   - **`split_data_prop`**: Numeric. The proportion of data to use for computing screening coefficients.
#'     The remaining data is used for fitting the marginal models. Default: `1` (use all data).
#'   - **`reuse_in_rp`**: Logical. If `TRUE`, the screening coefficients are reused in the construction of the random projection.
#'     Default: `FALSE`.
#'   - **`importance`**: Numeric vector. The screening coefficients computed from the data.
#'   - **`inc_prob`**: Numeric vector. The inclusion probabilities for probabilistic screening (normalized screening coefficients).
#'
#' @section Screening Methods:
#' The following predefined screening methods are available:
#'   - **`screen_marglik`**: Screening based on marginal likelihood in univariate GLMs.
#'     Uses [`stats::glm`] to fit a separate GLM for each predictor and extracts the coefficients.
#'   - **`screen_cor`**: Screening based on correlation between predictors and the response.
#'     Uses [`stats::cor`] to compute correlation coefficients.
#'   - **`screen_glmnet`**: Screening based on penalized regression coefficients from [`glmnet::glmnet`].
#'     Uses Lasso or Ridge regression to compute coefficients.
#'
#' Users can also define custom screening methods by providing their own `generate_fun` and `update_fun` functions
#' to `constructor_screencoef`.
#'
#' @examples
#' screen_marglik_obj <- screen_marglik(nscreen = 500, type = "prob")
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(
#'   example_data$x,
#'   example_data$y,
#'   xval = example_data$xtest,
#'   yval = example_data$ytest,
#'   screencoef = screen_marglik_obj
#' )
#'
#' @seealso
#' [`constructor_screencoef`], [`screen_marglik`], [`screen_cor`], [`screen_glmnet`], [`spar`]
#'
#' @name screencoef-class
NULL

#' @keywords internal
update_screen_default <- function(object, x, y, family, ...) {
  object
}


#' Constructor Function for Building [`screencoef-class`] Objects
#'
#' Creates an object class [`screencoef-class`]  using arguments passed by user.
#' @param name character
#' @param generate_fun function responsible for computing the screening coefficients. This
#'    function should have arguments
#'    `x`, `y` (the predictor matrix and response vector supplied by `spar()`
#'      where both have already been standardized by `spar()`),
#'    \code{'screencoef'} object and `...`, whereas all other potentially
#'    relevant arguments of `spar()` are passed internally to this function through the ellipsis.
#' @param update_fun optional function for updating the [`screencoef-class`]  object with
#' information from the data passed to `spar()`. This
#' function should have arguments \code{object}, which is a [`screencoef-class`]
#' object, `x`, `y` (the predictor matrix and response vector supplied by `spar()`
#' where both have already been standardized by `spar()`),
#' `family` and `...`, whereas all other potentially relevant arguments of `spar()`
#' are passed internally to this function through the ellipsis.
#' If `update_fun` is not provided, the object remains unchanged.
#' @return Returns a function that, when called, creates and returns an object of class `'screencoef'`.
#' @description
#' The created function will return a object of class [`screencoef-class`]  which
#' constitutes of a list. The attributes of the generating object will include by
#' default \code{type}, which can take one of two values \code{"prob"} (indicating
#' probabilistic screening should be employed),
#' \code{"fixed"} (indicating that the top \code{nscreen} variables should be employed).
#' @examples
#' generate_scr_sirs <- function(y, x, object, ...) {
#'   ctrl <- object$control[names(object$control) %in%
#'                           names(formals(VariableScreening::screenIID))]
#'  res_screen <- do.call(function(...)
#'    VariableScreening::screenIID(x, y, ...), ctrl)
#'  res_screen$measurement
#'}
#' screen_sirs <- constructor_screencoef("screen_sirs",
#'   generate_fun = generate_scr_sirs)
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_example <- spar(example_data$x, example_data$y,
#'   screencoef = screen_sirs(control = list(method = "SIRS")),
#'   rp = rp_sparse())
#' spar_example
#' @export
constructor_screencoef <- function(name = NULL,
                                   generate_fun,
                                   update_fun = update_screen_default) {
  ## Checks
  args_generate_fun <- formals(generate_fun)
  stopifnot("Function generate_fun should contain four arguments: x, y, an object
            of class \"screencoef\" and ... (ellipsis)." =
              length(args_generate_fun) == 4)
  stopifnot("Function generate_fun should contain argument 'y', the vector of standardized responses." =
              "y" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'x', the matrix of standardized predictors." =
              "x" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'object', an object
            of class \"screencoef\"" = "object" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument '...', which allows to pass all other arguments of spar() internally to this function."
            = "..." %in% names(args_generate_fun))
  stopifnot("Function update_fun should have as argument object, x, y, family, ..." =
              names(formals(update_fun)) %in% c("object","x", "y", "family", "..."))
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                generate_fun = generate_fun,
                update_fun = update_fun,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    if (is.null(attr(out, "type"))) {
      attr(out, "type") <- "prob"
    } else {
      stopifnot(
        "'type' must be either 'prob' or 'fixed'." =
          (attr(out, "type") == "prob" | attr(out, "type") == "fixed")
      )
    }
    if (is.null(attr(out, "reuse_in_rp"))) attr(out, "reuse_in_rp") <- FALSE
    class(out) <- c("screencoef")
    return(out)
  }
}

#' @keywords internal
update_screen_marglik <- function(object, x, y, family, ...) {
  if (is.null(object$control$family)) object$control$family <- family
  object
}
#'
#' Generate screening coefficient based  on marginal likelihood in univariate GLMs
#' @param y vector of responses
#' @param x matrix of predictors
#' @param object  \code{'screencoef'} object
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_marglik <- function(object, x, y, ...) {
  ctrl <- object$control[names(object$control) %in%
                           names(formals(glm))]
  coefs <- apply(x, 2, function(xj){
    glm_res <- do.call(function(...) glm(y ~ xj,  ...),
                       ctrl)
    glm_res$coefficients[2]
  })
  coefs
}

#' Screening Coefficient Based on Marginal GLMs
#'
#' @param ... includes arguments which can be passed as attributes to the
#' \code{'screencoef'} object
#' @param control list of controls to be passed to the screening function
#' @return object of class \code{'screencoef'} which is a list with elements:
#'
#' \itemize{
#'  \item \code{name} (character, optional, used for printing)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient.
#'  This function should have arguments  and   \code{y} (vector of (standardized for Gaussian) responses),
#'  \code{x} (the matrix of standardized predictors) and a \code{'screencoef'} object.
#'
#' }
#' @description
#' Creates an object class \code{'screencoef'} using arguments passed by user,
#' where the screening coefficient should be computed based on the marginal
#' likelihood of the univariate GLM where the response is regressed on
#' each predictor separately.
#'
#' @details
#' The function \code{generate_fun} relies on \link[stats]{glm}.
#'
#' Arguments related to the screening procedure can
#' be passed to the \code{screen_marglik()} function through \code{...}, and
#' will be saved as attributes of the \code{'screencoef'} object.
#' Note that if \code{family} is not provided in \code{control},
#' the \code{family} used in [spar]  or [spar.cv]  will be used.
#'
#' @seealso
#' [`constructor_screencoef`], [`screencoef-class`]
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   screencoef = screen_marglik(nscreen = 500))
#'
#' @export
#'
screen_marglik <- constructor_screencoef(
  name = "screen_marglik",
  generate_fun = generate_scrcoef_marglik,
  update_fun = update_screen_marglik)


#'
#' Generate screening coefficient based  on correlation
#'
#' @param y vector of responses
#' @param x matrix of predictors
#' @param object  \code{'screencoef'} object
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_cor <- function(object, x, y, ...) {
  coefs <- apply(x, 2, function(xj) {
    do.call(function(...) cor(y, xj, ...),
            object$control)
  })
  coefs
}



#' Screening Coefficient Based on Correlation
#'
#' Creates an object class \code{'screencoef'} using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' \code{'screencoef'} object
#' @param control list of controls to be passed to the screening function
#' @return object of class \code{'screencoef'} which is a list with elements
#' \itemize{
#'  \item \code{name} (character, optional, used for printing)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient.
#'  This function should have arguments  and   \code{y} (vector of (standardized for Gaussian) responses),
#'  \code{x} (the matrix of standardized predictors) and a \code{'screencoef'} object.
#' }
#'
#' @description
#' Creates an object class \code{'screencoef'} using arguments passed by user,
#' where the screening coefficient should be computed based on the correlation
#' coefficient of response and each predictor separately.
#'
#' @details
#' The function \code{generate_fun} relies on \link[stats]{cor}.
#'
#' Arguments related to the screening procedure can
#' be passed to the \code{screen_cor()} function through \code{...}, and
#' will be saved as attributes of the \code{'screencoef'} object.
#'
#' @seealso
#' [`constructor_screencoef`], [`screencoef-class`]
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   screencoef = screen_cor(control = list(method = "kendall")))
#'
#' @export
#'
screen_cor <- constructor_screencoef(
  name = "screen_cor",
  generate_fun = generate_scrcoef_cor)


#'
#' Update screening glmnet object
#' @param object  \code{'screencoef'} object
#' @param y vector of responses
#' @param x matrix of predictors
#' @param family family object to be passed from \code{spar()}
#' @return  vector of screening coefficients of length p
#' @keywords internal
update_screen_glmnet <- function(object, x, y, family, ...) {
  # Family compatibility with glmnet
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
  object$control$fit_family <- fit_family
  # Set alpha by default to 0
  if (is.null(object$control$alpha)) object$control$alpha <- 0
  # Set default for lambda.min.ratio
  if (is.null(object$control$lambda.min.ratio)) {
    object$control$lambda.min.ratio <-
      compute_default_lambda_min_ratio(
        x, y, family = fam)
  }

  # Set cutoff ratio for deviance
  if (is.null(object$control$dev.ratio_cutoff)) {
    object$control$dev.ratio_cutoff <- ifelse(fam$family == "gaussian", 0.999, 0.8)
  }

  object
}

#'
#' Screening coefficient based  on glmnet coefficients
#' @param y vector of responses
#' @param x matrix of predictors
#' @param object  \code{'screencoef'} object
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_glmnet <- function(object, x, y, ...) {
  control_glmnet <-
    object$control[names(object$control) %in% names(formals(glmnet))]

  # Obtain penalized coefs GLMNET
  glmnet_res <- do.call(function(...) glmnet(x = x, y = y,
                                             family = object$control$fit_family,
                                             ...),
                        control_glmnet)

  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= object$control$dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res, s = lam)[-1]
  scr_coef
}

#'
#' Screening Coefficient Based on \link[glmnet]{glmnet} Coefficients
#'
#' Creates an object class \code{'screencoef'} using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' \code{'screencoef'} object
#' @param control list of controls to be passed to the screening function
#' @return object of class \code{'screencoef'} which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character, optional, used for printing)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient.
#'  This function should have arguments  and   \code{y} (vector of (standardized for Gaussian) responses),
#'  \code{x} (the matrix of standardized predictors) and a \code{'screencoef'} object.
#' }
#'
#' @description
#' Creates an object class \code{'screencoef'} using arguments passed by user,
#' where the screening coefficient should be computed based on penalized coefficients.
#'
#' @details
#' The function \code{generate_fun} relies on \link[glmnet]{glmnet}.
#'
#' Arguments related to the screening procedure can
#' be passed to the \code{screen_glmnet()} function through \code{...}, and
#' will be saved as attributes of the \code{'screencoef'} object.
#'
#' @seealso
#' [`constructor_screencoef`], [`screencoef-class`]
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30),
#'   screencoef = screen_glmnet(control = list(alpha = 0.1)))
#'
#' @export
#'
screen_glmnet <- constructor_screencoef(
  name = "screen_glmnet",
  generate_fun = generate_scrcoef_glmnet,
  update_fun = update_screen_glmnet)


#' Print Method for \code{'screencoef'} Object
#'
#' Print method for a \code{'screencoef'} object
#' @param x description
#' @param ... further arguments passed to or from other methods
#' @return text summary
#'
#' @export
print.screencoef <- function(x, ...) {
  if (!is.null(x$name)) cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  cat("* proportion of data used for screening:",
      ifelse(is.null(attr(x, "split_data_prop")),
             1, attr(x, "split_data_prop")), "\n")
  cat("* number of screened variables:",
      ifelse(is.null(attr(x, "nscreen")),
             "not provided, will default to 2n",
             attr(x, "nscreen")), "\n")
  cat("* type:",  ifelse(attr(x, "type") == "prob",
                         "probabilistic screening",
                         "screening top nscreen variables"), "\n")
  imp_vals <- attr(x, "importance")
  out_imp <-  ifelse(!is.null(imp_vals),
                     sprintf("num [1:%d] %s ...", length(imp_vals),
                             paste(round(imp_vals[1:5], 3),
                                   collapse = " ")),
                     "not (yet) computed from the data.")
  cat("* screening coefficients:", out_imp,  "\n")
}
