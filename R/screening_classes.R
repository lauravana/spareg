#' @keywords internal
update_screen_default <- function(object, x, y, family, ...) {
  object
}


#' Constructor Function for Building \code{'screencoef'} Objects
#'
#' Creates an object class \code{'screencoef'} using arguments passed by user.
#' @param name character
#' @param generate_fun function for generating the screening coefficient. This
#'    function should have arguments  and   \code{y} (vector of responses -- standardized
#'    for Gaussian family), \code{x} (the matrix of standardized predictors) and a
#'    \code{'screencoef'} object.
#' @param update_fun optional function for updating the \code{'screencoef'} object with
#' information from the data passed to `spar()`. This
#' function should have arguments \code{object}, which is a \code{'randomprojection'}
#' object, `x`, y` (the predictor matrix and response vector supplied by `spar()`.
#' Both have already been standardized by `spar()`),
#' `family` and `...`, whereas all other potentially relevant arguments of `spar()`
#' are passed internally to this function through `...`.
#' If `update_fun` is not provided, the object remains unchanged.
#' @return a function which in turn creates an object of class \code{'screencoef'}
#' @description
#' The created function will return a object of class \code{'screencoef'} which
#' constitutes of a list. The attributes of the generating object will include by
#' default \code{type}, which can take one of two values \code{"prob"} (indicating
#' probabilistic screening should be employed),
#' \code{"fixed"} (indicating that the top \code{nscreen} variables should be employed).
#' @examples
#' generate_scr_sirs <- function(y, x, object) {
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
constructor_screencoef <- function(name = NULL, generate_fun,
                                   update_fun = update_screen_default) {
  ## Checks
  args_generate_fun <- formals(generate_fun)
  stopifnot("Function generate_fun should contain three arguments: x, y and an object
            of class \"screencoef\"." =
              length(args_generate_fun) == 3)
  stopifnot("Function generate_fun should contain argument 'y', the vector of responses." =
              "y" %in% names(args_generate_fun))
  stopifnot("Function generate_fun should contain argument 'x', the matrix of predictors." =
              "x" %in% names(args_generate_fun))
  stopifnot("Function update_fun should have as argument object, x, y, family, ..." = names(formals(update_fun)) %in% c("object","x", "y", "family", "..."))
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
generate_scrcoef_marglik <- function(y, x, object) {
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
#' The following attributes are relevant for  [spar] and [spar.cv]:
#' \itemize{
#' \item \code{nscreen} integer giving the number of variables to be retained
#' after screening; if not specified, defaults to $2n$.
#' \item \code{split_data_prop}, double between 0 and 1 which indicates the
#' proportion of the data that should be used for computing the screening
#' coefficient. The remaining data will be used for estimating the marginal
#' models in the SPAR algorithm; if not specified, the whole data will be used
#' for estimating the screening coefficient and the marginal models.
#' \item \code{type} character - either \code{"prob"} (indicating that
#' probabilistic screening should be employed)  or \code{"fixed"} (indicating
#' that a fixed set of \code{nscreen} variables should be employed across the
#' ensemble); defaults to \code{type = "prob"}.
#' \item \code{reuse_in_rp} logical - indicates whether the screening
#' coefficient should be reused at a later stage in the construction of the random
#' projection. Defaults to \code{FALSE}.
#' }
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
generate_scrcoef_cor <- function(y, x, object) {
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
#' The following attributes are relevant for [spar] and [spar.cv]:
#' \itemize{
#' \item \code{nscreen} integer giving the number of variables to be retained
#' after screening; if not specified, defaults to $2n$.
#' \item \code{split_data_prop}, double between 0 and 1 which indicates the
#' proportion of the data that should be used for computing the screening
#' coefficient. The remaining data will be used for estimating the marginal
#' models in the SPAR algorithm; if not specified, the whole data will be used
#' for estimating the screening coefficient and the marginal models.
#' \item \code{type} character - either \code{"prob"} (indicating that
#' probabilistic screening should be employed)  or \code{"fixed"} (indicating
#' that a fixed set of \code{nscreen} variables should be employed across the
#' ensemble); defaults to \code{type = "prob"}.
#' \item \code{reuse_in_rp} logical - indicates whether the screening
#' coefficient should be reused at a later stage in the construction of the random
#' projection. Defaults to \code{FALSE}.
#' }
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
  if (is.null(object$control$family)) {
    object$control$family <-  eval(parse(text=attr(object, "family_string")))
  }
#   if (!is.null(object$control$family)) {
#     fam <- object$control$family
#   } else {
#     fam <- family
#   }
#   if (fam$family == "gaussian" & fam$link=="identity") {
#     fit_family <- "gaussian"
#   } else {
#     if (fam$family=="binomial" & fam$link=="logit") {
#       fit_family <- "binomial"
#     } else if (fam$family=="poisson" & fam$link=="log") {
#       fit_family <- "poisson"
#     } else {
#       fit_family <- fam
#     }
#   }
#   object$control$family <- family# fit_family
  # Set alpha by default to 0
  if (is.null(object$control$alpha)) object$control$alpha <- 0
  # Set default for lambda.min.ratio
  fam <- object$control$family
  if (is.null(object$control$lambda.min.ratio)) {
    n <- NROW(x)
    tmp_sc <- apply(x, 2, function(col) sqrt(var(col)*(n-1)/n))
    # TODO: allow for weights in the future
    x2 <- scale(x, center = colMeans(x), scale = tmp_sc)
    mu0 <- glm(y ~ 1, family = fam)$fitted.values
    r <- y - mu0
    eta <- fam$linkfun(mu0)
    v <- fam$variance(mu0)
    me <- fam$mu.eta(eta)
    rv <- 1/n * r / v * me
    lam_max <- 1000 * max(abs(crossprod(rv, x2[,tmp_sc > 0])))
    object$control$lambda.min.ratio <- min(0.01, 1e-4 / lam_max)
  }

  # Set cutoff ratio for deviance
   if (is.null(object$control$dev.ratio_cutoff)) {
     object$control$dev.ratio_cutoff <- ifelse(fam$family == "gaussian", 0.999, 0.8)
   }
#
   object
#
}

#'
#' Screening coefficient based  on glmnet coefficients
#' @param y vector of responses
#' @param x matrix of predictors
#' @param object  \code{'screencoef'} object
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_glmnet <- function(y, x, object) {
  control_glmnet <-
    object$control[names(object$control) %in% names(formals(glmnet))]

  # Obtain penalized coefs GLMNET
  glmnet_res <- do.call(function(...) glmnet(x = x, y = y, ...),
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
#' The following attributes are relevant for [spar] and [spar.cv]:
#' \itemize{
#' \item \code{nscreen} integer giving the number of variables to be retained
#' after screening; if not specified, defaults to $2n$.
#' \item \code{split_data_prop}, double between 0 and 1 which indicates the
#' proportion of the data that should be used for computing the screening
#' coefficient. The remaining data will be used for estimating the marginal
#' models in the SPAR algorithm; if not specified, the whole data will be used
#' for estimating the screening coefficient and the marginal models.
#' \item \code{type} character - either \code{"prob"} (indicating that
#' probabilistic screening should be employed)  or \code{"fixed"} (indicating
#' that a fixed set of \code{nscreen} variables should be employed across the
#' ensemble); defaults to \code{type = "prob"}.
#' \item \code{reuse_in_rp} logical - indicates whether the screening
#' coefficient should be reused at a later stage in the construction of the random
#' projection. Defaults to \code{FALSE}.
#' }
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
