#' Retrieve the Best Model from a Fitted Object
#'
#' This generic function retrieves the best model from a fitted object of class `spar` or `spar.cv`.
#' The method dispatched depends on the class of the input object.
#'
#' @title get_model
#' @param object An object of class `spar` or `spar.cv`, typically the result of a model fitting function.
#' @param ... Additional arguments passed to the specific method.
#' @return The modified object with the best model selected.
#' @examples
#' # Example usage with a fitted object
#' fitted_obj <- spar(x, y, family = gaussian(), ...)
#' best_model <- get_model(fitted_obj)
#' @export
get_model <- function(object, ...){
  UseMethod("get_model")
}

#' @title Get Best Model for `spar` Objects
#' @description Extracts the best model from a fitted `spar` object based on validation results.
#' The best model is determined by the minimum validation measure.
#'
#' @param object An object of class `spar`, typically the result of `spar()`.
#' @param ... Additional arguments (currently unused).
#' @return An updated `spar` object containing:
#'   - `betas`: Coefficients of the best model, with values below `nu` set to 0.
#'   - `intercepts`: Intercepts for the selected number of models (`nummod`).
#'   - `val_res`: Validation results filtered for the best `nummod` and `nu`.
#' @examples
#' fitted_spar <- spar(x, y, family = gaussian(), ...)
#' best_spar <- get_model(fitted_spar)
#' @export
get_model.spar <- function(object, ...) {
  val_table <- object$val_res
  best_ind <- which.min(val_table$measure)
  parbest <- val_table[best_ind,]
  nummod <- parbest$nummod
  nu <- parbest$nu
  final_coef <- object$betas[, seq_len(nummod), drop=FALSE]
  final_coef[abs(final_coef) < nu] <- 0
  intercepts <- object$intercepts[seq_len(nummod)]

  object$betas <- final_coef
  object$intercepts <- intercepts
  object$val_res <- object$val_res[
    object$val_res$nummod == nummod & object$val_res$nu == nu, ,
    drop = FALSE]

  return(object)
}


#' @title Get Best Model for `spar.cv` Objects
#' @description Extracts the best or 1SE model from a fitted `spar.cv` object.
#' The selection is based on cross-validation results, and the model can be re-estimated if required.
#'
#' @param object An object of class `spar.cv`, typically the result of `spar.cv()`.
#' @param opt_par A character string specifying the selection criterion:
#'   - `"best"`: Selects the model with the minimum validation measure.
#'   - `"1se"`: Selects the simplest model within 1 standard error of the best model.
#' @param x A matrix or data frame of predictors. Required if `fast_fit != "fix_rpm_and_inds"`.
#' @param y A response vector. Required if `fast_fit != "fix_rpm_and_inds"`.
#' @param xval A matrix or data frame of validation predictors. If `NULL`, `x` is used.
#' @param yval A validation response vector. If `NULL`, `y` is used.
#' @param ... Additional arguments passed to `spar()` for re-estimation.
#' @return An updated `spar` or `spar.cv` object containing:
#'   - `betas`: Coefficients of the selected model.
#'   - `intercepts`: Intercepts for the selected model.
#'   - `val_res`: Validation results filtered for the selected `nummod` and `nu`.
#'   - Additional fields like `xscale`, `yscale`, `xcenter`, and `ycenter` if `fast_fit == "fix_rpm_and_inds"`.
#' @examples
#' fitted_spar_cv <- spar.cv(x, y, family = gaussian(), ...)
#' best_spar_cv <- get_model(fitted_spar_cv, opt_par = "best")
#' @export
get_model.spar.cv <- function(object, opt_par = c("best", "1se"),
                              x = NULL, y = NULL, xval = NULL, yval = NULL, ...) {
  opt_nunum <- match.arg(opt_par)
  val_table <- compute_val_summary(object$val_res)
  best_ind <- which.min(val_table$mean_measure)
  parbest <- val_table[best_ind,]
  allowed_ind <- val_table$mean_measure <
    (val_table$mean_measure + val_table$sd_measure)[best_ind]

  ind_1cv <- which.min(val_table$mean_numactive[allowed_ind])
  par1se <- val_table[allowed_ind,][ind_1cv,]
  nummod <- ifelse(opt_nunum == "best", parbest$nummod,
                   par1se$nummod)
  nu <- ifelse(opt_nunum == "best", parbest$nu, par1se$nu)

  if (object$fast_fit == "fix_rpm_and_inds" && is.null(xval) && is.null(yval)) {
    betas_std <- object$fitted_objects[[1]]$betas_std
    final_coef <- betas_std[, seq_len(nummod), drop=FALSE]
    final_coef[abs(final_coef) < nu] <- 0
   # betas <- Matrix(0, ncol(x), nummod, sparse = TRUE)
   # betas[object$fitted_objects[[1]]$xscale > 0,] <- final_coef
    object$betas <- final_coef
    object$xscale <- object$fitted_objects[[1]]$xscale
    object$yscale <- object$fitted_objects[[1]]$yscale
    object$xcenter <- object$fitted_objects[[1]]$xcenter
    object$ycenter <- object$fitted_objects[[1]]$ycenter

    intercepts <- object$fitted_objects[[1]]$intercepts[seq_len(nummod)]
    object$intercepts <- object$fitted_objects[[1]]$intercepts

    val_res <- object$val_res
    val_res <- val_res[val_res$fold == 0 &
      val_res$nummod == nummod & val_res$nu == nu, ,
      drop = FALSE]
    # val_res <- val_res[val_res$fold == 0, !grepl("sd_", colnames(val_res)), drop = FALSE]
    val_res <- val_res[, -c(1:2), drop = FALSE]
    #colnames(val_res) <- c("nu", "nummod", object$measure, "numactive")
    object$val_res <- val_res

    object$fitted_objects <- NULL
    class(object) <- c("spar", "spar.cv")
  } else {
    if (is.null(x) | is.null(y)) {
      stop(sprintf("If fast_fit != 'fix_rpm_and_inds', x and y must be provided for re-estimating the model with the %s (nu, M) combination.",
                   opt_nunum))
    }
    if (is.null(xval) | is.null(yval)) {
      xval <- x; yval <- y
      message("Using x and y as validation data.")
    }
    # Re-estimation of the model with the selected (nu, M) combination
    final_model <- spar(x = x, y = y,
                        family = eval(parse(text = object$family)),
                        model = object$model,
                        rp = object$rp, screencoef = object$screencoef,
                        xval = xval, yval = yval,
                        nus = nu, nummods = nummod,
                       measure = object$measure, avg_type = object$avg_type,
                       seed = object$seed)
    object <- final_model
  }

  return(object)
}



#' Extractor for Model Coefficients from \code{'coefspar'} Object
#' @param x A `\code{coefspar}' object.
#' @return A numeric vector or matrix of coefficients.
#' @seealso [coef.spar], [coef.spar.cv], [print.coefspar], [summary.coefspar]
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' coefs <- coef(spar_res)
#' get_coef(coefs)
#
#' @export
get_coef <- function(x) {
  stopifnot(inherits(x, "coefspar"))
  x$beta
}

#' Extractor for Model Intercept from \code{'coefspar'} Object
#' @param x A `\code{coefspar}' object.
#' @return Intercept (numeric or vector).
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' coefs <- coef(spar_res)
#' get_coef(coefs)
#' @export
get_intercept <- function(x) {
  stopifnot(inherits(x, "coefspar"))
  x$intercept
}


#' Extractor for (Cross-)Validation Measure from '\code{spar}' or '\code{spar.cv}' Object
#'
#' @param object A fitted '\code{spar}' or '\code{spar.cv}'  model
#' @return data.frame containing the (cross-)validation measure for the considered threshold and number of model combinations.
#' For '\code{spar}' objects it contains information about the measure  calculated on the validation set (or on the training sample if
#' xval and yval are missing) and the number of active variables. For '\code{spar.cv}' objects it contains information
#' on the average measure obtained across folds together with the standard deviation across the folds and the average number of active variables.
#' the \code{nfolds} of the training set.
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' get_measure(spar_res)
#'
#' @seealso [spar], [spar.cv], [get_model]
#' @export
get_measure <- function(object) {
  stopifnot(inherits(object, "spar") || inherits(object, "spar.cv"))
  if(inherits(object, "spar")) {
    val_table <- object$val_res
    colnames(val_table)[match("measure", colnames(val_table))] <- object$measure
  } else {
    if(inherits(object, "spar.cv")) {
    val_table <- compute_val_summary(object$val_res)
    colnames(val_table)[4] <- paste0("mean_", object$measure)
    colnames(val_table)[5] <- paste0("sd_", object$measure)
    colnames(val_table)[6] <- "mean_numactive"
    }
  }

  val_table[, !(colnames(val_table) %in% c("nnu"))]
}


