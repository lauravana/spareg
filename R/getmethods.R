#' get_model
#'
#' A details of get_model
#'
#' @title get_model: The get_model function
#' @param object numeric number
#' @param ... other arguments
#' @examples
#' get_model(a)
#' @export
get_model <- function(object, ...){
  UseMethod("get_model")
}

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


