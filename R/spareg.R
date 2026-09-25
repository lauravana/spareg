###########################################
### Main implementation of sparse projected averaged regression (SPAR)
##########################################
#' Sparse Projected Averaged Regression (SPAR)
#'
#' Fits a **Sparse Projected Averaged Regression (SPAR)** model to high-dimensional data.
#' SPAR builds an ensemble of generalized linear models (GLMs) where high-dimensional predictors
#' are first screened (using a screening coefficient) and then projected using random projection matrices.
#' This function evaluates the model over a grid of thresholds (`nus`) and a grid of the number of marginal models (`nummods`).
#' It is also used internally by the cross-validated procedure [`spar.cv`].
#'
#' @param x An `n x p` numeric matrix of predictor variables.
#' @param y A quantitative response vector of length `n`.
#' @param family  A \link[stats]{family}  object used for the marginal generalized linear model,
#'        default \code{gaussian("identity")}.
#' @param model A function that creates a \code{'sparmodel'} object. Defaults to:
#'   - \code{spar_glm()} for gaussian family with identity link.
#'   - \code{spar_glmnet()} for all other family-link combinations.
#' @param rp A function that creates a \code{'randomprojection'} object. Defaults to `NULL`.
#' If `NULL`, \code{rp_cw(data = TRUE)} is used.
#' @param screencoef function that creates a \code{'screencoef'} object. Defaults to `NULL`.
#' If `NULL`, no screening is used.
#' @param xval An optional matrix of predictor variables used for
#'        validation. If `NULL`, \code{x} is used.
#' @param yval An optional response vector for validation.
#' If `NULL`, \code{y} is used.
#' @param nnu number of different threshold values \eqn{\nu} for thresholding.
#'        Ignored if `nus` is provided. Defaults to `20L`.
#' @param nus An optional vector of thresholds \eqn{\nu}.
#'         If `NULL`, `nnu` values between `0` and the maximum absolute marginal
#'         coefficient equally spaced  on the probability scale are used (i.e.,
#'         we use as grid for the thresholds
#'         zero and \code{nnu-1} quantiles of the
#'         absolute values of the estimated non-zero coefficients
#'         from the marginal models).
#' @param nummods A vector of integers specifying the number of marginal models to consider for validation.
#'         Defaults to \code{20L}.
#' @param measure The loss function for validation. Options:
#'  - \code{"deviance"} (default, available for all families).
#'  - \code{"mse"} or \code{"mae"} (mean squared/absolute error, for all families).
#'  - \code{"class"} (misclassification error, for binomial family only).
#'  - \code{"1-auc"} (one minus area under the ROC curve for binomial family only).
#' @param avg_type The type of averaging for marginal models. Options:
#'   - `"link"` (default): Averaging on the link scale.
#'   - `"response"`: Averaging on the response scale.
#' @param parallel A logical indicating whether to use parallel
#'        estimation of the marginal models. Defaults to `FALSE`.
#' @param inds An optional list of index vectors corresponding to variables retained after screening for each marginal model.
#'   Must have length `max(nummods)`.
#' @param RPMs An optional list of projection matrices for each marginal model.
#'   Must have length `max(nummods)`.
#' @param seed An optional integer seed for reproducibility. Default: `NULL`.
#' @param ... Additional arguments for backward compatibility.
#'
#' @returns
#' An object of class [`spar`] with the following components:
#'   - **`betas`**: A `p x max(nummods)` sparse matrix of standardized coefficients for each marginal model.
#'   - **`intercepts`**: Intercepts for each marginal model.
#'   - **`scr_coef`**: A vector of length `p` with screening coefficients for standardized predictors.
#'   - **`inds`**: A list of index vectors for variables retained after screening.
#'   - **`RPMs`**: A list of projection matrices for each marginal model.
#'   - **`val_res`**: A `data.frame` with validation results (measure and number of active variables) for each \eqn{M} and \eqn{\nu}.
#'   - **`val_set`**: A logical flag indicating whether validation data were provided.
#'   - **`family`**: The family object used for the GLM.
#'   - **`nus`**: The vector of thresholds considered.
#'   - **`nummods`**: The vector of numbers of marginal models considered.
#'   - **`ycenter`**: The empirical mean of the initial response vector.
#'   - **`yscale`**: The empirical standard deviation of the initial response vector.
#'   - **`xcenter`**: A vector of empirical means for the initial predictors.
#'   - **`xscale`**: A vector of empirical standard deviations for the initial predictors.
#'   - **`avg_type`**: The averaging type used for validation.
#'   - **`measure`**: The validation measure used.
#'   - **`rp`**: The `'randomprojection'` object.
#'   - **`screencoef`**: The `'screencoef'` object.
#'   - **`x_rows_for_fitting_marginal_models`**: A vector of row indices from `x` used for fitting marginal models (if screening splits data).

#' @details
#' If a parallel backend (e.g., `doParallel`) is registered and `parallel = TRUE`,
#' the [`foreach`] package is used to parallelize the estimation of marginal models.
#' If a parallel backend is registered and \code{parallel = TRUE},
#' the \link[foreach]{foreach} function
#' is used to estimate the marginal models in parallel.
#'
#' @references{
#'   \insertRef{parzer2024lm}{spareg}
#'
#'   \insertRef{parzer2024glms}{spareg}
#'
#'   \insertRef{Clarkson2013LowRankApprox}{spareg}
#'
#'   \insertRef{ACHLIOPTAS2003JL}{spareg}
#' }
#' @examples
#' example_data <- simulate_spareg_data(n = 200, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30))
#' coefs <- coef(spar_res)
#' pred <- predict(spar_res, xnew = example_data$x)
#' plot(spar_res)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nu", nummod = 10)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nu", nummod = 10)
#' plot(spar_res, plot_type = "res_vs_fitted",  xfit = example_data$xtest,
#'   yfit = example_data$ytest)
#' plot(spar_res, plot_type = "coefs", prange = c(1,400))
#'
#' @seealso [spar.cv], [coef.spar], [predict.spar], [plot.spar], [print.spar]
#' @aliases spareg
#' @export
#'
#' @import methods
#' @importFrom stats median reshape glm.fit coef fitted gaussian predict rnorm quantile
#'  residuals sd var cor glm aggregate
#' @importFrom utils head
#' @importFrom Matrix Matrix solve crossprod tcrossprod rowMeans rowSums
#' @importFrom Rdpack reprompt
#' @importFrom rlang list2
#' @importFrom glmnet glmnet
#' @importFrom ROCR prediction performance
#'
spar <- function(x, y, family = gaussian("identity"), model = NULL, rp = NULL,
                 screencoef = NULL, xval = NULL, yval = NULL, nnu = 20L, nus = NULL,
                 nummods = 20L, measure = c("deviance","mse","mae","class","1-auc"),
                 avg_type = c("link", "response"),
                 parallel = FALSE, inds = NULL, RPMs = NULL, seed = NULL, ...) {
  # Set up and checks ----
  measure <- match.arg(measure)
  stopifnot("Length of y does not fit nrow(x)." = length(y) == nrow(x))
  stopifnot("Response y must be numeric." = is.numeric(y))
  # Ensure back compatibility ----
  args <- list(...)
  arg_list <- check_and_set_args(args, x, y, family, model,
                                 screencoef, rp,  measure)
  model <- arg_list$model; rp <- arg_list$rp
  screencoef <- arg_list$screencoef; measure <- arg_list$measure
  avg_type <- match.arg(avg_type)

  # Call SPAR algorithm ----
  res <- spar_algorithm(x = x, y = y,
                        family = family,
                        model = model, rp = rp, screencoef = screencoef,
                        xval = xval, yval = yval,
                        nnu = nnu, nus = nus,
                        nummods = nummods,
                        measure = measure,
                        avg_type = avg_type,
                        inds = inds, RPMs = RPMs,
                        parallel = parallel,
                        seed = seed)
  return(res)

}

#' @rdname spar
#' @examples
#' spar_res <- spareg(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10, 15, 20, 25, 30))
#' @aliases spar
#' @export
spareg <- spar

fit_spar_models <- function(x, y, family, model, rp, screencoef,
                            nnu, nus, nummods, measure, avg_type,
                            inds = NULL, RPMs = NULL, parallel = FALSE, seed = NULL) {
  # Set up and checks -----
  p <- ncol(x)
  n <- nrow(x)

  # Scaling the x matrix -----
  xcenter <- colMeans(x)
  xscale <- apply(x, 2, sd)
  if (is.null(inds) || is.null(RPMs)) {
    actual_p <- sum(xscale > 0)
    xz <- scale(x[, xscale > 0], center = xcenter[xscale > 0], scale = xscale[xscale > 0])
  } else {
    actual_p <- p
    xscale[xscale == 0] <- 1
    xz <- scale(x, center = xcenter, scale = xscale)
  }
  # Scaling the y vector ----
  if (family$family == "gaussian" && family$link == "identity") {
    ycenter <- mean(y)
    yscale <- sd(y)
  } else {
    ycenter <- 0
    yscale <- 1
  }
  yz <- scale(y, center = ycenter, scale = yscale)

  # Argument names ----
  formal_names <- names(formals(fit_spar_models))
  formal_names_wo_x_y <- setdiff(formal_names, c("x", "y"))
  all_args_wo_x_y <- mget(formal_names_wo_x_y, envir = environment())

  # Set up seed ----
  if (!is.null(seed)) {
    if (parallel && requireNamespace("doRNG", quietly = TRUE)) {
      registerDoRNG <- getNamespace("doRNG")$registerDoRNG
      registerDoRNG(seed = seed)
    } else {
      set.seed(seed)
    }
  }

  # Update model object ----
  all_args_wo_x_y <- mget(formal_names_wo_x_y, envir = environment())
  model <- do.call(function(...)
    model$update_fun(object = model, x = xz, y = yz, ...),
    all_args_wo_x_y)

  # Setup screening ----
  if (!is.null(attr(screencoef, "split_data_prop"))) {
    scr_inds <- sample(n, ceiling(n * attr(screencoef, "split_data_prop")))
    mar_inds <- seq_len(n)[-scr_inds]
  } else {
    mar_inds <- scr_inds <- seq_len(n)
  }
  if (2 * n > p) {
    message("Screening is not performed by default, as 2 * n, the default number of screened variables, is larger than the number of predictors. For performing screening, adjust nscreen in screen_*().")
  }
  if (is.null(attr(screencoef, "nscreen"))) {
    nscreen <- attr(screencoef, "nscreen") <- min(p, 2 * n)
  } else {
    nscreen <- attr(screencoef, "nscreen")
  }

  # Checks for mslow, msup, nscreen ----
  mslow <- attr(rp, "mslow")
  if (is.null(mslow)) mslow <- attr(rp, "mslow") <- ceiling(log(p))
  msup <- attr(rp, "msup")
  if (is.null(msup)) msup <- attr(rp, "msup") <- ceiling(n/2)
  if (!(msup <= nscreen)) {
    message("Provided upper bound on goal dimension of random projection (msup) or its default value (n/2) is larger than nscreen. Setting msup to nscreen.")
    msup <- nscreen
  }
  stopifnot("Provided lower bound on goal dimension of random projection (mslow) or its default value (log(p)) is larger than upper bound (msup)." = mslow <= msup)
  # Perform screening ----
  if (nscreen < p) {
    # First update the object
    all_args_wo_x_y <- mget(formal_names_wo_x_y, envir = environment())
    screencoef <- do.call(function(...)
      screencoef$update_fun(object = screencoef, x = xz, y = yz, ...),
      all_args_wo_x_y)
    # Then compute the screening coefs
    scr_coef <- do.call(function(...)
      screencoef$generate_fun(object = screencoef,
                              x = xz[scr_inds, ],
                              y = yz[scr_inds, ], ...), all_args_wo_x_y)
    inc_probs <- abs(scr_coef)
    max_inc_probs <- max(inc_probs)
    inc_probs <- inc_probs / max_inc_probs
    attr(screencoef, "inc_prob") <- inc_probs
    if (attr(screencoef, "type") == "prob" && sum(inc_probs > 0) < nscreen) {
      warning(
        sprintf("The number of variables with non-zero screening coefficients (%i) is less than the number of variables to screen (%i). Probabilistic screening with nscreen variables is performed anyway, but some of the variables with a zero inclusion probability will be randomly added to the set of screened variables. Alternatively, nscreen can be lowered in screen_*().",
                sum(inc_probs > 0), nscreen)
      )
    }
  } else {
    scr_coef <- NULL
  }
  attr(screencoef, "importance") <- scr_coef

  # Update RP -----
  all_args_wo_x_y <- mget(formal_names_wo_x_y, envir = environment())
  rp <- do.call(function(...)
    rp$update_fun(object = rp, x = xz, y = yz, ...), all_args_wo_x_y)

  # Flags for draw RPMs -----
  max_num_mod <- max(nummods)

  drawRPMs <- FALSE
  if (is.null(RPMs)) {
    RPMs <- vector("list", length = max_num_mod)
    drawRPMs <- TRUE
    ms <- sample(seq(floor(mslow), ceiling(msup)), max_num_mod, replace = TRUE)
  }

  drawinds <- FALSE
  if (is.null(inds)) {
    inds <- vector("list", length = max_num_mod)
    drawinds <- TRUE
  }

  all_args_wo_x_y <- mget(formal_names_wo_x_y, envir = environment())
  # SPAR algorithm -----
  marginal_model_function <- function(i) {
    out <- list()
    if (drawinds) {
      if (nscreen < p) {
        ind_use <- switch(attr(screencoef, "type"),
                          "fixed" = order(inc_probs, decreasing = TRUE)[seq_len(nscreen)],
                          "prob" = c(
                            sample(seq_len(actual_p)[inc_probs > 0], min(sum(inc_probs > 0), nscreen), prob = inc_probs[inc_probs > 0]),
                            sample(seq_len(actual_p)[inc_probs == 0], nscreen - min(sum(inc_probs > 0), nscreen))
                          ),
                          stop("Type of screening coef should be fixed or prob.")
        )
      } else {
        ind_use <- seq_len(actual_p)
      }
    } else {
      ind_use <- inds[[i]]
    }
    out$inds <- ind_use
    p_use <- length(ind_use)

    ## RP step -----
    if (drawRPMs) {
      m <- ms[i]
      if (p_use < m) {
        m <- p_use
        RPM <- Matrix::Matrix(diag(1, m), sparse = TRUE)
      } else {
        RPM <- do.call(function(...)
          rp$generate_fun(object = rp, x = xz, y = yz, m = m,
                          included_vector = ind_use, ...), all_args_wo_x_y)
      }
    } else {
      RPM <- RPMs[[i]]
      if (drawinds) RPM <- RPM[, c(ind_use)]
      ## Update RPM w data
      RPM <- do.call(function(...)
        rp$update_rpm_w_data(rpm = RPM, object = rp, x = xz, y = yz,
                             included_vector = ind_use, ...), all_args_wo_x_y)
    }
    out$RPMs <- RPM

    # Marginal model
    znew <- Matrix::tcrossprod(xz[mar_inds, ind_use], RPM)
    res <- do.call(function(...)
      model$generate_fun(y = yz[mar_inds], z = znew, object = model, ...),
      all_args_wo_x_y)
    out$intercepts <- res$intercept
    out$betas_std_m <- as(numeric(actual_p), "sparseMatrix")
    out$betas_std_m[ind_use] <- crossprod(RPM, res$gammas)
    out
  }

  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel execution. Please install it using install.packages('foreach').")
    }
    foreach <- getNamespace("foreach")$foreach
    `%dopar%` <- getNamespace("foreach")$`%dopar%`
    `%do%` <- getNamespace("foreach")$`%do%`
    getDoParRegistered <- getNamespace("foreach")$getDoParRegistered
    getDoParName <- getNamespace("foreach")$getDoParName
    getDoParWorkers <- getNamespace("foreach")$getDoParWorkers

    if (!getDoParRegistered()) {
      message('Warning: No doPar backend. Executing SPAR algorithm sequentially. For using parallelization, please register backend and rerun.')
      `%d%` <- `%do%`
    } else {
      message('Using ', getDoParName(), ' with ', getDoParWorkers(), ' workers')
      `%d%` <- `%dopar%`
    }
    i <- NULL
    res_all <- foreach(i = seq_len(max_num_mod), .verbose = FALSE, .packages = "spareg", .errorhandling = "stop") %d% {
      marginal_model_function(i = i)
    }
  } else {
    res_all <- lapply(seq_len(max_num_mod), marginal_model_function)
  }

  if (drawRPMs) RPMs <- lapply(res_all, "[[", "RPMs")
  if (drawinds) inds <- lapply(res_all, "[[", "inds")
  intercepts <- sapply(res_all, "[[", "intercepts")
  betas_std <- Reduce("cbind2", lapply(res_all, "[[", "betas_std_m"))
  if (is.null(colnames(x))) {
    rownames(betas_std) <- paste0("V", seq_len(ncol(x[, xscale > 0, drop = FALSE])))
  } else {
    rownames(betas_std) <- colnames(x[, xscale > 0, drop = FALSE])
  }
  if (is.null(nus)) {
    if (nnu > 1) {
      nus <- unname(c(0, quantile(abs(betas_std@x), probs = seq_len(nnu - 1) / (nnu - 1))))
    } else {
      nus <- 0
    }
  }
  # Return fitted objects
  return(list(
    betas_std = betas_std,
    intercepts = intercepts,
    scr_coef = scr_coef,
    inds = inds,
    RPMs = RPMs,
    nus = nus,
    xcenter = xcenter,
    xscale = xscale,
    ycenter = ycenter,
    yscale = yscale,
    avg_type = avg_type,
    measure = measure,
    family = family,
    model = model,
    rp = rp,
    screencoef = screencoef,
    x_rows_for_fitting_marginal_models = if (!is.null(attr(screencoef, "split_data_prop"))) mar_inds else NULL
  ))
}

validate_spar <- function(fitted_objects, xval, yval, nus, nummods, measure, avg_type) {
  p <- length(fitted_objects$xscale)
  n <- nrow(xval)
  # Get validation measure function
  val.meas <- get_val_measure_function(measure, fitted_objects$family)

  # Initialize validation results
  val_res <- data.frame(nnu = NULL, nu = NULL, nummod = NULL,
                        numactive = NULL, measure = NULL)

  # Loop over nummods
  tabnummodres <- lapply(nummods, function(nummod) {
    tabres <- lapply(seq_along(nus), function(l) {
      thresh <- nus[l]
      tmp_coef <- fitted_objects$betas_std[, seq_len(nummod), drop = FALSE]
      tmp_coef[abs(tmp_coef) <= thresh] <- 0
      tmp_beta <- Matrix(0, nrow = p, ncol = nummod)
      tmp_beta[fitted_objects$xscale > 0, ] <- fitted_objects$yscale * tmp_coef / (fitted_objects$xscale[fitted_objects$xscale > 0])
      if (avg_type == "link") {
        beta_hat <- rowMeans(tmp_beta)
        alpha_hat <- mean(fitted_objects$intercepts[seq_len(nummod)]) + (fitted_objects$ycenter - sum(fitted_objects$xcenter * beta_hat))
        eta_hat <- xval %*% beta_hat + alpha_hat
        val_measure <- val.meas(yval, eta_hat = eta_hat)
        numactive <- sum(beta_hat != 0)
      } else {
        tmp_intercept <- fitted_objects$intercepts[seq_len(nummod)] + drop(fitted_objects$ycenter - crossprod(fitted_objects$xcenter, tmp_beta))
        eta_hat <- sweep((xval %*% tmp_beta), tmp_intercept, MARGIN = 2, FUN = "+")
        y_hat <- rowMeans(fitted_objects$family$linkinv(as.matrix(eta_hat)))
        val_measure <- val.meas(yval, y_hat = y_hat)
        numactive <- sum(rowSums(tmp_beta != 0) > 0)
      }
      c(nnu = l, nu = unname(thresh), nummod = nummod, measure = val_measure, numactive = numactive)
    })
    out <- do.call("rbind", tabres)
    colnames(out) <- c("nnu", "nu", "nummod", "measure", "numactive")
    out
  })

  val_res <- do.call("rbind.data.frame", tabnummodres)
  return(val_res)
}

spar_algorithm <- function(x, y, family, model, rp, screencoef,
                           xval = NULL, yval = NULL,
                           nnu, nus,
                           nummods, measure,
                           avg_type,
                           inds = NULL, RPMs = NULL,
                           parallel = FALSE,
                           seed = NULL){
  # Start fitting SPAR algorithm -----
  res <- fit_spar_models(x = x, y = y, family = family, model = model, rp = rp, screencoef = screencoef,
                         nnu = nnu, nus = nus, nummods = nummods, measure = measure, avg_type = avg_type,
                         inds = inds, RPMs = RPMs, parallel = parallel, seed = seed)

  betas <- Matrix(0, ncol(x), max(nummods), sparse = TRUE)
  betas[res$xscale > 0, ] <- res$betas_std
  if (is.null(colnames(x))) {
    rownames(betas) <- paste0("V", seq_len(ncol(x)))
  } else {
    rownames(betas) <- colnames(x)
  }
  # Compute validation measures -----
  if (is.null(xval)) xval <- x
  if (is.null(yval)) yval <- y
  val_res <- validate_spar(fitted_objects = res,
                           xval = xval, yval = yval,
                           nus = res$nus,
                           nummods = nummods,
                           measure = measure,
                           avg_type = avg_type)
  res[["val_res"]] <- val_res
  res[["betas"]] <- betas
  res[["betas_std"]] <- NULL
  attr(res,"class") <- "spar"

  return(res)
}


#' Coef Method for \code{'spar'} Object
#'
#' Extracts coefficients from \code{'spar'} object
#' @param object result of [spar] of class \code{'spar'}.
#' @param nummod number of models used to form coefficients; value with minimal
#'        validation \code{measure} is used if not provided.
#' @param nu threshold level used to compute the coefficients; value with minimal
#'        validation \code{measure} is used if not provided.
#' @param aggregate character, one of \code{c("mean", "median", "none")},
#'        giving the method of aggregating
#'        the coefficients over the marginal models. If set to \code{"none"},
#'        the coefficients are not aggregated over the marginal models and a
#'        matrix of coefficients, one column for each marginal model, is returned.
#'        Otherwise
#'        the coefficients are aggregated using the specified method (mean or median).
#'        Defaults to mean aggregation.
#' @param ... further arguments passed to or from other methods.
#' @return Returns an object of class  \code{'coefspar'}, which is a list with elements
#' \itemize{
#'  \item \code{intercept}: The average intercept value or vector intercepts (one for
#'  each marginal model) if \code{aggregate = "none"}.
#'  \item \code{beta}: Vector of length \eqn{p} of averaged coefficients or a
#'        \eqn{p} x \code{max(nummods)} matrix of coefficients if \code{agregate = "none"}.
#'  \item \code{nummod}: Number of models based on which the coefficients are computed.
#'  \item \code{nu}:  Threshold value based on which the coefficients are computed.
#' }
#' @seealso [print.coefspar], [summary.coefspar]
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' coef(spar_res)
#' coef(spar_res, aggregate = "median")
#' coef(spar_res, aggregate = "none")
#' coef(spar_res, nummod = 5, nu = 0)
#' @export
#' @method coef spar
coef.spar <- function(object,
                      nummod = NULL,
                      nu = NULL,
                      aggregate = c("mean", "median", "none"),
                      ...) {
  aggregate <- match.arg(aggregate)
  if (is.null(nummod) & is.null(nu)) {
    best_ind <- which.min(object$val_res$measure)
    par <- object$val_res[best_ind,]
    nummod <- par$nummod
    nu <- par$nu
  } else if (is.null(nummod)) {
    if (!nu %in% object$val_res$nu) {
      stop("Nu needs to be among the previously considered values when nummod is not provided!")
    }
    tmp_val_res <- object$val_res[object$val_res$nu==nu,]
    nummod <- tmp_val_res$nummod[which.min(tmp_val_res$measure)]
  } else if (is.null(nu)) {
    if (!nummod %in% object$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when nu is not provided!")
    }
    tmp_val_res <- object$val_res[object$val_res$nummod==nummod,]
    nu <- tmp_val_res$nu[which.min(tmp_val_res$measure)]
  } else {
    if (length(nummod)!=1 | length(nu)!=1) {
      stop("Length of nummod and nu must be 1!")
    }
  }

  if (nummod > ncol(object$betas)) {
    warning("Number of models is too high, maximum of number of models use to fit the models is used instead!")
    nummod <- ncol(object$betas)
  }

  # calc for chosen parameters
  final_coef <- object$betas[object$xscale>0, seq_len(nummod), drop=FALSE]
  final_coef[abs(final_coef) < nu] <- 0
  p <- length(object$xscale)
  if (aggregate == "none") {
    beta <- matrix(0, nrow = p, ncol = nummod)
    beta_std <- final_coef
    beta[object$xscale>0,] <- as.matrix(object$yscale *
                                          beta_std/(object$xscale[object$xscale>0]))
    rownames(beta) <- rownames(object$betas)
    colnames(beta) <- paste0("Model_", seq_len(nummod))
    intercept <- drop(object$ycenter + object$intercepts[seq_len(nummod)] -
                        crossprod(object$xcenter, beta))
    names(intercept) <- colnames(beta)
  } else {
    avg_fun <- switch(aggregate,
                      "mean" = function(x) mean(x, na.rm = TRUE),
                      "median" = function(x) median(x, na.rm = TRUE),
                      "Aggregration method not implemented.")
    beta <- numeric(p)
    beta_std <- apply(final_coef, 1, avg_fun)
    beta[object$xscale>0] <- object$yscale * beta_std/(object$xscale[object$xscale>0])
    names(beta) <- rownames(object$betas)
    intercept <- object$ycenter + avg_fun(object$intercepts[seq_len(nummod)]) - sum(object$xcenter*beta)
    names(intercept) <- "(Intercept)"
  }
  res <- list(intercept = intercept,
              beta = beta,
              nummod = nummod,
              nu = nu)
  class(res) <- "coefspar"
  best_ind <- which.min(object$val_res$measure)
  par <- object$val_res[best_ind,]
  attr(res, "M_best") <- par$nummod
  attr(res, "nu_best") <- par$nu
  attr(res, "M_nu_combination") <- ifelse(nummod == par$nummod &  nu == par$nu,
                                          "best", "given")
  attr(res, "aggregate") <- aggregate
  attr(res, "parent_object") <- class(object)
  return(res)
}

#' Print Method for \code{'coefspar'} Object
#'
#' Prints a summary of coefficients from a `'coefspar'` object, including the selected \eqn{M} and \eqn{\nu},
#' and the number of active variables.
#'
#' @param x A  \code{'coefspar'} object.
#' @param digits The number of significant digits for numeric output. Default: `4L`.
#' @param show The number of coefficients to display. Default: `6L`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 2000, ntest = 100)
#' spar_res <- spareg(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' coef(spar_res)
#' coef(spar_res, aggregate = "median")
#' coef(spar_res, aggregate = "none")
#' print(coef(spar_res), show = 10L, digits = 6L)
#' @export
#' @method print coefspar
print.coefspar <- function(x, digits = 4L, show = 6L, ...) {
  cat(sprintf("Coefficients from %s object:\n", attr(x, "parent_object")))
  cat(sprintf("Based on the %s selection\n\n",
              switch(attr(x, "M_nu_combination"),
                     "best" = "*best rule* (min error)",
                     "1se" = "*1se rule*",
                     "given"= "*given* nummod and nu")))

  if (!is.null(attr(x, "aggregate"))) {
    cat(sprintf("Aggregation method over models: %s\n",
                attr(x, "aggregate")))
  }
  cat("Selected combination: ", x$nummod, " models, threshold = ",
      x$nu, "\n")
  cat("\n")
  if (attr(x, "M_nu_combination") != "best") {
    cat("Best combination overall (min error):",
        attr(x, "M_best"), "models, threshold =",
        attr(x, "nu_best"), "\n")
  }
  if (attr(x, "M_nu_combination") != "1se") {
    if (!is.null(attr(x, "M_1se")) & !is.null(attr(x, "nu_1se"))) {
      cat("1se combination:",  attr(x, "M_1se"), "models, threshold =",
          attr(x, "nu_1se"), "\n")
    }
  }
  # cat("\n")
  # Coefficient vectors or matrices
  cat("Coefficients:\n")
  if (attr(x, "aggregate") == "none") {
    ## No aggregation ----
    coefs <- rbind("(Intercept)" = x$intercept, x$beta)
    shown <- head(coefs, show)
    # Assign default names if unnamed
    # if (any(is.null(rownames(shown))) || any(rownames(shown) == "")) {
    #   rownames(shown)[-1] <- paste0("V", seq_along(shown - 1))
    # }
    # Print header and values
    print(noquote(format(round(shown, digits), nsmall = digits,
                         justify = "right")))
    # Add inline ...
    if (nrow(coefs) > show) {
      cat("...", sprintf("(%d rows not shown)\n\n", nrow(coefs) - show))
    }
    cat("Number of active variables: \n")
    no_non_zero_coefs <- paste0(colSums(x$beta != 0), "/", nrow(x$beta))
    names(no_non_zero_coefs) <- colnames(x$beta)
    print(noquote(format(no_non_zero_coefs, nsmall = digits,
                         justify = "right")))
    cat("\n")
  } else {
    ## Aggregation ----
    coefs <- c(x$intercept, x$beta)
    shown <- head(coefs, show)
    # Assign default names if unnamed
    if (is.null(names(shown)) || any(names(shown) == "")) {
      names(shown) <- paste0("V", seq_along(shown) - 1)
      names(shown)[1] <- "(Intercept)"
    }

    # Print header and values
    print(noquote(format(round(shown, digits), nsmall = digits,
                         justify = "right")))

    # Add inline ...
    if (length(coefs) > show) {
      cat("...", sprintf("(%d coefficients not shown)\n\n",
                         length(coefs) - show))
    }

    cat("Number of active variables: ",
        paste0(sum(x$beta != 0), "/", length(x$beta)), "\n\n")
  }

  invisible(x)
}

#' Summary Method for \code{'coefspar'} Object
#'
#' Provides a detailed summary of a \code{coefspar} object, including coefficient statistics.
#'
#' @param object An object of class \code{coefspar}.
#' @param digits Number of digits to be printed. Defaults to `4L`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return
#' Invisibly returns input \code{object}.
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 2000, ntest = 100)
#' spar_res <- spareg(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' summary(coef(spar_res))
#' summary(coef(spar_res, aggregate = "none"))
#'
#' @export
#' @method summary coefspar
summary.coefspar <- function(object, digits = 4L, ...) {
  stopifnot(inherits(object, "coefspar"))
  cat(sprintf(
    "Summary of coefficients from %s object:\n", attr(object, "parent_object")))
  cat(sprintf("Based on the %s selection\n\n",
              switch(attr(object, "M_nu_combination"),
                     "best" = "*best rule* (min error)",
                     "1se" = "*1se rule*",
                     "given"= "*given* nummod and nu")))

  if (!is.null(attr(object, "aggregate"))) {
    cat(sprintf("Aggregation method over models: %s\n",
                attr(object, "aggregate")))
  }
  cat("Selected combination: ", object$nummod, " models, threshold = ",
      object$nu, "\n")
  cat("\n")
  if (attr(object, "M_nu_combination") != "best") {
    cat("Best combination overall (min error):",
        attr(object, "M_best"), "models, threshold =",
        attr(object, "nu_best"), "\n\n")
  }
  if (attr(object, "M_nu_combination") != "1se") {
    if (!is.null(attr(object, "M_1se")) & !is.null(attr(object, "nu_1se"))) {
      cat("1se combination:",  attr(object, "M_1se"), "models, threshold =",
          attr(object, "nu_1se"), "\n\n")
    }
  }
  if (attr(object, "aggregate") == "none") {
    cat("Number of coefficients equal to zero across all models:",
        paste0(sum(rowMeans(object$beta) == 0), "/", nrow(object$beta)), "\n")
    cat("Number of coefficients non-zero across all models:",
        paste0(sum(rowMeans(object$beta != 0) == 1), "/", nrow(object$beta)), "\n\n")
  } else {
    cat("Number of active coefficients:",
        paste0(sum(object$beta != 0), "/", length(object$beta)), "\n\n")
  }
  if (attr(object, "aggregate") == "none") {
    cat("Intercept:\n  ")
  }
  print(round(object$intercept, digits))
  cat("Coefficient summary (beta):\n")
  print(summary(object$beta, digits = digits))
  invisible(object)
}




#' Predict Method for \code{'spar.cv'} Object
#'
#' Predict responses for new predictors from \code{'spar'} object
#' @param object result of spar function of class  \code{'spar'}.
#' @param xnew matrix of new predictor variables; must have same number of columns as  \code{x}.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param nummod number of models used to form coefficients; value with minimal validation measure is used if not provided.
#' @param nu threshold level used to form coefficients; value with minimal validation measure is used if not provided.
#' @param aggregate character one of c("mean", "median");
#'        the aggregation over the ensembles is done using the specified method (mean or median).
#'        Defaults to mean aggregation.
#' @param ... further arguments passed to or from other methods
#' @return Returns a vector of predictions.
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' pred <- predict(spar_res, xnew = example_data$xtest)
#' @export
predict.spar <- function(object,
                         xnew = NULL,
                         type     = c("response", "link"),
                         avg_type = c("link","response"),
                         nummod   = NULL,
                         nu       = NULL,
                         aggregate = c("mean", "median"),
                         ...) {
  if (is.null(xnew)) {
    stop("No 'xnew' provided. This 'spar' object does not retain training data.",
         "Please provide xnew explicitly.",
         "If you want to predict in-sample, use the original data used for fitting the model as xnew.\n")
  }

  if (ncol(xnew) != length(object$xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  aggregate <- match.arg(aggregate)
  # object$family <- eval(parse(text = object$family))

  if (avg_type != object$avg_type && object$family$link != "identity") {
    warning("The best model combination was selected for ",
            paste0(object$avg_type, " averaging, but avg_type = ", avg_type,
                   " is used for prediction. This may lead to suboptimal results."))
  }

  coefs_avg <- coef(object, nummod, nu, aggregate = aggregate)
  eta <- as.numeric(xnew %*% coefs_avg$beta + coefs_avg$intercept)
  if (avg_type == "link") {
    res <- if (type == "link") eta else object$family$linkinv(eta)
  } else {
    if (type == "link") {
      res <- eta
    } else {
      avg_fun <- switch(aggregate,
                        "mean" = function(x) mean(x, na.rm = TRUE),
                        "median" = function(x) median(x, na.rm = TRUE),
                        "Aggregration method not implemented.")
      coefs_all <- coef(object, nummod, nu, aggregate = "none")
      eta_all <- sweep(xnew %*% coefs_all$beta, coefs_all$intercept, MARGIN = 2, FUN = "+")
      preds <- object$family$linkinv(eta_all)
      res <- apply(preds, 1, avg_fun)
    }
  }
  return(res)
}
#'
#' Plot Method for \code{'spar'} Object
#'
#' @description
#' Creates diagnostic plots for a `'spar'`] object, including:
#'   - Validation measure vs. `nus` or `nummods`.
#'   - Number of active variables vs. `nus` or `nummods`.
#'   - Residuals vs. fitted values.
#'   - Coefficient heatmap.
#'
#'
#' @param x A \code{'spar'} object.
#' @param plot_type The type of plot. Options:
#'   - `"val_measure"` (default): Validation measure vs. `nus` or `nummods`.
#'   - `"val_numactive"`: Number of active variables vs. `nus` or `nummods`.
#'   - `"res_vs_fitted"`: Residuals vs. fitted values.
#'   - `"coefs"`: Heatmap of coefficients.
#' @param plot_along The variable to plot along the x-axis. Options:
#'   - `"nu"` (default): Threshold values.
#'   - `"nummod"`: Number of marginal models.
#'   Ignored if `plot_type = "res_vs_fitted"` or `plot_type = "coefs"`.
#' @param nummod The number of models to fix when `plot_along = "nu"`. If `NULL`, the optimal value is used.
#'   The value will be used in \code{\link{predict.spar}} when \code{plot_type="res_vs_fitted"}.
#' @param nu The threshold value to fix when `plot_along = "nummod"`. If `NULL`, the optimal value is used.
#'   The value will be used in \code{\link{predict.spar}} when \code{plot_type="res_vs_fitted"}.
#' @param xfit The predictor data for fitted values (required if `plot_type = "res_vs_fitted"`).
#' @param yfit The response data for fitted values (required if `plot_type = "res_vs_fitted"`).
#' @param prange A vector of length 2 specifying the range of predictors to plot (for `plot_type = "coefs"`).
#'   Default: `c(1, p)`.
#' @param coef_order An optional vector specifying the order of coefficients for `plot_type = "coefs"`.
#'   Default: `1:p` (original order).
#' @param digits The number of significant digits for axis labels. Default: `2L`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Returns a \code{'\link[ggplot2:ggplot]{ggplot2::ggplot}'}  object.
#'
#' @import ggplot2
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' plot(spar_res)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nu", nummod = 10)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nu", nummod = 10)
#' plot(spar_res, plot_type = "res_vs_fitted",  xfit = example_data$xtest,
#'   yfit = example_data$ytest)
#' plot(spar_res, plot_type = "coefs", prange = c(1,400))
#
#' @export
#'
plot.spar <- function(x,
                      plot_type = c("val_measure","val_numactive","res_vs_fitted","coefs"),
                      plot_along = c("nu","nummod"),
                      nummod = NULL,
                      nu = NULL,
                      xfit = NULL,
                      yfit = NULL,
                      prange = NULL,
                      coef_order = NULL,
                      digits = 2L, ...) {
  spar_res <- x
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  mynummod <- nummod
  if (plot_type == "res_vs_fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res_vs_fitted plot!")
    }
    pred <- predict(spar_res, xnew = xfit, nummod = nummod, nu = nu,
                    type = "response")
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,
                                             residuals=yfit-pred),
                           ggplot2::aes(x=.data$fitted,y=.data$residuals)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  } else if (plot_type == "val_measure") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$measure)]
        tmp_title <- "Optimal~number~of~models~M[best]=="
      } else {
        tmp_title <- "Given~number~of~models~M=="
      }
      tmp_df <- spar_res$val_res[spar_res$val_res$nummod==mynummod, ]
      ind_min <- which.min(tmp_df$measure)

      tmp_title_text <- paste0(tmp_title,mynummod)

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x=.data$nu,y=.data$measure)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(x=expression(nu),
                      y=spar_res$measure) +
        ggplot2::geom_vline(xintercept = tmp_df$nu[ind_min],linetype=2,linewidth=0.5)+
        ggplot2::geom_point(data=data.frame(x=tmp_df$nu[ind_min],
                                            y=tmp_df$measure[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red") +
        ggplot2::ggtitle(parse(text = tmp_title_text))
    } else {
      if (is.null(nu)) {
        nu <- spar_res$val_res$nu[which.min(spar_res$val_res$measure)]
        tmp_title <- "Optimal~threshold~nu[best]=="
      } else {
        tmp_title <- "Given~threshold~nu=="
      }

      nu_grid <- max(spar_res$nus[spar_res$nus <= nu])
      tmp_df <- spar_res$val_res[which(spar_res$val_res$nu==nu_grid), ]
      ind_min <- which.min(tmp_df$measure)
      tmp_title_text <- paste0(tmp_title, round(nu, 3))


      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x=.data$nummod,y=.data$measure)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(y=spar_res$measure) +
        ggplot2::geom_vline(xintercept = tmp_df$nummod[ind_min],linetype=2,linewidth=0.5)+
        ggplot2::geom_point(data = data.frame(x = tmp_df$nummod[ind_min],
                                              y = tmp_df$measure[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        scale_x_continuous(breaks=seq(min(tmp_df$nummod), max(tmp_df$nummod),1),
                           minor_breaks = NULL)+
        ggplot2::ggtitle(parse(text = tmp_title_text))
      }
  } else if (plot_type=="val_numactive") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$measure)]
        tmp_title <- "Optimal~number~of~models~M[best]=="
      } else {
        tmp_title <- "Given~number~of~models~M=="
      }
      tmp_df <- spar_res$val_res[spar_res$val_res$nummod==mynummod, ]
      ind_min <- which.min(tmp_df$measure)
      tmp_title_text <- paste0(tmp_title,mynummod)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nu,y=.data$numactive)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        # ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),labels=round(spar_res$val_res$nu,3)) +
        #ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),
        #                            labels=formatC(spar_res$val_res$nu[seq(1,nrow(spar_res$val_res),1)],
        #                                           format = "e", digits = digits)) +
        ggplot2::labs(x=expression(nu)) +
        ggplot2::geom_vline(xintercept = tmp_df$nu[ind_min],linetype=2,linewidth=0.5)+
        ggplot2::geom_point(data=data.frame(x=tmp_df$nu[ind_min],y=tmp_df$numactive[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(parse(text = tmp_title_text))
    } else {
      if (is.null(nu)) {
        nu <- spar_res$val_res$nu[which.min(spar_res$val_res$measure)]
        tmp_title <- "Optimal~threshold~nu[best]=="
      } else {
        tmp_title <- "Given~threshold~nu=="
      }

      nu_grid <- max(spar_res$nus[spar_res$nus <= nu])
      tmp_df <- spar_res$val_res[which(spar_res$val_res$nu==nu_grid), ]
      ind_min <- which.min(tmp_df$measure)
      tmp_title_text <- paste0(tmp_title, round(nu, 3))

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nummod,y=.data$numactive)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::geom_point(
          data=data.frame(x=tmp_df$nummod[ind_min],
                          y=tmp_df$numactive[ind_min]),
          ggplot2::aes(x = .data$x,y=.data$y),col="red")+
        ggplot2::geom_vline(xintercept = tmp_df$nummod[ind_min],linetype=2,linewidth=0.5)+
        scale_x_continuous(breaks=seq(min(tmp_df$nummod), max(tmp_df$nummod),1),
                           minor_breaks = NULL)+
        ggplot2::ggtitle(parse(text = tmp_title_text))
    }
  } else if (plot_type=="coefs") {
    p <- nrow(spar_res$betas)
    nummod <- ncol(spar_res$betas)
    if (is.null(prange)) {
      prange <- c(1,p)
    }
    if (is.null(coef_order)) {
      coef_order <- 1:p
    }

    tmp_mat <- data.frame(t(apply(as.matrix(spar_res$betas)[coef_order,],1,
                                  function(row)row[order(abs(row),decreasing = TRUE)])),
                          predictor=1:p)
    colnames(tmp_mat) <- c(seq_len(nummod),"predictor")
    tmp_df <- reshape(tmp_mat, idvar = "predictor",
                      varying = seq_len(nummod),
                      v.names = "value",
                      timevar = "marginal model",
                      direction = "long")

    tmp_df$`marginal model` <- as.numeric(tmp_df$`marginal model`)

    mrange <- max(Matrix::rowSums(spar_res$betas != 0))
    res <- ggplot2::ggplot(tmp_df,ggplot2::aes(x=.data$predictor,
                                               y=.data$`marginal model`,
                                               fill=.data$value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2() +
      ggplot2::coord_cartesian(xlim=prange,ylim=c(1,mrange)) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Index of marginal model") +
      ggplot2::theme(panel.border = ggplot2::element_blank())

  } else {
    res <- NULL
  }
  return(res)
}

#' Print summary of \code{'spar'} Object
#'
#' @description
#' Prints a summary of a [`spar`] object, including the validation
#' measure, optimal \eqn{M} and \eqn{\nu}, and the number of active predictors.
#'
#' @param x A \code{'spar'} object.
#' @param digits The number of significant digits for numeric output. Default: `4L`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return
#' Invisibly returns the input object `x`.
#'
#' @examples
#' example_data <- simulate_spareg_data(n = 100, p = 400, ntest = 100)
#' spar_res <- spar(example_data$x, example_data$y, xval = example_data$xtest,
#'   yval = example_data$ytest, nummods=c(5, 10))
#' print(spar_res)
#' @export
print.spar <- function(x, digits = 4L,...) {
  mycoef <- coef(x)
  beta <- mycoef$beta
  measure <- x$val_res$measure[mycoef$nu == x$val_res$nu &
                                 mycoef$nummod == x$val_res$nummod ]
  if (nrow(x$val_res) == 1) {
    cat(sprintf("spar object: \nValidation measure (%s) of %s reached for nummod=%d,
              nu=%s leading to %d / %d active predictors.\n",
                x$measure,
                formatC(measure,digits = 2,format = "e"),
                mycoef$nummod, formatC(mycoef$nu,digits = digits,format = "e"),
                sum(beta!=0),length(beta)))
  } else {
    cat(sprintf("spar object:\nSmallest validation measure (%s) of %s reached for nummod=%d,
              nu=%s leading to %d / %d active predictors.\n",
                x$measure,
                formatC(measure,digits = 2,format = "e"),
                mycoef$nummod, formatC(mycoef$nu,digits = digits,format = "e"),
                sum(beta!=0),length(beta)))
  }
  cat("Summary of those non-zero coefficients:\n")
  print(summary(beta[beta!=0]))
}
