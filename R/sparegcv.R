#' Sparse Projected Averaged Regression with Cross-Validation
#'
#' Apply Sparse Projected Averaged Regression to High-Dimensional Data, where the
#' number of models and the threshold parameter is chosen using a cross-validation
#' procedure.
#'
#' @param x n x p numeric matrix of predictor variables.
#' @param y quantitative response vector of length n.
#' @param family  a \code{'\link[stats]{family}'} object used for the marginal
#'        generalized linear model; defaults to \code{gaussian("identity")}.
#' @param model function creating a \code{'sparmodel'} object;
#'   defaults to \code{spar_glm()} for gaussian family with identity link and to
#'   \code{spar_glmnet()} for all other family-link combinations.
#' @param rp function creating a \code{'randomprojection'} object. If not specified by user it
#'   defaults to \code{rp_cw(data = TRUE)}.
#' @param screencoef function creating a \code{'screeningcoef'} object.  If not specified by user it
#'   defaults to no screening.
#' @param nfolds number of folds to use for cross-validation; should be at least 2, defaults to 10.
#' @param nnu number of different threshold values \eqn{\nu} to consider for thresholding;
#'        ignored when \code{nus} is provided; defaults to 20.
#' @param nus optional vector of \eqn{\nu}'s to consider for thresholding;
#'         if not provided, \code{nnu} values ranging from 0 to the maximum absolute
#'         marginal coefficient are used.
#' @param nummods vector of numbers of marginal models to consider for
#'        validation; defaults to \code{c(20)}.
#' @param measure loss to use for validation; defaults to \code{"deviance"}
#'        available for all families. Other options are \code{"mse"} or \code{"mae"}
#'         (between responses and predicted means, for all families),
#'         \code{"class"} (misclassification error) and
#'         \code{"1-auc"} (one minus area under the ROC curve) both just for
#'         binomial family.
#' @param avg_type type of averaging the marginal models; either on link (default)
#'        or on response level. This is used in computing the validation measure.
#' @param parallel assuming a parallel backend is loaded and available, a
#'        logical indicating whether the function should use it in parallelizing the
#'        estimation of the marginal models. Defaults to FALSE.
#' @param seed integer seed to be set at the beginning of the SPAR algorithm. Default to NULL, in which case no seed is set.
#' @param fast_fit character, one of \code{c("fix_rpm_and_inds", "fix_rpm", "none")}, indicating whether to use the same random projection matrix and/or the same indices for screening across folds. Defaults to \code{"fix_rpm_and_inds"}.
#' @param ... further arguments mainly to ensure back-compatibility
#' @returns object of class \code{'spar.cv'} with elements
#' \itemize{
#'  \item \code{val_res} a \code{data.frame} with CV results for each fold and for each element of nus and nummods
#'  \item \code{fitted_objects} list of fitted objects for the whole data (first element) and for each fold, each element is a list with elements
#'    \itemize{
#'  \item \code{betas_std} p x  \code{max(nummods)} sparse matrix of class
#'   \code{'\link[Matrix:dgCMatrix-class]{Matrix::dgCMatrix}'} containing the
#'   standardized coefficients from each marginal model computed with the SPAR
#'   algorithm on the training data.
#'  \item \code{intercepts} in each marginal model, vector of length \code{max(nummods)}
#'    computed with the SPAR algorithm on the  training data.
#'  \item \code{scr_coef} p-vector of coefficients used for screening for standardized predictors
#'  \item \code{inds} list of index-vectors corresponding to variables kept after
#'  screening in each marginal model of length  \code{max(nummods)}. If kept fixed, only the first element corresponding to the indicators computed on the whole data is populated.
#'  \item \code{RPMs} list of projection matrices used in each marginal model of length \code{max(nummods)}.  If kept fixed, only the first element corresponding to the projection matrices generated once on the whole data is populated.
#'  \item \code{ycenter} empirical mean of response vector in training data
#'  \item \code{yscale} empirical standard deviation of response vector in training data
#'. \item \code{xcenter} p-vector of empirical means of predictor variables  in training data
#'  \item \code{xscale} p-vector of empirical standard deviations of  predictor variables  in training data
#'  \item \code{seed} integer seed used at the beginning of the algorithm. In each fold the seed is modified as \code{seed + fold_id}.
#'  }
#'  \item \code{nus} vector of \eqn{\nu}'s considered for thresholding
#'  \item \code{nummods} vector of numbers of marginal models considered for validation
#'  \item \code{family}  a character corresponding to \link[stats]{family}  object used for the marginal generalized linear model e.g.,
#'  \code{"gaussian(identity)"}
#'  \item \code{measure} character, type of validation measure used
#'  \item \code{avg_type} character, averaging type for computing the validation measure
#'  \item \code{rp} an object of class \code{'randomprojection'}
#'  \item \code{screencoef} an object of class \code{'screeningcoef'}
#'  \item \code{model} an object of class \code{'sparmodel'}
#'  \item \code{fast_fit} character, one of \code{c("fix_rpm_and_inds", "fix_rpm", "none")}, indicating whether the same random projection matrix and/or the same indices for screening across folds were used.
#'  \item \code{seed} integer seed used at the beginning of the algorithm.
#'  }
#' @examples
#' \donttest{
#' example_data <- simulate_spareg_data(n = 80, p = 200, ntest = 100)
#' spar_res <- spar.cv(example_data$x, example_data$y, nfolds = 3L,
#'   rp = rp_gaussian(), nummods = c(5, 10))
#' spar_res
#' }
#' @seealso [spar], [coef.spar.cv], [predict.spar.cv], [plot.spar.cv], [print.spar.cv]
#' @aliases spareg.cv
#' @export
spar.cv <- function(x, y, family = gaussian("identity"), model = spar_glmnet(),
                    rp = NULL, screencoef = NULL, nfolds = 10,
                    nnu = 20, nus = NULL, nummods = c(20),
                    measure = c("deviance","mse","mae","class","1-auc"),
                    avg_type = c("link","response"),
                    parallel = FALSE, seed = NULL,
                    fast_fit = c("fix_rpm_and_inds", "fix_rpm", "none"),
                    ...) {
  # Set up and checks ----
  n <- length(y)
  stopifnot("Length of y does not fit nrow(x)." = n == nrow(x))
  stopifnot("Response y must be numeric." = is.numeric(y))
  stopifnot("nfolds must be at least 2." = nfolds >= 2)
  stopifnot("matrix" %in% class(x) |"data.frame" %in% class(x))
  x <- as.matrix(x)
  if (!is.numeric(x[1,1])) {
    stop("There are non-numeric data entries, numerical matrix needed!")
  }

  measure <- match.arg(measure)
  fast_fit <- match.arg(fast_fit)
  avg_type <- match.arg(avg_type)

  # Ensure back compatibility ----
  args <- list(...)
  arg_list <- check_and_set_args(args, x, y, family, model,
                                 screencoef, rp,  measure)
  model <- arg_list$model; rp <- arg_list$rp
  screencoef <- arg_list$screencoef; measure <- arg_list$measure

  # Run initial spar algorithm ----
  ## Precompute nus if not provided ----
  #if (is.null(nus) | !(fast_fit == "none")) {
  if (fast_fit == "fix_rpm") {
    temp_screencoef <- NULL
    temp_arg_list <- check_and_set_args(
      args, x, y, family, model,
      temp_screencoef, rp,  measure)
    model <- temp_arg_list$model; rp <- temp_arg_list$rp
    temp_screencoef <- temp_arg_list$screencoef;
    measure <- temp_arg_list$measure
  } else {
    temp_screencoef <- screencoef
  }
  temp_fit <- fit_spar_models(
    x = x, y = y, family = family, model = model, rp = rp,
    screencoef = temp_screencoef, nnu = nnu, nus = nus,
    nummods = nummods, measure = measure, avg_type = avg_type,
    parallel = parallel, seed = seed
  )
  temp_val_res <- validate_spar(temp_fit, x, y, temp_fit$nus,
                                nummods, measure, avg_type)

  if (is.null(nus)) nus <- temp_fit$nus
  if (!(fast_fit == "none")) {
    RPMs <- temp_fit$RPMs
  } else {
    RPMs <- NULL
  }
  if (fast_fit == "fix_rpm_and_inds") {
    inds <- temp_fit$inds
  } else {
    inds <- NULL
  }
  #}

  # Folds
  folds <- sample(cut(seq_len(n), breaks = nfolds, labels=FALSE))

  # Initialize storage for results
  all_val_res <- list()
  ## First element is the initial fit
  all_val_res[[1]] <- cbind(fold = 0, temp_val_res)

  all_fitted_objects <- vector(length = nfolds + 1, mode = "list")
  ## First element is the initial fit
  all_fitted_objects[[1]] <-  list(
    betas_std = temp_fit$betas_std,
    intercepts = temp_fit$intercepts,
    scr_coef = temp_fit$scr_coef,
    inds = temp_fit$inds, RPMs = temp_fit$RPMs,
    xcenter = temp_fit$xcenter, xscale = temp_fit$xscale,
    ycenter = temp_fit$ycenter, yscale = temp_fit$yscale
  )

  # Loop over folds ----
  for (fold in seq_len(nfolds)) {
    # Split data
    if (is.null(seed)) {
      seed_fold <- NULL
    } else {
      seed_fold <- seed + fold
    }
    x_train <- x[folds != fold, temp_fit$xscale>0]
    y_train <- y[folds != fold]
    x_val <- x[folds == fold, temp_fit$xscale>0]
    y_val <- y[folds == fold]

    # Fit models on training data
    fitted_objects <- fit_spar_models(
      x_train, y_train, family, model, rp, screencoef,
      nnu = nnu, nus = nus, nummods = nummods,
      inds = inds, RPMs = RPMs,
      measure = measure, avg_type = avg_type,
      parallel = parallel, seed = seed_fold
    )

    # Validate on held-out data
    val_res <- validate_spar(fitted_objects, x_val, y_val, nus,
                             nummods, measure, avg_type)
    all_val_res[[fold + 1]] <- cbind(fold, val_res)
    all_fitted_objects[[fold + 1]] <-
      list("betas_std"  = fitted_objects$"betas_std",
           "intercepts" = fitted_objects$"intercepts",
           "scr_coef"   = fitted_objects$"scr_coef",
           "inds"       = if (fast_fit == "fix_rpm_and_inds") NULL else fitted_objects$"inds",
           "RPMs"       = if (fast_fit %in% c("fix_rpm_and_inds", "fix_rpm")) NULL else fitted_objects$"RPMs",
           "xcenter"    = fitted_objects$"xcenter",
           "xscale"     = fitted_objects$"xscale",
           "ycenter"    = fitted_objects$"ycenter",
           "yscale"     = fitted_objects$"yscale",
           "seed" = fitted_objects$"seed",
           "x_rows_for_fitting_marginal_models" = fitted_objects$"x_rows_for_fitting_marginal_models")
  }

  # Aggregate validation results across folds
  combined_val_res <- do.call(rbind, all_val_res)

  ## We do not refit on whole data with best parameters
  res <- list(
    val_res = combined_val_res,
    fitted_objects = all_fitted_objects,
    nus = nus,
    nummods = nummods,
    family = paste0(family$family, "(", family$link, ")"),
    measure = measure,
    avg_type = avg_type,
    nfolds = nfolds,
    rp = rp,
    screencoef = screencoef,
    model = model,
    fast_fit = fast_fit,
    seed = seed
  )
  attr(res,"class") <- "spar.cv"
  return(res)
}

#' @rdname spar.cv
#' @examples
#' \donttest{
#' spar_res <- spareg.cv(example_data$x, example_data$y,
#'   nummods=c(5, 10, 15, 20, 25, 30))
#' }
#' @aliases spar.cv
#' @export
spareg.cv <- spar.cv



#' Plot Method for \code{'spar.cv'} Object
#'
#' @description
#' Plot cross-validation measure or number of active variables over different thresholds or number
#' of models of \code{'spar.cv'} object, produce a residuals vs fitted plot,
#' or a plot of the estimated coefficients in each marginal model, sorted by their absolute value.
#'
#' @param x result of [spar.cv] function of class  \code{'spar.cv'}.
#' @param plot_type one of  \code{c("val_measure","val_numactive", "res_vs_fitted", "coefs")}.
#' @param plot_along one of  \code{c("nu","nummod")}.
#' @param nummod fixed value for  \code{nummod} when  \code{plot_along="nu"} for
#'  \code{plot_type="val_measure"} or  \code{"val_numactive"};
#' @param nu fixed value for \eqn{\nu} when  \code{plot_along="nummod"}
#' for  \code{plot_type="val_measure"} or  \code{"val_numactive"}; same as for \code{\link{predict.spar.cv}} when  \code{plot_type="res_vs_fitted"}.
#' @param xfit optional vector of fitted values to be used for  \code{plot_type="res_vs_fitted"}; if not provided, the fitted values are computed using the best parameters. Argument is used only if \code{fast_fit = "fix_rpm_and_inds"} and \code{plot_type = "res_vs_fitted"}.
#' @param yfit optional vector of response values to be used for  \code{plot_type="res_vs_fitted"}; if not provided, the response values are computed using the best parameters. Argument is used only if \code{fast_fit = "fix_rpm_and_inds"} and \code{plot_type = "res_vs_fitted"}.
#' @param opt_par one of \code{c("1se","best")}, chooses whether to select the best pair of \code{nus} and \code{nummods} according to cross-validated (CV) measure, or the sparsest solution within one sd of that optimal CV measure. Argument is used only if \code{fast_fit = "fix_rpm_and_inds"} and \code{plot_type = "res_vs_fitted"}.
#' @param prange optional vector of length 2 indicating the range of predictors to be plotted for \code{plot_type = "coefs"}; defaults to \code{c(1, p)} where \code{p} is the number of predictors. Argument can be used only if \code{fast_fit = "fix_rpm_and_inds"}.
#' @param coef_order optional vector of length \code{p} indicating the order of predictors to be plotted for \code{plot_type = "coefs"}; defaults to \code{seq_len(p)} where \code{p} is the number of predictors. Argument can be used only if \code{fast_fit = "fix_rpm_and_inds"}.
#' @param digits number of significant digits to be displayed in the axis; defaults to 2L.
#' @param ... further arguments passed to or from other methods
#' @return \code{'\link[ggplot2:ggplot]{ggplot2::ggplot}'}  object
#' @import ggplot2
#' @examples
#' \donttest{
#' example_data <- simulate_spareg_data(n = 80, p = 200, ntest = 100)
#' spar_res <- spar.cv(example_data$x, example_data$y, nfolds = 3L,
#'   screencoef = screen_cor(), rp = rp_gaussian(), nummods = c(5, 10))
#' plot(spar_res)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_measure", plot_along = "nu", nummod = 10)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nummod", nu = 0)
#' plot(spar_res, plot_type = "val_numactive",  plot_along = "nu", nummod = 10)
#' }
#' @export
plot.spar.cv <- function(x,
                         plot_type = c("val_measure","val_numactive","res_vs_fitted","coefs"),
                         plot_along = c("nu","nummod"),
                         nummod = NULL,
                         nu = NULL,
                         xfit = NULL,
                         yfit = NULL,
                         opt_par = c("best","1se"),
                         prange = NULL,
                         coef_order = NULL, digits = 2L, ...)  {
  spar_res <- x
  plot_type <- match.arg(plot_type)
  if (plot_type %in% c("res_vs_fitted", "coefs") && spar_res$fast_fit != "fix_rpm_and_inds") {
    stop("Plot types 'res_vs_fitted' and 'coefs' are not yet implemented!")
  }
  plot_along <- match.arg(plot_along)
  opt_par <- match.arg(opt_par)
  mynummod <- nummod
  my_val_sum <- compute_val_summary(spar_res$val_res)
  colnames(my_val_sum)[match(c("mean_measure", "mean_numactive"),colnames(my_val_sum))] <- c("Meas", "numactive")


  if (plot_type=="val_measure") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- my_val_sum[my_val_sum$nummod==mynummod, ]
      ind_min <- which.min(tmp_df$Meas)
      allowed_ind <- tmp_df$Meas<=tmp_df$Meas[ind_min]+
        tmp_df$sd_measure[ind_min]
      ind_1se <- which(min(tmp_df$numactive[allowed_ind]) ==
                         tmp_df$numactive)

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x = .data$nu,y = .data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(x=expression(nu),y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nu[ind_min],
                                            y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$Meas-.data$sd_measure,
                                          ymax=.data$Meas+.data$sd_measure),
                             alpha=0.2,linetype=2,show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color="red",show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nu[ind_min],
                                                  tmp_df$nu[ind_1se]),
                                            y = c(tmp_df$Meas[ind_min],tmp_df$Meas[ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nu[ind_min],
                          y = tmp_df$Meas[ind_min] + tmp_df$sd_measure[ind_min],
                          xend = tmp_df$nu[ind_1se],
                          yend = tmp_df$Meas[ind_min] + tmp_df$sd_measure[ind_min],
                          color = 2, linetype = 2)
    } else {
      if (is.null(nu)) {
        nu <- my_val_sum$nu[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- my_val_sum[my_val_sum$nu==nu, ]
      ind_min <- which.min(tmp_df$Meas)
      allowed_ind <- tmp_df$Meas<=tmp_df$Meas[ind_min]+
        tmp_df$sd_measure[ind_min]
      ind_1se <- match(min(tmp_df$numactive[allowed_ind]),
                       tmp_df$numactive)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nummod,y=.data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),list(txt=tmp_title,v=round(nu,3)))) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$Meas-.data$sd_measure,
                                          ymax=.data$Meas+.data$sd_measure),
                             alpha=0.2,linetype=2,show.legend = FALSE)+
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color="red",show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$Meas[ind_min],tmp_df$Meas[allowed_ind][ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nummod[ind_min],
                          y = tmp_df$Meas[ind_min] + tmp_df$sd_measure[ind_min],
                          xend = tmp_df$nummod[allowed_ind][ind_1se],
                          yend = tmp_df$Meas[ind_min] + tmp_df$sd_measure[ind_min],
                          color=2,linetype=2) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=seq(min(tmp_df$nummod), max(tmp_df$nummod),1),
                                    minor_breaks = NULL)
    }
  }

  if (plot_type=="val_numactive") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- my_val_sum[my_val_sum$nummod==mynummod, ]
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sd_measure[ind_min]
      ind_1se <- which.min(tmp_df$numactive[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nu,y=.data$numactive)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x=expression(nu)) +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nu[ind_min],tmp_df$nu[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numactive[ind_min],tmp_df$numactive[allowed_ind][ind_1se]))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(nu)) {
        nu <- my_val_sum$nu[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- my_val_sum[my_val_sum$nu==nu, ]
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sd_measure[ind_min]
      ind_1se <- which.min(tmp_df$numactive[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x=.data$nummod,y=.data$numactive)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() + ggplot2::theme_bw() +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numactive[ind_min],tmp_df$numactive[allowed_ind][ind_1se]))) +
        ggplot2::scale_x_continuous(breaks=seq(min(tmp_df$nummod), max(tmp_df$nummod),1),
                                    minor_breaks = NULL)+
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),list(txt=tmp_title,v=round(nu,3))))

    }
  }

  if (plot_type == "coefs") {
    betas_std <- spar_res$fitted_objects[[1]]$betas_std
    nummod <- ncol(betas_std)

    xscale <- spar_res$fitted_objects[[1]]$xscale
    p <- length(xscale)
    beta <- matrix(0, nrow = p, ncol = nummod)
    beta[xscale>0,] <- as.matrix(betas_std)

    if (is.null(prange)) {
      prange <- c(1, p)
    }
    if (is.null(coef_order)) {
      coef_order <- seq_len(p)
    }

    tmp_mat <- data.frame(t(apply(beta[coef_order,],1,
                                  function(row)row[order(abs(row),decreasing = TRUE)])),
                          predictor=seq_len(p))
    colnames(tmp_mat) <- c(1:nummod,"predictor")

    tmp_df <- reshape(tmp_mat, idvar = "predictor",
                      varying = seq_len(nummod),
                      v.names = "value",
                      timevar = "marginal model",
                      direction = "long")

    tmp_df$`marginal model` <- as.numeric(tmp_df$`marginal model`)

    mrange <- max(Matrix::rowSums(beta != 0))

    res <- ggplot2::ggplot(tmp_df,ggplot2::aes(x=.data$predictor,
                                               y=.data$`marginal model`,
                                               fill=.data$value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2() +
      ggplot2::coord_cartesian(xlim=prange,ylim=c(1,mrange)) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Index of marginal model") +
      ggplot2::theme(panel.border = ggplot2::element_blank())
  }

  if (plot_type=="res_vs_fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res_vs_fitted plot!")
    }
    pred <- predict(spar_res, xfit, opt_par=opt_par, nummod=nummod, nu=nu)
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,residuals=yfit-pred),
                           ggplot2::aes(x=.data$fitted,y=.data$residuals)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  }
  return(res)
}


#' Print Method for \code{'spar.cv'} Object
#'
#' Print summary of \code{'spar.cv'} object
#' @param x result of  [spar.cv] function of class  \code{'spar.cv'}.
#' @param digits integer digits to be printed, defaults to 4L.
#' @param ... further arguments passed to or from other methods
#' @return text summary
#' @examples
#' \donttest{
#' example_data <- simulate_spareg_data(n = 80, p = 200, ntest = 100)
#' spar_res <- spareg.cv(example_data$x, example_data$y, nfolds = 3L,
#'   screencoef = screen_cor(), rp = rp_gaussian(), nummods = c(5, 10))
#' print(spar_res)
#' }
#' @export
print.spar.cv <- function(x, digits = 4L, ...) {
  val_sum <- compute_val_summary(x$val_res)
  if (nrow(val_sum) == 1) {
    cat(sprintf(
      "spar.cv object: \nCV measure (%s) of  %s (averaged over the %i fold) reached for nummod=%d, nu=%s.",
      #leading to %d / %d active predictors.\n",
      x$measure, formatC(min(val_sum$mean_measure),digits = digits,format = "e"),
      x$nfolds,
      val_sum$nummod,
      formatC(val_sum$nu,digits = digits,format = "e")
    )
    )
  } else {
    tmp_df <- val_sum
    ind_min <- which.min(tmp_df$mean_measure)

    allowed_ind <- tmp_df$mean_measure<tmp_df$mean_measure[ind_min]+tmp_df$sd_measure[ind_min]
    ind_1se <- match(min(tmp_df$mean_numactive[allowed_ind]),
                     tmp_df$mean_numactive)
    cat(sprintf(
      "spar.cv object:\n\nSmallest CV measure (%s) of %s reached for nummod=%d, nu=%s. \n",
      # leading  to %d / %d active predictors (averaged over folds).",
      x$measure,  formatC(min(val_sum$mean_measure),digits = digits,format = "e"),
      tmp_df$nummod[ind_min],
      formatC(tmp_df$nu[ind_min],digits = digits,format = "e")
      #      sum(my_best$beta!=0),length(my_best$beta)
    )
    )
    cat(sprintf(
      "\nSparsest coefficient within one standard error of best CV measure (%s)  \
        of %s  reached for nummod=%d, nu=%s.\n",
      x$measure, formatC(tmp_df$mean_measure[ind_1se],digits = digits,format = "e"),
      tmp_df$nummod[ind_1se],
      formatC(tmp_df$nu[ind_1se],digits = digits,format = "e")
    ))
  }
}

#' Coef Method for \code{'spar.cv'} Object
#'
#' Extract coefficients from \code{'spar.cv'} object
#' Only allowed if \code{fast_fit = "fix_rpm_and_inds"} is used when calling \code{spar.cv}.
#' In this case the ensemble obtained on the whole data is used for the
#' combination of threshold and number of models using the best or the 1se identified through cross-validation.
#' Otherwise coefficients cannot be extracted without refitting on whole data with
#'  best or 1se parameters. To be able to extract coefficients in case \code{fast_fit = "fix_rpm"}
#'  or \code{fast_fit = "none"}, use \code{spar()} to refit with the desired (nu, M) combination.
#' @param object result of [spar.cv] function of class \code{'spar.cv'}. Refitting
#' can be done using \code{get_model().}
#' @param nummod optional number of models used to form coefficients
#' @param nu optional threshold level used to form coefficients
#' @param opt_par one of \code{c("1se","best")}, chooses whether to select the
#'        best pair of \code{nus} and \code{nummods} according to cross-validated
#'        (CV) measure, or the sparsest solution within one sd of that optimal
#'        CV measure;
#'        ignored when \code{nummod} and \code{nu} are given
#' @param aggregate character one of c("mean", "median", "none"). If set to "none"
#'        the coefficients are not aggregated over the marginal models, otherwise
#'        the coefficients are aggregated using the specified method (mean or median).
#'        Defaults to mean aggregation.
#' @param ... further arguments passed to or from other methods
#' @return List with elements
#' \itemize{
#'  \item \code{intercept} intercept value
#'  \item \code{beta} vector of length p of averaged coefficients
#'  \item \code{nummod} number of models based on which the coefficient is computed
#'  \item \code{nu}  threshold based on which the coefficient is computed
#' }
#' @examples
#' \donttest{
#' example_data <- simulate_spareg_data(n = 80, p = 200, ntest = 100)
#' spar_res <- spar.cv(example_data$x, example_data$y, nfolds = 3L,
#'   nummods = c(5, 10))
#' coef(spar_res)
#' }
#' @seealso [predict.spar.cv], [get_model.spar.cv]
#' @export

coef.spar.cv <- function(object,
                         nummod = NULL,
                         nu = NULL,
                         opt_par = c("best","1se"),
                         aggregate = c("mean", "median", "none"),
                         ...) {
  if (object$fast_fit != "fix_rpm_and_inds") {
    stop("Coefficients cannot be extracted without refitting on whole data with best or 1se parameters. To be able to extract coefficients, set fast_fit to 'fix_rpm_and_inds' when calling spar.cv or use spar() to refit with the desired (nu, M) combination. ")
  }
  opt_nunum <- match.arg(opt_par)
  aggregate <- match.arg(aggregate)
  given_pars <- !is.null(nummod) & !is.null(nu)

  # best model
  val_sum <- compute_val_summary(object$val_res)
  best_ind <- which.min(val_sum$mean_measure)
  parbest <- val_sum[best_ind,]

  # 1se model
  allowed_ind <- val_sum$mean_measure<val_sum$mean_measure[best_ind]+
    val_sum$sd_measure[best_ind]
  ind_1cv <- which.min(val_sum$mean_numactive[allowed_ind])
  par1se <- val_sum[allowed_ind,][ind_1cv,]

  if (is.null(nummod) & is.null(nu)) {

    if (opt_nunum == "1se") {
      nummod <- par1se$nummod
      nu <- par1se$nu
    } else {
      nummod <- parbest$nummod
      nu <- parbest$nu
    }

  } else if (is.null(nummod)) {
    if (!is.null(nu) & !nu %in% val_sum$nu) {
      stop("nu needs to be among the previously considered values when nummod is not provided!")
    }
    tmp_val_sum <- val_sum[val_sum$nu==nu,]
    if (opt_nunum=="1se") {
      allowed_ind <- tmp_val_sum$mean_measure<tmp_val_sum$mean_measure[best_ind]+tmp_val_sum$sd_measure[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mean_numactive[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mean_measure),]
    }
    nummod <- par$nummod
  } else if (is.null(nu)) {
    if (!is.null(nummod) & !nummod %in% object$val_res$nummod) {
      stop("Number of models needs to be among the previously considered values when nu is not provided!")
    }
    tmp_val_sum <- val_sum[val_sum$nummod==nummod,]
    if (opt_nunum=="1se") {
      allowed_ind <- tmp_val_sum$mean_measure<tmp_val_sum$mean_measure[best_ind]+tmp_val_sum$sd_measure[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mean_numactive[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mean_measure),]
    }
    nu <- par$nu
  } else {
    if (length(nummod)!=1 | length(nu)!=1) {
      stop("Length of nummod and nu must be 1!")
    }
  }



  # calc for chosen parameters
  betas_std <- object$fitted_objects[[1]]$betas_std
  if (nummod > ncol(betas_std)) {
    warning("Number of models is too high, maximum of fitted is used instead!")
    nummod <- ncol(betas_std)
  }

  xcenter <- object$fitted_objects[[1]]$xcenter
  ycenter <- object$fitted_objects[[1]]$ycenter
  xscale <- object$fitted_objects[[1]]$xscale
  yscale <- object$fitted_objects[[1]]$yscale
  intercepts <- object$fitted_objects[[1]]$intercepts
  val_res <- object$val_res[object$val_res$fold == 0, ]

  final_coef <- betas_std[, seq_len(nummod), drop=FALSE]
  final_coef[abs(final_coef) < nu] <- 0
  final_coef[abs(final_coef) < nu] <- 0
  p <- length(object$fitted_objects[[1]]$xscale)
  if (aggregate == "none") {
    beta <- matrix(0, nrow = p, ncol = nummod)
    beta_std <- final_coef
    beta[xscale>0,] <- as.matrix(yscale * beta_std/(xscale[xscale>0]))
    rownames(beta) <- rownames(final_coef)
    colnames(beta) <- paste0("Model_", seq_len(nummod))
    intercept <- drop(ycenter + intercepts[seq_len(nummod)] -
                        crossprod(xcenter, beta))
    names(intercept) <- colnames(beta)
  } else {
    avg_fun <- switch(aggregate,
                      "mean" = function(x) mean(x, na.rm = TRUE),
                      "median" = function(x) median(x, na.rm = TRUE),
                      "Aggregration method not implemeneted")
    beta <- numeric(p)
    beta_std <- apply(final_coef, 1, avg_fun)
    beta[xscale>0] <- yscale * beta_std/(xscale[xscale>0])
    names(beta) <- rownames(final_coef)
    intercept <- ycenter + avg_fun(intercepts[seq_len(nummod)]) - sum(xcenter*beta)
    names(intercept) <- "(Intercept)"
  }
  res <- list(intercept = intercept,
              beta = beta,
              nummod = nummod,
              nu = nu)

  class(res) <- "coefspar"
  best_ind <- which.min(val_res$measure)
  par <- val_res[best_ind,]
  attr(res, "M_best") <- parbest$nummod
  attr(res, "nu_best") <- parbest$nu
  attr(res, "M_1se") <- par1se$nummod
  attr(res, "nu_1se") <- par1se$nu
  attr(res, "M_nu_combination") <-
    ifelse(given_pars, "given", opt_par)
  attr(res, "aggregate") <- aggregate
  attr(res, "parent_object") <- class(object)
  return(res)
}

#' Predict Method for \code{'spar.cv'} Object
#'
#' Predict responses for new predictors from \code{'spar.cv'} object.
#' Only allowed if \code{fast_fit = "fix_rpm_and_inds"} is used when calling \code{spar.cv}.
#' In this case the ensemble obtained on the whole data is used for the
#' combination of threshold and number of models using the best or the 1se identified through cross-validation.
#' Otherwise predictions cannot be generated without refitting on whole data with
#'  best or 1se parameters. To be able to generate predictions in case \code{fast_fit = "fix_rpm"}
#'  or \code{fast_fit = "none"}, use \code{spar()} to refit with the desired (nu, M) combination.
#' @param object result of spar function of class \code{'spar.cv'}.
#' @param xnew matrix of new predictor variables; must have same number of columns as  \code{x}.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param opt_par one of  \code{c("best","1se")}, chooses whether to select the
#'  best pair of  \code{nus} and  \code{nummods} according to CV measure, or the
#'  sparsest solution within one sd of that optimal CV measure;
#'  ignored when  \code{nummod} and  \code{nu}, or  \code{coef} are given
#' @param nummod number of models used to form coefficients; value with
#' minimal validation  \code{Meas} is used if not provided.
#' @param nu threshold level used to form coefficients; value with minimal
#'  validation  \code{Meas} is used if not provided.
#' @param aggregate character one of c("mean", "median");
#'        the aggregation over the ensembles is done using the specified method (mean or median).
#'        Defaults to mean aggregation.
#' @param ... further arguments passed to or from other methods
#' @return Vector of predictions
#' @examples
#' \donttest{
#' example_data <- simulate_spareg_data(n = 80, p = 200, ntest = 100)
#' spar_res <- spar.cv(example_data$x, example_data$y, nfolds = 3L,
#'   rp = rp_gaussian(), nummods = c(5, 10))
#' pred <- predict(spar_res, example_data$x)
#' }
#' @seealso [coef.spar.cv], [get_model.spar.cv]
#' @export
predict.spar.cv <- function(object,
                            xnew = NULL,
                            type = c("response","link"),
                            avg_type = c("link","response"),
                            opt_par = c("best","1se"),
                            nummod = NULL,
                            nu = NULL,
                            aggregate = c("mean", "median"),
                            ...) {
  if (object$fast_fit != "fix_rpm_and_inds") {
    stop("Predictions cannot be generated without refitting on whole data with best or 1se parameters. To be able to generate predictions, set fast_fit to 'fix_rpm_and_inds' when calling spar.cv or use spar() to refit with the desired (nu, M) combination. ")
  }

  if (is.null(xnew)) {
    stop("No 'xnew' provided. This 'spar.cv' object does not retain training data. ",
         "Please provide xnew explicitly. ",
         "If you want to predict in-sample, use the original data used for fitting the model as xnew.")
  }

  xscale <- object$fitted_objects[[1]]$xscale
  if (ncol(xnew)!=length(xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  aggregate <- match.arg(aggregate)
  opt_par <- match.arg(opt_par)

  object$family <- eval(parse(text = object$family))

  coefs_avg <- coef(object, nummod, nu, opt_par = opt_par,
                    aggregate = aggregate)
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
      eta_all <- sweep(xnew %*% coefs_all$beta, coefs_all$intercept,
                       MARGIN = 2, FUN = "+")
      preds <- object$family$linkinv(eta_all)
      res <- apply(preds, 1, avg_fun)
    }
  }

  return(res)
}



