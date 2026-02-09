#' Evaluate Heterogeneous Treatment Effects
#' @param treatment Treatment variable
#' @param form a formula object that takes the form \code{y ~ D + x1 + x2 + ...}.
#' @param data A data frame that contains the outcome \code{y} and the treatment \code{D}.
#' @param algorithms List of machine learning algorithms to be used.
#' @param n_folds Number of cross-validation folds. Default is 5.
#' @param split_ratio Split ratio between train and test set under sample splitting. Default is 0.
#' @param ngates The number of groups to separate the data into. The groups are determined by tau. Default is 5.
#' @param preProcess caret parameter
#' @param weights caret parameter
#' @param trControl caret parameter
#' @param tuneGrid caret parameter
#' @param tuneLength caret parameter
#' @param user_model A user-defined function to estimate heterogeneous treatment effects.
#' @param SL_library A list of machine learning algorithms to be used in the super learner.
#' @param meta_learner The type of meta-learner to use (e.g., "slearner", "tlearner"). Default is "slearner".
#' @param ... Additional arguments passed to \code{caret::train}.
#' @import dplyr
#' @importFrom rlang !! sym
#' @export
#' @return An object of \code{hte} class
estimate_hte <- function(
    treatment,
    form,
    data,
    algorithms,
    n_folds = 5,
    split_ratio = 0,
    ngates = 5,
    preProcess = NULL,
    weights = NULL,
    trControl = caret::trainControl(method = "none"),
    tuneGrid = NULL,
    tuneLength = ifelse(trControl$method == "none", 1, 3),
    user_model = NULL,
    SL_library = NULL,
    meta_learner = "slearner",
    ...
) {

  # estimate HTE
  fit <- evalITR::estimate_itr(
    treatment = treatment,
    form = form,
    data = data,
    algorithms = algorithms,
    budget = 1, # no budget constraint
    n_folds = n_folds,
    split_ratio = split_ratio,
    ngates = ngates,
    preProcess = preProcess,
    weights = weights,
    trControl = trControl,
    tuneGrid = tuneGrid,
    tuneLength = tuneLength,
    user_model = user_model,
    SL_library = SL_library,
    meta_learner = meta_learner,
    ...)

  # return the fit
  return(fit)
}


#' Evaluate Heterogeneous Treatment Effects
#' @param fit Fitted model. Usually an output from \code{estimate_hte}.
#' @param ... Additional arguments passed to the function.
#' @import dplyr
#' @importFrom rlang !! sym
#' @export
#' @return An object of \code{hte} class
evaluate_hte <- function(
    fit,
    ...
) {

  # parameters
  user_model <- fit$user_model

  # initialize output lists
  out_algs <- list()
  out_user  <- list()

  # estimate HTE from ML algorithms
  if(!is.null(fit)){

    # estimate HTE from the fitted model
    estimates  <- fit$estimates
    cv         <- estimates$params$cv
    df         <- fit$df
    algorithms <- fit$df$algorithms
    outcome    <- fit$df$outcome

    # compute qoi
    qoi <- vector("list", length = length(outcome))
    qoi <- compute_qoi(estimates, algorithms)

    # store the results
    out_algs <- list(
      qoi = qoi,
      cv = cv,
      df = df,
      estimates = estimates)
  }

  # get HTE from the user-defined function
  if(!is.null(user_model)){
    # Corrected variable mapping from the fit object
    data_df    <- fit$df$data
    outcome    <- fit$df$outcome
    treatment  <- fit$df$treatment
    ngates     <- fit$estimates$params$ngates
    budget     <- fit$estimates$budget

    # Extract CV vectors
    Ycv <- data_df[[outcome]]
    Tcv <- data_df[[treatment]]

    estimates <- list(
      Ycv = Ycv,
      Tcv = Tcv,
      algorithms = "user-defined")
    cv   <- FALSE

    # compute qoi
    qoi   <- vector("list", length = length(outcome))
    qoi   <- compute_qoi_user(
      user_model, Tcv, Ycv, data_df, ngates, budget, ...)

    # store the results
    out_user <- list(
      qoi = qoi,
      cv = cv,
      df = df,
      estimates = estimates)
  }

  # store the results
  out <- list(
    fit = fit,
    out_algs = out_algs,
    out_user = out_user)

  class(out) <- c("hte", class(out))

  return(out)

}

#' Conduct hypothesis tests
#' @param model Fitted model. Usually an output from \code{evaluate_hte}.
#' @param nsim Number of Monte Carlo simulations used to simulate the null distributions. Default is 1000.
#' @param ... Further arguments passed to the function.
#' @return An object of \code{test_itr} class
#' @export
test_itr <- function(
    model,
    nsim = 1000,
    ...
) {

  fit <- model$fit

  # test parameters
  estimates  <- fit$estimates
  cv         <- estimates$params$cv
  fit_ml     <- estimates$fit_ml
  Tcv        <- estimates$Tcv
  Ycv        <- estimates$Ycv
  indcv      <- estimates$indcv
  n_folds    <- estimates$params$n_folds
  ngates     <- estimates$params$ngates
  algorithms <- fit$df$algorithms
  outcome    <- fit$df$outcome

  # caret and rlearner parameters
  caret_algorithms <- estimates$params$caret_algorithms
  rlearner_algorithms <- estimates$params$rlearner_algorithms

  # super learner library
  SL_library <- estimates$params$SL_library

  # run tests

  ## =================================
  ## sample splitting
  ## =================================

  if(cv == FALSE){
    message('Conduct hypothesis tests for GATEs under sample splitting ...\n')

    # create empty lists to for consistcv and hetcv
    consist <- list()
    het <- list()

    # run consistency and heterogeneity tests for each model
    for (i in algorithms) {

      consist[[i]] <- consist.test(
        D   = Tcv,
        tau = fit_ml[[i]]$tau,
        Y   = Ycv,
        ngates = ngates)

      het[[i]] <- het.test(
        D   = Tcv,
        tau = fit_ml[[i]]$tau,
        Y   = Ycv,
        ngates = ngates)
    }


    # return a list of consist and het
    out <- list(consist = consist, het = het)

  }

  ## =================================
  ## cross validation
  ## =================================

  if(cv == TRUE){
    message('Conduct hypothesis tests for GATEs under cross-validation ...\n')

    get_tau_cv_mat <- function(fit, alg, indcv) {
      obj <- fit$estimates$fit_ml[[alg]]
      K <- max(indcv)
      
      # Case 1: obj is a list of fold objects, each has $tau_cv (your case)
      if (is.list(obj) && length(obj) >= K && !is.null(obj[[1]]$tau_cv)) {
        tau_mat <- sapply(1:K, function(k) obj[[k]]$tau_cv)
        return(as.matrix(tau_mat))
      }
      
      # Case 2: obj itself has $tau_cv (sometimes happens in other implementations)
      if (is.list(obj) && !is.null(obj$tau_cv)) {
        tau_mat <- obj$tau_cv
        if (is.list(tau_mat)) tau_mat <- do.call(cbind, tau_mat)
        return(as.matrix(tau_mat))
      }
      
      stop(paste0("Cannot locate tau_cv for algorithm = ", alg,
                  ". Inspect str(fit$estimates$fit_ml[[alg]])."))
    }
    
    consistcv <- list()
    hetcv <- list()
    
    for (alg in algorithms) {
      
      tau_mat <- get_tau_cv_mat(fit, alg, indcv)
      
      # basic sanity checks (fail fast if something is misaligned)
      stopifnot(nrow(tau_mat) == length(Ycv))
      stopifnot(ncol(tau_mat) == max(indcv))
      
      consistcv[[alg]] <- consistcv.test(
        D = Tcv,
        tau = tau_mat,
        Y = Ycv,
        ind = indcv,
        ngates = ngates
      )
      
      hetcv[[alg]] <- hetcv.test(
        D = Tcv,
        tau = tau_mat,
        Y = Ycv,
        ind = indcv,
        ngates = ngates
      )
    }

  # return a list of consistcv and hetcv
  out <- list(consistcv = consistcv, hetcv = hetcv)
  }

  class(out) <- c("test_hte", class(out))

  return(out)

}


utils::globalVariables(c("D", "aupec", "sd", "pval", "Pval", "aupec.y", "fraction", "AUPECmin", "AUPECmax", ".", "fit", "out", "pape", "alg", "papep", "papd", "type", "gate", "group", "qnorm", "vec", "Y", "algorithm", "statistic", "p.value", "GATEcv", "RATEmin", "RATEpoint", "Type", "best_ind", "best_rate", "consist.test", "consistcv.test", "est", "gettaucv", "het.test", "hetcv.test", "rate", "value", "map", "model.matrix", "quantile","rnorm","estimate", "std.deviation", "lower", "upper", "z.score", "conf.low.uniform"))

