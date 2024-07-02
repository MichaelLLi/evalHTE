#' Evaluate Heterogeneous Treatment Effects
#' @param treatment Treatment variable
#' @param form a formula object that takes the form \code{y ~ T + x1 + x2 + ...}.
#' @param data
#'   A data frame that contains the outcome \code{y} and the treatment \code{T}.
#' @param algorithms
#'   List of machine learning algorithms to be used.
#' @param budget The maximum percentage of population that can be treated under the budget constraint.
#' @param c_threshold The threshold for determining the treatment assignment under the ITR. Default is 0, which assumes all units with a score greater than 0 will receive treatment.
#' @param n_folds
#'   Number of cross-validation folds. Default is 5.
#' @param split_ratio
#'   Split ratio between train and test set under sample splitting. Default is 0.
#' @param ngates
#'   The number of groups to separate the data into. The groups are determined by tau. Default is 5.
#' @param preProcess caret parameter
#' @param weights caret parameter
#' @param trControl caret parameter
#' @param tuneGrid caret parameter
#' @param tuneLength caret parameter
#' @param user_model A user-defined function to estimate heterogeneous treatment effects. The function should take the data as input and return a model to estimate the HTE.
#' @param SL_library A list of machine learning algorithms to be used in the super learner.
#' @param ... Additional arguments passed to \code{caret::train}.
#' @import dplyr
#' @importFrom rlang !! sym
#' @export
#' @return An object of \code{hte} class
evaluate_hte <- function(
    treatment,
    form,
    data,
    algorithms,
    budget,
    n_folds = 5,
    split_ratio = 0,
    ngates = 5,
    c_threshold = 0,
    preProcess = NULL,
    weights = NULL,
    trControl = caret::trainControl(method = "none"),
    tuneGrid = NULL,
    tuneLength = ifelse(trControl$method == "none", 1, 3),
    user_model = NULL,
    SL_library = NULL,
    ...
) {

  # estimate HTE
  fit <- evalITR::estimate_itr(
    treatment = treatment,
    form = form,
    data = data,
    algorithms = algorithms,
    budget = budget,
    n_folds = n_folds,
    split_ratio = split_ratio,
    ngates = ngates,
    c_threshold = c_threshold,
    preProcess = preProcess,
    weights = weights,
    trControl = trControl,
    tuneGrid = tuneGrid,
    tuneLength = tuneLength,
    user_model = user_model,
    SL_library = SL_library,
    ...)

  # parameters
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
    df  = data
    Ycv = data[, outcome]
    Tcv = data[, treatment]

    estimates <- list(
      Ycv = Ycv,
      Tcv = Tcv,
      algorithms = "user-defined")
    cv   <- FALSE

    # compute qoi
    qoi   <- vector("list", length = length(outcome))
    qoi   <- compute_qoi_user(
      user_model, Tcv, Ycv, data, ngates, budget, ...)

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
    cat('Conduct hypothesis tests for GATEs unde sample splitting ...\n')

    # create empty lists to for consistcv and hetcv
    consist <- list()
    het <- list()

    # run consistency and heterogeneity tests for each model
    for (i in algorithms) {

      consist[[i]] <- consist.test(
        T   = Tcv,
        tau = fit_ml[[i]]$tau,
        Y   = Ycv,
        ngates = ngates)

      het[[i]] <- het.test(
        T   = Tcv,
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
    cat('Conduct hypothesis tests for GATEs unde cross-validation ...\n')

    # create empty lists to for consistcv and hetcv
    consistcv <- list()
    hetcv <- list()

    # run consistency and heterogeneity tests for each model
    for (i in algorithms) {

      consistcv[[i]] <- consistcv.test(
        T   = Tcv,
        tau = gettaucv(fit)[[i]],
        Y   = Ycv,
        ind = indcv,
        ngates = ngates)

      hetcv[[i]] <- hetcv.test(
        T   = Tcv,
        tau = gettaucv(fit)[[i]],
        Y   = Ycv,
        ind = indcv,
        ngates = ngates)

    }

  # return a list of consistcv and hetcv
  out <- list(consistcv = consistcv, hetcv = hetcv)
  }

  class(out) <- c("test_hte", class(out))

  return(out)

}


utils::globalVariables(c("T", "aupec", "sd", "pval", "Pval", "aupec.y", "fraction", "AUPECmin", "AUPECmax", ".", "fit", "out", "pape", "alg", "papep", "papd", "type", "gate", "group", "qnorm", "vec", "Y", "algorithm", "statistic", "p.value", "GATEcv", "RATEmin", "RATEpoint", "Type", "best_ind", "best_rate", "consist.test", "consistcv.test", "est", "gettaucv", "het.test", "hetcv.test", "rate", "value", "map", "model.matrix", "quantile","rnorm"))

