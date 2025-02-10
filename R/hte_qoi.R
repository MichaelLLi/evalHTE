#' Compute Quantities of Interest (GATE, GATEcv, URATE)
#' @param fit_obj An output object of \code{estimate_itr} function from the \code{evalITR} package.
#' @param algorithms Machine learning algorithms
#' @importFrom rlang .data
#' @importFrom furrr future_map

compute_qoi <- function(fit_obj, algorithms) {

  ## extract objects
  fit_ml <- fit_obj$fit_ml
  params <- fit_obj$params
  Ycv    <- fit_obj$Ycv
  if (is.factor(Ycv)) {
    Ycv <- as.numeric(as.character(Ycv))
  } 
  Tcv    <- fit_obj$Tcv
  indcv  <- fit_obj$indcv
  cv     <- fit_obj$params$cv

  ## -----------------------------------------
  ## compute quantities under cross validation
  ## -----------------------------------------
  if (cv == TRUE) {

    ## GATE
    GATE <- vector("list", length = length(algorithms))
    for (i in seq_along(algorithms)) {
      tau <- furrr::future_map(fit_ml[[i]], ~.x$tau) %>% do.call(cbind, .)
      tau_cv <- furrr::future_map(fit_ml[[i]], ~.x$tau_cv) %>% do.call(cbind, .)

      ## Compute GATE
      GATE[[i]] <- GATEcv(Tcv, tau_cv, Ycv, indcv, params$ngates)

      ## indicate algorithm
      GATE[[i]]$alg <- algorithms[i]

      ## indicate group number
      GATE[[i]]$group <- seq_along(GATE[[i]]$gate)
    }

    ## URATE not supported
    URATE <- NULL

  }


  ## -----------------------------------------
  ## compute quantities under sample splitting
  ## -----------------------------------------
  if (cv == FALSE) {

    ## GATE
    GATE <- vector("list", length = length(algorithms))
    for (i in seq_along(algorithms)) {

      ## Compute GATE
      GATE[[i]] <- GATE(Tcv, fit_ml[[i]][["tau"]], Ycv, params$ngates)

      ## indicate algorithm
      GATE[[i]]$alg <- algorithms[i]

      ## indicate group number
      GATE[[i]]$group <- seq_along(GATE[[i]]$gate)
    }

    ## URATE
    URATE <- vector("list", length = length(algorithms))
    for (i in seq_along(algorithms)) {

      ## Compute URATE
      URATE[[i]] <- URATE(Tcv, fit_ml[[i]][["tau"]], Ycv)

      ## indicate algorithm
      URATE[[i]]$alg <- algorithms[i]

    }

  }

  out <- list(
        GATE = GATE,
        URATE = URATE)

  return(out)
}



#' Compute Quantities of Interest (GATE, GATEcv, URATE) with user defined functions
#' @param user_hte A user-defined function to estimate heterogeneous treatment effects (HTE). The function should take the data as input and return an unit-level continuous score for treatment assignment. We assume those that have score less than 0 should not have treatment. The default is \code{NULL}, which means the heterogeneous treatment effects will be estimated from by the package.
#' @param Tcv A vector of the unit-level binary treatment.
#' @param Ycv A vector of the unit-level continuous outcome.
#' @param data A data frame containing the variables of interest.
#' @param ngates The number of gates to be used in the GATE function.
#' @param ... Additional arguments to be passed to the user-defined function.
#' @importFrom rlang .data
compute_qoi_user <- function(user_hte, Tcv, Ycv, data, ngates, ...) {

  # parameters
  function_name <- as.character(substitute(user_hte))

  # HTE
  tau <- do.call(user_hte, list(data))

  ## GATE
  GATE <- vector("list", length = length(user_hte))
  for (i in seq_along(length(user_hte))) {

    ## Compute GATE
    GATE[[i]] <- GATE(Tcv, tau, Ycv, ngates)

    ## indicate algorithm
    GATE[[i]]$alg <- function_name[i]

    ## indicate group number
    GATE[[i]]$group <- seq_along(GATE[[i]]$gate)
  }

  ## URATE
  URATE <- vector("list", length = length(user_hte))
  for (i in seq_along(length(user_hte))) {

    ## Compute URATE
    URATE[[i]] <- URATE(Tcv, tau, Ycv)

    ## indicate algorithm
    URATE[[i]]$alg <- function_name[i]

  }
  
  out <- list(
        GATE = GATE,
        URATE = URATE)

  return(out)
}
