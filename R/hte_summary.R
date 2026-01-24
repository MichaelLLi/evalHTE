#' Summarize Heterogeneity and Consistency Tests
#' @param object An object of \code{evaluate_hte} class (typically an output of \code{evaluate_hte()} function).
#' @param ... Other parameters.
#' @import purrr
#' @importFrom stats pnorm
#' @return An object of class \code{summary.hte}, which is a list containing:
#'   \describe{
#'     \item{GATE}{A tibble with group average treatment effect estimates, including columns: group, algorithm, estimate, std.deviation, lower, upper, z.score, and p.value.}
#'     \item{URATE}{A tibble with uplift rate estimates for exceptional responders, including columns: algorithm, estimate, std.deviation, conf.low.uniform, z.score, and p.value. Returns NULL when cross-validation is used.}
#'   }
#' @export
summary.hte <- function(object, ...) {

  # parameters
  out <- urate_algs_vec <- urate_user_vec <- list()

  estimate_algs = object$out_algs
  estimate_user = object$out_user

  # algorithm   <- object$df$algorithm
  # cv          <- object$cv
  # fit         <- object$qoi

# -----------------------------------------
# estimate HTE from ML algorithms
# -----------------------------------------

if(length(estimate_algs) != 0){

  # parameters
  algorithm   <- estimate_algs$df$algorithm
  cv          <- estimate_algs$cv
  fit         <- estimate_algs$qoi

  # compute quantities under sample splitting -----------------------------------------

  if (cv == FALSE) {

    # group HTE
    gate_algs_vec <- fit$GATE %>%
      purrr::map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(
        statistic = gate / sd,
        p.value = 2 * pnorm(abs(gate / sd), lower.tail = FALSE),
        upper = gate + qnorm(0.975) * sd,
        lower = gate - qnorm(0.975) * sd
      ) %>%
      rename(
        estimate = gate,
        std.deviation = sd,
        algorithm = alg,
        group = group,
        z.score = statistic
      ) %>%
      select(
        group, algorithm, estimate, std.deviation, lower, upper, z.score, p.value)

    urate_algs_vec <- fit$URATE %>%
      purrr::map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(
        statistic = rate / sd,
        p.value = 2 * pnorm(abs(rate / sd), lower.tail = FALSE)
      ) %>%
      group_by(alg) %>%
      mutate(
        fraction = seq(1,length(rate))/length(rate),
        est = rate - 1.2*sd - 0.68*sd[length(sd)] * length(sd)/seq(1, length(sd)),
        best_ind = which.max(est),
        # estimate
        proportion = fraction[best_ind],
        best_rate = rate[best_ind],
        conf.low.uniform = rate[best_ind] - 1.2*sd[best_ind] - 0.68* sd[length(sd)] * length(sd)/best_ind
      ) %>%
      filter(rate == best_rate) %>%
      rename(
        estimate = rate,
        std.deviation = sd,
        algorithm = alg,
        z.score = statistic
      ) %>%
      select(
        algorithm, estimate, std.deviation, conf.low.uniform, z.score, p.value)

    out <- list(
      GATE = gate_algs_vec,
      URATE = urate_algs_vec)
  }

  # compute quantities under cross-validation -----------------------------------------

  if (cv == TRUE) {

    # group HTE
    gate_algs_vec <- fit$GATE %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(
        statistic = gate / sd,
        p.value = 2 * pnorm(abs(gate / sd), lower.tail = FALSE),
        upper = gate + qnorm(0.975) * sd,
        lower = gate - qnorm(0.975) * sd
      ) %>%
      rename(
        estimate = gate,
        std.deviation = sd,
        algorithm = alg,
        group = group,
        z.score = statistic
      ) %>%
      select(
        group, algorithm, estimate, std.deviation, lower, upper, z.score, p.value)

    # exceptional reponders not supported for CV
    urate_algs_vec <- NULL

  }

  out <- list(
    GATE = gate_algs_vec,
    URATE = urate_algs_vec)

}

if(length(estimate_user) != 0){

  # parameters
  algorithm   <- estimate_user$df$algorithm
  cv          <- estimate_user$cv
  fit         <- estimate_user$qoi

    # group HTE
    gate_user_vec <- fit$GATE %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(
        statistic = gate / sd,
        p.value = 2 * pnorm(abs(gate / sd), lower.tail = FALSE),
        upper = gate + qnorm(0.975) * sd,
        lower = gate - qnorm(0.975) * sd
      ) %>%
      rename(
        estimate = gate,
        std.deviation = sd,
        algorithm = alg,
        group = group,
        z.score = statistic
      ) %>%
      select(
        group, algorithm, estimate, std.deviation, lower, upper, z.score, p.value)

    # exceptional reponders
    urate_user_vec <- fit$URATE %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(
        statistic = rate / sd,
        p.value = 2 * pnorm(abs(rate / sd), lower.tail = FALSE)
      ) %>%
      group_by(alg) %>%
      mutate(
        fraction = seq(1,length(rate))/length(rate),
        est = rate - 1.2*sd - 0.68*sd[length(sd)] * length(sd)/seq(1, length(sd)),
        best_ind = which.max(est),
        # estimate
        proportion = fraction[best_ind],
        best_rate = rate[best_ind],
        conf.low.uniform = rate[best_ind] - 1.2*sd[best_ind] - 0.68* sd[length(sd)] * length(sd)/best_ind
      ) %>%
      filter(rate == best_rate) %>%
      rename(
        estimate = rate,
        std.deviation = sd,
        algorithm = alg,
        z.score = statistic
      ) %>%
      select(
        algorithm, estimate, std.deviation, conf.low.uniform, z.score, p.value)

  out <- list(
    GATE = bind_rows(gate_algs_vec, gate_user_vec),
    URATE = bind_rows(urate_algs_vec, urate_user_vec)
  )

}


  class(out) <- c("summary.hte", class(out))

  return(out)
}


#' Print
#' @importFrom cli cat_rule
#' @param x An object of \code{summary.hte} class. This is typically an output of \code{summary.hte()} function.
#' @param ... Other parameters. Currently not supported.
#' @return No return value, called for side effects (prints summary tables to console).
#' @export
print.summary.hte <- function(x, ...) {
  # GATE
  cli::cat_rule(left = "GATE")
  print(as.data.frame(x[["GATE"]]), digits = 2)
  cli::cat_line("")

  # URATE
  cli::cat_rule(left = "URATE")
  if (is.null(x[["URATE"]]) || ncol(x[["URATE"]]) == 0) {
    cli::cat_line("Not supported with cross-validation")
  } else {
    print(as.data.frame(x[["URATE"]]), digits = 2)
  }
  cli::cat_line("")

}


#' Summarize Heterogeneity and Consistency Tests
#' @param object An object of \code{test_hte} class (typically an output of \code{test_hte()} function).
#' @param ... Other parameters.
#' @importFrom stats pnorm
#' @return An object of class \code{summary.test_hte}, which is a list containing:
#'   \describe{
#'     \item{Consistency}{A tibble with consistency test results, including columns: algorithm, statistic, and p.value (for sample splitting).}
#'     \item{Heterogeneity}{A tibble with heterogeneity test results, including columns: algorithm, statistic, and p.value (for sample splitting).}
#'     \item{Consistency_cv}{A tibble with consistency test results for cross-validation.}
#'     \item{Heterogeneity_cv}{A tibble with heterogeneity test results for cross-validation.}
#'   }
#'   Note: The output contains either the first two or last two elements depending on whether cross-validation was used.
#' @export
summary.test_hte <- function(object, ...) {
  out            <- list()
  consist_tibble <- tibble()
  het_tibble     <- tibble()

  ## -----------------------------------------
  ## hypothesis tests
  ## -----------------------------------------
  if (names(object[1]) == "consist") {

    # parameters for test_hte object
    consist        <- object$consist
    het            <- object$het
    consist_names <- names(consist)
    het_names <- names(het)

    # reformat
    out[["Consistency"]] <- consist %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(algorithm = consist_names) %>%
      rename(statistic = stat,
            p.value = pval) %>%
      select(algorithm, statistic, p.value)


    out[["Heterogeneity"]] <- het %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(algorithm = het_names) %>%
      rename(statistic = stat,
            p.value = pval) %>%
      select(algorithm, statistic, p.value)
  }


  if (names(object[1]) == "consistcv") {

    # parameters for test_hte object
    consist <- object$consistcv
    het <- object$hetcv
    consist_names <- names(consist)
    het_names <- names(het)

    # reformat
    out[["Consistency_cv"]] <- consist %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(algorithm = consist_names) %>%
      rename(statistic = stat,
            p.value = pval) %>%
      select(algorithm, statistic, p.value)

    out[["Heterogeneity_cv"]] <- het %>%
      map(., ~ as_tibble(.)) %>%
      bind_rows() %>%
      mutate(algorithm = het_names) %>%
      rename(statistic = stat,
            p.value = pval) %>%
      select(algorithm, statistic, p.value)
  }

  class(out) <- c("summary.test_hte", class(out))

  return(out)
}

#' Print
#' @importFrom cli cat_rule
#' @param x An object of \code{summary.test_hte} class. This is typically an output of \code{summary.test_hte()} function.
#' @param ... Other parameters.
#' @return No return value, called for side effects (prints test results to console).
#' @export
print.summary.test_hte <- function(x, ...) {

  # Rank Consistency Test
  cli::cat_rule(left = "Rank Consistency Test Results")
  print(as.data.frame(x[["Consistency"]], digits = 2))
  cli::cat_line("")

  # Group Heterogeneity Test
  cli::cat_rule(left = "Group Heterogeneity Test Results")
  print(as.data.frame(x[["Heterogeneity"]], digits = 2))
  cli::cat_line("")
}

