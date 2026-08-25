#' The Heterogeneity Test for Grouped Average Treatment Effects (GATEs) in Randomized Experiments
#'
#' This function calculates statistics related to the test of heterogeneous treatment effects across groups.
#'
#' The details of the methods for this design are given in Imai and Li (2025).
#' The supplied outcome is centered by default. Optionally, arm-specific
#' outcome regressions can be used to construct a cross-fitted residualized
#' outcome before conducting the test.
#'
#' @importFrom stats cov pchisq
#' @param T A vector of the unit-level binary treatment receipt variable for each sample.
#' @param tau A vector of the unit-level continuous score. Conditional Average Treatment Effect is one possible measure.
#' @param Y A vector of the outcome variable of interest for each sample.
#' @param ngates The number of groups to separate the data into. The groups are determined by \code{tau}. Default is 5.
#' @param center Whether to center the outcome before constructing the test
#'   statistic. Default is `TRUE`.
#' @param residualize Whether to residualize `Y` before testing. Default is
#'   `FALSE`, which uses the supplied outcome without nuisance adjustment.
#' @param X Optional covariate matrix or data frame used for automatic nuisance
#'   cross-fitting.
#' @param p Treatment probability used in residualization. The default is the
#'   empirical treatment proportion. Supply the known assignment probability
#'   when available.
#' @param n_folds Number of nuisance cross-fitting folds. Default is 5.
#' @param nuisance_learner Either `"SuperLearner"` or a custom prediction
#'   function; see [residualize_outcome()].
#' @param SL.library SuperLearner library used for automatic residualization.
#'   The default includes mean, GLM, GAM, glmnet, neural-network, earth, and
#'   pruned-tree learners.
#' @param SL_folds Number of inner SuperLearner folds. Default is 5.
#' @param learner_args Additional arguments for the nuisance learner.
#' @return A list that contains the following items: \item{stat}{The estimated
#' statistic for the test of heterogeneity.} \item{pval}{The p-value of the null
#' hypothesis (that the treatment effects are homogeneous)}
#' @examples
#' T <- rep(c(1, 0), 20)
#' tau <- seq_len(40) / 40
#' Y <- tau * T + rnorm(40)
#' hettestlist <- het.test(T,tau,Y,ngates=5)
#' hettestlist$stat
#' hettestlist$pval
#' @author Michael Lingzhi Li, Technology and Operations Management, Harvard Business School
#' \email{mili@hbs.edu}, \url{https://www.michaellz.com/};
#' @references Imai and Li (2025). \dQuote{Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
#' @keywords evaluation
#' @export het.test
het.test <- function(
    T,
    tau,
    Y,
    ngates = 5,
    center = TRUE,
    residualize = TRUE,
    X = NULL,
    p = mean(T),
    n_folds = 5,
    nuisance_learner = "SuperLearner",
    SL.library = c(
      "SL.mean", "SL.glm", "SL.gam", "SL.glmnet",
      "SL.nnet", "SL.earth", "SL.rpartPrune"
    ),
    SL_folds = 5,
    learner_args = list()
) {
  if ((!is.numeric(T) && !is.logical(T)) || length(T) == 0L ||
      anyNA(T) || !all(T %in% c(0, 1))) {
    stop("T should be binary.")
  }
  if (length(T) != length(tau) || length(tau) != length(Y)) {
    stop("All the data should have the same length.")
  }
  if (length(T) == 0L) {
    stop("The data should have positive length.")
  }
  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop("tau must be a finite numeric vector.")
  }
  if (!is.numeric(Y) || any(!is.finite(Y))) {
    stop("Y must be a finite numeric vector.")
  }
  if (sum(T == 1) < 2L || sum(T == 0) < 2L) {
    stop("Each treatment arm must contain at least two observations.")
  }

  if (isTRUE(residualize)) {
    Y <- residualize_outcome(
      Y = Y,
      T = T,
      X = X,
      p = p,
      n_folds = n_folds,
      nuisance_learner = nuisance_learner,
      SL.library = SL.library,
      SL_folds = SL_folds,
      learner_args = learner_args
    )
  }
  if (isTRUE(center)) {
    Y <- Y - mean(Y)
  }

  T <- as.numeric(T)
  K <- ngates
  n <- length(Y)
  n1 <- sum(T)
  n0 <- n - n1
  fd_label <- .gate_labels(tau, K, "tau")
  papes <- numeric(K)
  kf1s <- numeric(K)
  kf0s <- numeric(K)
  Sfp1s <- vector("list", K)
  Sfp0s <- vector("list", K)

  for (k in seq_len(K)) {
    That <- as.numeric(fd_label == k)
    plim <- 1 / K
    papes[k] <- K * (
      1 / n1 * sum(T * That * Y) +
        1 / n0 * sum(Y * (1 - T) * (1 - That)) -
        plim / n1 * sum(Y * T) -
        (1 - plim) / n0 * sum(Y * (1 - T))
    )
    Sfp1s[[k]] <- ((That - 1 / K) * Y)[T == 1]
    Sfp0s[[k]] <- ((That - 1 / K) * Y)[T == 0]
    kf1s[k] <- mean(Y[T == 1 & That == 1]) -
      mean(Y[T == 0 & That == 1])
    kf0s[k] <- mean(Y[T == 1 & That == 0]) -
      mean(Y[T == 0 & That == 0])
  }

  mcov <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      term1 <- K^2 * (
        stats::cov(Sfp1s[[i]], Sfp1s[[j]]) / n1 +
          stats::cov(Sfp0s[[i]], Sfp0s[[j]]) / n0
      )
      mcov[i, j] <- term1 + .kblock(
        kf1s[i], kf0s[i], kf1s[j], kf0s[j], K, n
      )
    }
  }
  mcov <- diag(diag(mcov), nrow = K, ncol = K)
  if (!is.finite(determinant(mcov)$modulus) || any(diag(mcov) <= 0)) {
    return(list(stat = NA_real_, pval = NA_real_))
  }
  stat <- as.numeric(t(papes) %*% solve(mcov) %*% papes)
  list(stat = stat, pval = stats::pchisq(stat, K, lower.tail = FALSE))
}
