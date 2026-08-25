#' The Heterogeneity Test for Grouped Average Treatment Effects (GATEs) under Cross Validation in Randomized Experiments
#'
#' This function calculates statistics related to the test of heterogeneous treatment effects across groups under cross-validation.
#'
#' The details of the methods for this design are given in Imai and Li (2025).
#' The supplied outcome is centered by default. Optionally, arm-specific
#' outcome regressions can be cross-fitted using the same outer folds before
#' conducting the test.
#'
#' @importFrom stats cov pchisq
#' @param T A vector of the unit-level binary treatment receipt variable for each sample.
#' @param tau A vector of the unit-level continuous score. Conditional Average Treatment Effect is one possible measure.
#' @param Y A vector of the outcome variable of interest for each sample.
#' @param ind A vector of integers (between 1 and number of folds inclusive) indicating which testing set does each sample belong to.
#' @param ngates The number of groups to separate the data into. The groups are determined by \code{tau}. Default is 5.
#' @param center Whether to center the outcome before constructing the test
#'   statistic. Default is `TRUE`.
#' @param residualize Whether to residualize `Y` before testing. Default is
#'   `FALSE`, which uses the supplied outcome without nuisance adjustment.
#' @param X Optional covariate matrix or data frame used for automatic nuisance
#'   cross-fitting. Nuisance models use the folds supplied in `ind`.
#' @param p Treatment probability used in residualization. The default is the
#'   empirical treatment proportion. Supply the known assignment probability
#'   when available.
#' @param nuisance_learner Either `"SuperLearner"` or a custom prediction
#'   function; see [residualize_outcome()].
#' @param SL.library SuperLearner library used for automatic residualization.
#'   The default includes mean, GLM, GAM, glmnet, neural-network, earth, and
#'   pruned-tree learners.
#' @param SL_folds Number of inner SuperLearner folds. Default is 5.
#' @param learner_args Additional arguments for the nuisance learner.
#' @return A list that contains the following items: \item{stat}{The estimated
#' statistic for the test of heterogeneity under cross-validation.} \item{pval}{The p-value of the null
#' hypothesis (that the treatment effects are homogeneous)}
#' @examples
#' T <- rep(c(1, 0), 40)
#' ind <- rep(1:2, each = 40)
#' tau <- cbind(seq_len(80) / 80, seq_len(80) / 80)
#' Y <- tau[, 1] * T + rnorm(80)
#' hettestlist <- hetcv.test(T, tau, Y, ind, ngates = 2)
#' hettestlist$stat
#' hettestlist$pval
#' @author Michael Lingzhi Li, Technology and Operations Management, Harvard Business School
#' \email{mili@hbs.edu}, \url{https://www.michaellz.com/};
#' @references Imai and Li (2025). \dQuote{Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
#' @keywords evaluation
#' @export hetcv.test
hetcv.test <- function(
    T,
    tau,
    Y,
    ind,
    ngates = 5,
    center = TRUE,
    residualize = TRUE,
    X = NULL,
    p = mean(T),
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
  if (!is.matrix(tau) || length(T) != nrow(tau) ||
      nrow(tau) != length(Y)) {
    stop("All the data should have the same length.")
  }
  if (length(T) == 0L) {
    stop("The data should have positive length.")
  }
  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop("tau must be a finite numeric matrix.")
  }
  if (!is.numeric(Y) || any(!is.finite(Y))) {
    stop("Y must be a finite numeric vector.")
  }
  if (length(ind) != length(Y) || anyNA(ind)) {
    stop("ind must have one non-missing fold label per observation.")
  }
  fold_labels <- sort(unique(ind))
  if (!identical(fold_labels, seq_along(fold_labels))) {
    stop("ind must contain consecutive integer labels beginning at one.")
  }
  L <- length(fold_labels)
  if (L < 2L || ncol(tau) < L) {
    stop("tau must contain one score column for each of at least two folds.")
  }

  if (isTRUE(residualize)) {
    if ((!is.matrix(X) && !is.data.frame(X)) || NROW(X) != length(Y)) {
      stop("X must be a matrix or data frame with one row per outcome.")
    }
    if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
      stop("p must be a finite scalar strictly between zero and one.")
    }
    nuisance <- .crossfit_outcome_models(
      Y = Y,
      T = T,
      X = X,
      folds = ind,
      nuisance_learner = nuisance_learner,
      SL.library = SL.library,
      SL_folds = SL_folds,
      learner_args = learner_args
    )
    Y <- as.numeric(Y - (1 - p) * nuisance$q1 - p * nuisance$q0)
  }
  if (isTRUE(center)) {
    Y <- Y - mean(Y)
  }

  T <- as.numeric(T)
  K <- ngates
  papesm <- matrix(NA_real_, L, K)
  cov1s <- array(NA_real_, c(L, K, K))
  cov0s <- array(NA_real_, c(L, K, K))
  kf1t <- matrix(NA_real_, L, K)
  kf0t <- matrix(NA_real_, L, K)
  kf1cv <- matrix(NA_real_, L, K)
  m1s <- m0s <- ms <- numeric(L)
  constant_folds <- integer(0)

  for (fold in seq_len(L)) {
    fold_index <- ind == fold
    Tind <- T[fold_index]
    tauind <- tau[fold_index, fold]
    Yind <- Y[fold_index]
    tauind_full <- tau[, fold]
    m <- length(Yind)
    m1 <- sum(Tind)
    m0 <- m - m1
    if (m1 < 2L || m0 < 2L) {
      stop(sprintf(
        "Fold %d must contain at least two observations in each treatment arm.",
        fold
      ))
    }
    ms[fold] <- m
    m1s[fold] <- m1
    m0s[fold] <- m0
    if (length(unique(tauind)) == 1L) {
      constant_folds <- c(constant_folds, fold)
    }
    fd_label <- .gate_labels(tauind, K, sprintf("fold %d scores", fold))
    Ystar <- matrix(NA_real_, m, K)

    for (gate in seq_len(K)) {
      That <- as.numeric(fd_label == gate)
      plim <- sum(That) / m
      papesm[fold, gate] <- K * (
        1 / m1 * sum(Tind * That * Yind) +
          1 / m0 * sum(Yind * (1 - Tind) * (1 - That)) -
          plim / m1 * sum(Yind * Tind) -
          (1 - plim) / m0 * sum(Yind * (1 - Tind))
      )
      Ystar[, gate] <- (That - 1 / K) * Yind

      if (sum(Tind == 1 & That == 1) > 0 &&
          sum(Tind == 0 & That == 1) > 0) {
        kf1t[fold, gate] <- mean(Yind[Tind == 1 & That == 1]) -
          mean(Yind[Tind == 0 & That == 1])
      }
      if (sum(Tind == 1 & That == 0) > 0 &&
          sum(Tind == 0 & That == 0) > 0) {
        kf0t[fold, gate] <- mean(Yind[Tind == 1 & That == 0]) -
          mean(Yind[Tind == 0 & That == 0])
      }

      if (gate == 1L) {
        upper <- max(tauind[fd_label == gate])
        lower <- -Inf
      } else if (gate == K) {
        upper <- Inf
        lower <- min(tauind[fd_label == gate])
      } else {
        upper <- max(tauind[fd_label == gate])
        lower <- min(tauind[fd_label == gate])
      }
      That_full <- as.numeric(
        tauind_full <= upper & tauind_full >= lower
      )
      if (sum(T == 1 & That_full == 1) > 0 &&
          sum(T == 0 & That_full == 1) > 0) {
        kf1cv[fold, gate] <- mean(Y[T == 1 & That_full == 1]) -
          mean(Y[T == 0 & That_full == 1])
      }
    }
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        cov1s[fold, i, j] <- stats::cov(
          Ystar[Tind == 1, i], Ystar[Tind == 1, j]
        )
        cov0s[fold, i, j] <- stats::cov(
          Ystar[Tind == 0, i], Ystar[Tind == 0, j]
        )
      }
    }
  }

  papes <- colMeans(papesm)
  fold_spread <- stats::cov(papesm)
  guard <- logical(K)
  mcov <- matrix(0, K, K)

  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      term1 <- mean(
        K^2 * (cov1s[, i, j] / m1s + cov0s[, i, j] / m0s),
        na.rm = TRUE
      )
      kappa_block <- mean(
        .kblock(
          kf1t[, i], kf0t[, i], kf1t[, j], kf0t[, j], K, ms
        ),
        na.rm = TRUE
      )
      training_variance <- tryCatch(
        stats::cov(kf1cv[, i], kf1cv[, j], use = "complete.obs"),
        error = function(e) 0
      )
      if (i == j && is.finite(term1 + kappa_block + training_variance) &&
          is.finite(fold_spread[i, i]) &&
          fold_spread[i, i] > term1 + kappa_block + training_variance) {
        guard[i] <- TRUE
      }
      mcov[i, j] <-
        (if (i == j) max(term1 + kappa_block, 0) else
          term1 + kappa_block) / L +
        (if (i == j) max(training_variance, 0) else training_variance)
    }
  }

  if (any(guard)) {
    warning(
      sprintf(
        paste0(
          "Cross-fold spread exceeds the analytic variance for gate(s) %s. ",
          "The score may be unstable across folds."
        ),
        paste(which(guard), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mcov <- diag(diag(mcov), nrow = K, ncol = K)
  if (!is.finite(determinant(mcov)$modulus) || any(diag(mcov) <= 0)) {
    return(list(
      stat = NA_real_, pval = NA_real_, guard = guard,
      constant_folds = constant_folds
    ))
  }
  stat <- as.numeric(t(papes) %*% solve(mcov) %*% papes)
  list(
    stat = stat,
    pval = stats::pchisq(stat, K, lower.tail = FALSE),
    guard = guard,
    constant_folds = constant_folds
  )
}
