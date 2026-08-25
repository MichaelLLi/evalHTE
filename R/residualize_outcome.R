#' Residualize an outcome for heterogeneity testing
#'
#' Constructs
#' \deqn{Y^* = Y - (1-p)Q_1(X) - pQ_0(X),}
#' where \eqn{Q_a(X) = E(Y \mid T=a, X)}. The two outcome regressions are
#' estimated separately by treatment arm and cross-fitted so that every
#' observation receives out-of-fold predictions of both \eqn{Q_0(X)} and
#' \eqn{Q_1(X)}.
#'
#' By default, the outcome regressions are fitted with SuperLearner. A custom
#' learner may instead be supplied as a function with arguments `Y`, `X`, and
#' `newX`. The package calls the function once for each treatment arm in each
#' fold, so it only needs to fit one regression and return one prediction per
#' row of `newX`.
#'
#' @param Y Numeric outcome vector.
#' @param T Binary treatment vector.
#' @param X Matrix or data frame of covariates.
#' @param p Scalar treatment probability. The default is the empirical
#'   treatment proportion. Supply the known assignment probability when it is
#'   available.
#' @param n_folds Number of outer folds used to cross-fit the outcome models.
#'   Treatment-stratified folds are generated internally. Default is 5.
#' @param nuisance_learner Either `"SuperLearner"` or a user-supplied function
#'   with arguments `Y`, `X`, and `newX`. Default is `"SuperLearner"`.
#' @param SL.library Learner library passed to
#'   `SuperLearner::SuperLearner()`. The default includes mean, GLM, GAM,
#'   glmnet, neural-network, earth, and pruned-tree learners.
#' @param SL_folds Number of inner cross-validation folds used by SuperLearner.
#'   This argument is ignored for a custom learner. Default is 5.
#' @param learner_args Named list of additional arguments passed to
#'   `SuperLearner::SuperLearner()` or the custom learner.
#'
#' @return A numeric residualized-outcome vector with the same length as `Y`.
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' T <- rbinom(n, 1, 0.5)
#' Y <- X$x1 + T * (1 + X$x2) + rnorm(n)
#'
#' lm_learner <- function(Y, X, newX) {
#'   training_data <- data.frame(.outcome = Y, X, check.names = FALSE)
#'   fit <- stats::lm(.outcome ~ ., data = training_data)
#'   as.numeric(stats::predict(fit, newdata = as.data.frame(newX)))
#' }
#'
#' Y_res <- residualize_outcome(
#'   Y = Y,
#'   T = T,
#'   X = X,
#'   p = 0.5,
#'   n_folds = 5,
#'   nuisance_learner = lm_learner
#' )
#'
#' @examplesIf requireNamespace("grf", quietly = TRUE)
#' grf_learner <- function(Y, X, newX, num.trees = 2000) {
#'   fit <- grf::regression_forest(as.matrix(X), Y, num.trees = num.trees)
#'   as.numeric(predict(fit, as.matrix(newX))$predictions)
#' }
#'
#' Y_res_grf <- residualize_outcome(
#'   Y = Y,
#'   T = T,
#'   X = X,
#'   p = 0.5,
#'   n_folds = 5,
#'   nuisance_learner = grf_learner,
#'   learner_args = list(num.trees = 500)
#' )
residualize_outcome <- function(
    Y,
    T,
    X,
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
  if (!is.numeric(Y) || length(Y) == 0L || any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (!identical(as.numeric(T), as.numeric(as.logical(T))) ||
      length(T) != length(Y)) {
    stop("T must be a binary vector with the same length as Y.",
         call. = FALSE)
  }
  if ((!is.matrix(X) && !is.data.frame(X)) || NROW(X) != length(Y)) {
    stop("X must be a matrix or data frame with one row per outcome.",
         call. = FALSE)
  }
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be a finite scalar strictly between zero and one.",
         call. = FALSE)
  }
  if (length(n_folds) != 1L || !is.finite(n_folds) ||
      n_folds != as.integer(n_folds) || n_folds < 2L) {
    stop("n_folds must be an integer of at least two.", call. = FALSE)
  }
  if (!is.list(learner_args)) {
    stop("learner_args must be a list.", call. = FALSE)
  }

  n_folds <- as.integer(n_folds)
  if (any(table(factor(T, levels = c(0, 1))) < n_folds)) {
    stop("Each treatment arm must contain at least n_folds observations.",
         call. = FALSE)
  }

  folds <- integer(length(T))
  for (arm in c(0, 1)) {
    indices <- sample(which(T == arm))
    folds[indices] <- rep(seq_len(n_folds), length.out = length(indices))
  }

  nuisance <- .crossfit_outcome_models(
    Y = Y,
    T = T,
    X = X,
    folds = folds,
    nuisance_learner = nuisance_learner,
    SL.library = SL.library,
    SL_folds = SL_folds,
    learner_args = learner_args
  )

  as.numeric(Y - (1 - p) * nuisance$q1 - p * nuisance$q0)
}

.crossfit_outcome_models <- function(
    Y,
    T,
    X,
    folds,
    nuisance_learner,
    SL.library,
    SL_folds,
    learner_args
) {
  fold_labels <- sort(unique(folds))
  q0_hat <- q1_hat <- rep(NA_real_, length(Y))

  for (fold in fold_labels) {
    test_index <- folds == fold
    train_index <- !test_index

    for (arm in 0:1) {
      arm_train <- train_index & T == arm
      n_arm_train <- sum(arm_train)
      if (n_arm_train < 3L) {
        stop("Too few observations to fit an arm-specific outcome model.",
             call. = FALSE)
      }

      if (is.character(nuisance_learner) &&
          length(nuisance_learner) == 1L &&
          identical(tolower(nuisance_learner), "superlearner")) {
        if (!requireNamespace("SuperLearner", quietly = TRUE)) {
          stop("The SuperLearner package is required.", call. = FALSE)
        }
        inner_v <- max(2L, min(as.integer(SL_folds), n_arm_train))
        fit_args <- utils::modifyList(
          list(
            Y = Y[arm_train],
            X = data.frame(X[arm_train, , drop = FALSE]),
            newX = data.frame(X[test_index, , drop = FALSE]),
            family = stats::gaussian(),
            SL.library = SL.library,
            cvControl = list(V = inner_v, shuffle = TRUE),
            env = asNamespace("SuperLearner"),
            verbose = FALSE
          ),
          learner_args
        )
        prediction <- as.numeric(
          do.call(SuperLearner::SuperLearner, fit_args)$SL.predict
        )
      } else if (is.function(nuisance_learner)) {
        prediction <- as.numeric(do.call(
          nuisance_learner,
          utils::modifyList(
            list(
              Y = Y[arm_train],
              X = data.frame(X[arm_train, , drop = FALSE]),
              newX = data.frame(X[test_index, , drop = FALSE])
            ),
            learner_args
          )
        ))
      } else {
        stop(
          "nuisance_learner must be \"SuperLearner\" or a function.",
          call. = FALSE
        )
      }

      if (length(prediction) != sum(test_index) ||
          any(!is.finite(prediction))) {
        stop("The nuisance learner returned invalid predictions.",
             call. = FALSE)
      }
      if (arm == 0L) {
        q0_hat[test_index] <- prediction
      } else {
        q1_hat[test_index] <- prediction
      }
    }
  }

  if (any(!is.finite(q0_hat)) || any(!is.finite(q1_hat))) {
    stop("Cross-fitted outcome predictions are incomplete.", call. = FALSE)
  }
  list(q0 = q0_hat, q1 = q1_hat)
}
