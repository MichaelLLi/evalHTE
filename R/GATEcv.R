#' Estimation of the Grouped Average Treatment Effects (GATEs) in Randomized Experiments Under Cross Validation
#'
#' Estimates grouped average treatment effects from cross-validated scores.
#'
#' @param D Binary treatment indicator (0 or 1).
#' @param tau A matrix of scores with one column per fold. Column \code{i}
#' contains predictions for all observations from a model trained without fold \code{i}.
#' @param Y Outcome vector.
#' @param ind Integer validation-fold labels starting at 1.
#' @param ngates Number of groups (at least 2).
#' @param centered Whether to center outcomes before estimation.
#' @return A list with group estimates \code{gate} and standard errors \code{sd},
#' ordered by increasing score.
#' @details
#' Uses \code{\link[evalITR]{GATEcv}} for inference.
#' @examples
#' D = c(1,0,1,0,1,0,1,0)
#' tau = matrix(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9),nrow = 8, ncol = 2)
#' Y = c(4,5,0,2,4,1,-4,3)
#' ind = c(rep(1,4),rep(2,4))
#' gatelist <- GATEcv(D, tau, Y, ind, ngates = 2)
#' gatelist$gate
#' gatelist$sd
#' @author Michael Lingzhi Li, Technology and Operations Management, Harvard Business School
#' \email{mili@hbs.edu}, \url{https://www.michaellz.com/};
#' @references Imai and Li (2022). \dQuote{Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
#' @keywords evaluation
#' @export GATEcv
#'
#'
GATEcv <- function(D, tau, Y, ind, ngates = 5, centered = FALSE) {
  evalITR::GATEcv(T = D, tau = tau, Y = Y, ind = ind, ngates = ngates,
    centered = centered)
}
