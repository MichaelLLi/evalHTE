#' Estimation of the Grouped Average Treatment Effects (GATEs) in Randomized Experiments
#'
#' Estimates grouped average treatment effects from a continuous score.
#'
#' @param D Binary treatment indicator (0 or 1).
#' @param tau Continuous score vector.
#' @param Y Outcome vector.
#' @param ngates Number of groups (at least 2).
#' @param centered Whether to center outcomes before estimation.
#' @return A list with group estimates \code{gate} and standard errors \code{sd},
#' ordered by increasing score.
#' @details
#' Uses \code{\link[evalITR]{GATE}} for inference.
#' @examples
#' D = c(1,0,1,0,1,0,1,0)
#' tau = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#' Y = c(4,5,0,2,4,1,-4,3)
#' gatelist <- GATE(D,tau,Y,ngates=2)
#' gatelist$gate
#' gatelist$sd
#' @author Michael Lingzhi Li, Technology and Operations Management, Harvard Business School
#' \email{mili@hbs.edu}, \url{https://www.michaellz.com/};
#' @references Imai and Li (2022). \dQuote{Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
#' @keywords evaluation
#' @export GATE
GATE <- function(D, tau, Y, ngates = 5, centered = FALSE) {
  evalITR::GATE(T = D, tau = tau, Y = Y, ngates = ngates,
               centered = centered)
}
