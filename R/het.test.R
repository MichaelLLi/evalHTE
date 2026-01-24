#' The Heterogeneity Test for Grouped Average Treatment Effects (GATEs) in Randomized Experiments
#'
#' This function calculates statistics related to the test of heterogeneous treatment effects across groups.
#'
#' The details of the methods for this design are given in Imai and Li (2022).
#'
#'
#' @importFrom stats cov pchisq
#' @param D A vector of the unit-level binary treatment receipt variable for each sample.
#' @param tau A vector of the unit-level continuous score. Conditional Average Treatment Effect is one possible measure.
#' @param Y A vector of the outcome variable of interest for each sample.
#' @param ngates The number of groups to separate the data into. The groups are determined by \code{tau}. Default is 5.
#' @return A list that contains the following items: \item{stat}{The estimated
#' statistic for the test of heterogeneity.} \item{pval}{The p-value of the null
#' hypothesis (that the treatment effects are homogeneous)}
#' @examples
#' D = c(1,0,1,0,1,0,1,0)
#' tau = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#' Y = c(4,5,0,2,4,1,-4,3)
#' hettestlist <- het.test(D,tau,Y,ngates=5)
#' hettestlist$stat
#' hettestlist$pval
#' @author Michael Lingzhi Li, Technology and Operations Management, Harvard Business School
#' \email{mili@hbs.edu}, \url{https://www.michaellz.com/};
#' @references Imai and Li (2022). \dQuote{Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
#' @keywords evaluation
#' @export het.test
het.test <- function(D, tau, Y, ngates = 5) {
  if (!(identical(as.numeric(D),as.numeric(as.logical(D))))) {
    stop("D should be binary.")
  }
  if ((length(D)!=length(tau)) | (length(tau)!=length(Y))) {
    stop("All the data should have the same length.")
  }
  if (length(D)==0) {
    stop("The data should have positive length.")
  }
  n = length(Y)
  n1 = sum(D)
  n0 = n-n1
  fd_label = ntile(tau, ngates)
  vargts = numeric(ngates)
  papes = numeric(ngates)
  kf1s = numeric(ngates)
  kf0s = numeric(ngates)
  Sfp1s = as.list(numeric(ngates))
  Sfp0s = as.list(numeric(ngates))
  mcov = matrix(0, nrow = ngates, ncol = ngates)
  for (i in 1:ngates) {
    That = as.numeric(fd_label == i)
    plim = 1 / ngates
    papes[i] = ngates * (1/n1*sum(D*That*Y)+1/n0*sum(Y*(1-D)*(1-That))-plim/n1*sum(Y*D)-(1-plim)/n0*sum(Y*(1-D)))
    Sfp1s[[i]] = (That*Y)[D==1]
    Sfp0s[[i]] = (That*Y)[D==0]
    kf1s[i] = mean(Y[D==1 & That==1])-mean(Y[D==0 & That==1])
    kf0s[i] = mean(Y[D==1 & That==0])-mean(Y[D==0 & That==0])
  }
  for (i in 1:ngates) {
    for (j in 1:ngates) {
      mcov[i,j] = ngates ^ 2 * (cov(Sfp1s[[i]],Sfp1s[[j]]) / n1 + cov(Sfp0s[[i]],Sfp0s[[j]]) / n0) +
        1/ (ngates * (n - 1)) *((ngates - 1)*(kf1s[i]^2-kf1s[i]*kf0s[i]+kf1s[j]^2-kf1s[j]*kf0s[j]) - ngates * (ngates - 1) *kf1s[i] * kf1s[j])
    }
  }
  mcov = diag(diag(mcov), nrow = ngates, ncol = ngates)
  if (is.finite(determinant(mcov)$modulus)) {
    stat = t(papes) %*% solve(mcov) %*% papes
    return(list(stat=stat,pval=pchisq(stat, ngates, lower.tail = FALSE)))
  } else {
    return(list(stat=NA,pval=NA))
  }
}
