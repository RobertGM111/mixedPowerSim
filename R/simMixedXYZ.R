#' Simulate data for a three-way interaction of nested observations.
#'
#' `simMixedXYZ()` creates a simulated data set based upon the supplied arguments.
#'
#' By specifying effect sizes (\code{dtb0}, ..., \code{dtb7}) of betas as seen in Dawson and Richter (2006),
#' number of clusters, and number of observations within cluster, this function creates a data set with the
#' appropriate relationships. Variables X, Y, and Z are each drawn from N~(0,1).
#'
#' @param nCluster Integer. Number of clusters (level-2 observations).
#' @param nObs Integer. Number of observations within each cluster.
#' @param dtb0 Numeric. Effect size of b0.
#' @param dtb1 Numeric. Effect size of b1.
#' @param dtb2 Numeric. Effect size of b2.
#' @param dtb3 Numeric. Effect size of b3.
#' @param dtb4 Numeric. Effect size of b4.
#' @param dtb5 Numeric. Effect size of b5.
#' @param dtb6 Numeric. Effect size of b6.
#' @param SDb0 Numeric. Between cluster standard deviation of b0.
#' @param SDresid Numeric. Standard deviation of residuals.
#' @param XWithin Logical. If TRUE, variable X is simulated as a within-cluster variable. If False, variable X is simulated as a between-cluster variable.
#' @param ZWithin Logical. If TRUE, variable Z is simulated as a within-cluster variable. If False, variable Z is simulated as a between-cluster variable.
#' @param WWithin Logical. If TRUE, variable W is simulated as a within-cluster variable. If False, variable W is simulated as a between-cluster variable.
#' @return A data frame with the following columns:
#' \itemize{
#' \item{"Cluster"}{Cluster identification number.}
#' \item{"Observation"}{Obervcation number within cluster.}
#' \item{"Y"}{Outcome variable.}
#' \item{"X"}{X variable.}
#' \item{"Z"}{Z variable.}
#' \item{"W"}{W variable.}
#' }
#' @examples
#' nGroups <- 20
#' nObs <- 15
#' dtb0 <- .1
#' dtb1 <- .1
#' dtb2 <- .1
#' dtb3 <- .1
#' dtb4 <- .1
#' dtb5 <- .1
#' dtb6 <- .1
#' dtb7 <- .1
#' SDb0 <- 1
#' SDresid <- 5
#' simDat <- simMixedXYZ(nGroups = nGroups, nObs = nObs, dtb0 = dtb0, dtb1 = dtb1, dtb2 = dtb2,
#'                       dtb3 = dtb3, dtb4 = dtb4, dtb5 = dtb5, dtb6 = dtb6, dtb7 = dtb7,
#'                      SDb0 = SDb0, SDresid = SDresid, XWithin = TRUE, ZWithin = FALSE, WWithin = TRUE)

simMixedXYZ <- function(nCluster, nObs, dtb0, dtb1, dtb2, dtb3, dtb4, dtb5, dtb6, dtb7, SDb0, SDresid, XWithin = TRUE,
                        ZWithin = TRUE, WWithin = TRUE){
  group <- rep(1:nCluster, each = nObs)
  obs <- rep(1:nObs, nCluster)

  b0fixed <- dtb0*SDb0*SDresid
  b0rand <- rep(rnorm(nCluster,b0fixed,SDb0),each = nObs)
  b1 <- dtb1*sqrt(SDb0^2+SDresid^2)
  b2 <- dtb2*sqrt(SDb0^2+SDresid^2)
  b3 <- dtb3*sqrt(SDb0^2+SDresid^2)
  b4 <- dtb4*sqrt(SDb0^2+SDresid^2)
  b5 <- dtb5*sqrt(SDb0^2+SDresid^2)
  b6 <- dtb6*sqrt(SDb0^2+SDresid^2)
  b7 <- dtb7*sqrt(SDb0^2+SDresid^2)

  if(XWithin == TRUE){
    X <- rnorm(nCluster*nObs)
  } else {
    X <- rep(rnorm(nCluster), each = nObs)
  }
  if(WWithin == TRUE){
    W <- rnorm(nCluster*nObs)
  } else {
    W <- rep(rnorm(nCluster), each = nObs)
  }
  if(ZWithin == TRUE){
    Z <- rnorm(nCluster*nObs)
  } else {
    Z <- rep(rnorm(nCluster), each = nObs)
  }

  Y <- b0rand + b1*X + b2*W + b3*Z + b4*X*Z + b5*X*W + b6*Z*W + b7*X*Z*W + rnorm(nCluster*nObs,0,SDresid)

  simDat <- data.frame("Cluster" = group, "Observation" = obs,Y,X,W,Z)
  return(simDat)
}
