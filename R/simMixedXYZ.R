#' Simulate data for a three-way interaction of nested observations.
#'
#' Creates a simulated data set based upon the supplied arguments.
#'
#' By specifying effect sizes (\code{db0}, ..., \code{db7}) of betas as seen in Dawson and Richter (2006),
#' number of clusters, and number of observations within cluster, this function creates a data set with the
#' appropriate relationships. Variables X, Y, and Z are each drawn from N~(0,1).
#'
#' @usage simMixedXYZ(nCluster, nObs,
#'  db0, db1, db2, db3, db4, db5, db6, db7,
#'  SDb0, SDresid, XWithin = TRUE, ZWithin = TRUE, WWithin = TRUE)
#' @param nCluster Integer. Number of clusters (level-2 observations).
#' @param nObs Integer. Number of observations within each cluster.
#' @param db0 Numeric. Effect size of b0.
#' @param db1 Numeric. Effect size of b1.
#' @param db2 Numeric. Effect size of b2.
#' @param db3 Numeric. Effect size of b3.
#' @param db4 Numeric. Effect size of b4.
#' @param db5 Numeric. Effect size of b5.
#' @param db6 Numeric. Effect size of b6.
#' @param db7 Numeric. Effect size of b7.
#' @param SDb0 Numeric. Between cluster standard deviation of b0.
#' @param SDresid Numeric. Standard deviation of residuals.
#' @param XWithin Logical. If TRUE, variable X is simulated as a within-cluster variable. If FALSE, variable X is simulated as a between-cluster variable.
#' @param ZWithin Logical. If TRUE, variable Z is simulated as a within-cluster variable. If FALSE, variable Z is simulated as a between-cluster variable.
#' @param WWithin Logical. If TRUE, variable W is simulated as a within-cluster variable. If FALSE, variable W is simulated as a between-cluster variable.
#' @return A data frame with the following columns:
#' \itemize{
#' \item{Cluster:}{ Cluster identification number.}
#' \item{Observation:}{ Observation number within cluster.}
#' \item{Y:}{ Outcome variable.}
#' \item{X:}{ X variable.}
#' \item{Z:}{ Z variable.}
#' \item{W:}{ W variable.}
#' }
#' @examples
#' library(mixedPowerSim)
#' nCluster <- 20
#' nObs <- 15
#' db0 <- .1
#' db1 <- .1
#' db2 <- .1
#' db3 <- .1
#' db4 <- .1
#' db5 <- .1
#' db6 <- .1
#' db7 <- .1
#' SDb0 <- 1
#' SDresid <- 5
#' simDat <- simMixedXYZ(nCluster = nCluster, nObs = nObs,
#'  db0 = db0, db1 = db1, db2 = db2,
#'  db3 = db3, db4 = db4, db5 = db5, db6 = db6, db7 = db7,
#'  SDb0 = SDb0, SDresid = SDresid, XWithin = TRUE, ZWithin = FALSE, WWithin = TRUE)
#' @export
#'
simMixedXYZ <- function(nCluster, nObs, db0, db1, db2,
                        db3, db4, db5, db6, db7,
                        SDb0, SDresid, XWithin = TRUE,
                        ZWithin = TRUE, WWithin = TRUE){

  group <- rep(1:nCluster, each = nObs)
  obs <- rep(1:nObs, nCluster)

  b0fixed <- db0*SDb0*SDresid
  b0rand <- rep(stats::rnorm(nCluster,b0fixed,SDb0),each = nObs)
  b1 <- db1*sqrt(SDb0^2+SDresid^2)
  b2 <- db2*sqrt(SDb0^2+SDresid^2)
  b3 <- db3*sqrt(SDb0^2+SDresid^2)
  b4 <- db4*sqrt(SDb0^2+SDresid^2)
  b5 <- db5*sqrt(SDb0^2+SDresid^2)
  b6 <- db6*sqrt(SDb0^2+SDresid^2)
  b7 <- db7*sqrt(SDb0^2+SDresid^2)

  if(XWithin == TRUE){
    X <- stats::rnorm(nCluster*nObs)
  } else {
    X <- rep(stats::rnorm(nCluster), each = nObs)
  }
  if(WWithin == TRUE){
    W <- stats::rnorm(nCluster*nObs)
  } else {
    W <- rep(stats::rnorm(nCluster), each = nObs)
  }
  if(ZWithin == TRUE){
    Z <- stats::rnorm(nCluster*nObs)
  } else {
    Z <- rep(stats::rnorm(nCluster), each = nObs)
  }

  Y <- b0rand + b1*X + b2*W + b3*Z + b4*X*Z + b5*X*W + b6*Z*W + b7*X*Z*W + stats::rnorm(nCluster*nObs,0,SDresid)

  simDat <- data.frame("Cluster" = group, "Observation" = obs,Y,X,W,Z)
  return(simDat)
}
