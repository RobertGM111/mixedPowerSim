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
#' @param type Character. One of "raw", "effectSize", or "slopeDifference".
#' @param b0 Numeric. Required for all types. Value of b0.
#' @param b1 Numeric. Required for type = "raw". Value of b1.
#' @param b2 Numeric. Required for type = "raw". Value of b2.
#' @param b3 Numeric. Required for type = "raw". Value of b3.
#' @param b4 Numeric. Required for type = "raw". Value of b4.
#' @param b5 Numeric. Required for type = "raw". Value of b5.
#' @param b6 Numeric. Required for type = "raw". Value of b6.
#' @param b7 Numeric. Required for type = "raw". Value of b7.
#' @param a_dif Numeric. Required for type = "slopeDifference". Value of slope test "a".
#' @param b_dif Numeric. Required for type = "slopeDifference". Value of slope test "b".
#' @param c_dif Numeric. Required for type = "slopeDifference". Value of slope test "c".
#' @param d_dif Numeric. Required for type = "slopeDifference". Value of slope test "d".
#' @param e_dif Numeric. Required for type = "slopeDifference". Value of slope test "e".
#' @param f_dif Numeric. Required for type = "slopeDifference". Value of slope test "f".
#' @param db0 Numeric. Required for type = "effectSize". Effect size of b0.
#' @param db1 Numeric. Required for type = "effectSize". Effect size of b1.
#' @param db2 Numeric. Required for type = "effectSize". Effect size of b2.
#' @param db3 Numeric. Required for type = "effectSize". Effect size of b3.
#' @param db4 Numeric. Required for type = "effectSize". Effect size of b4.
#' @param db5 Numeric. Required for type = "effectSize". Effect size of b5.
#' @param db6 Numeric. Required for type = "effectSize". Effect size of b6.
#' @param db7 Numeric. Required for type = "effectSize". Effect size of b7.
#' @param SDb0 Numeric. Required for all types. Between cluster standard deviation of b0.
#' @param SDresid Numeric. Required for all types. Standard deviation of residuals.
#' @param XWithin Logical. Required for all types. If TRUE, variable X is simulated as a within-cluster variable. If FALSE, variable X is simulated as a between-cluster variable.
#' @param ZWithin Logical. Required for all types. If TRUE, variable Z is simulated as a within-cluster variable. If FALSE, variable Z is simulated as a between-cluster variable.
#' @param WWithin Logical. Required for all types. If TRUE, variable W is simulated as a within-cluster variable. If FALSE, variable W is simulated as a between-cluster variable.
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
simMixedXYZ <- function(nCluster, nObs, type = "raw",
                        b0 = 0, b1 = 0, b2 = 0, b3 = 0,
                        b4 = 0, b5 = 0, b6 = 0, b7 = 0,
                        a_Dif = NA, b_Dif = NA,c_Dif = NA,
                        d_Dif = NA, e_Dif = NA, f_Dif = NA,
                        db0 = NA, db1 = NA, db2 = NA,db3 = NA,
                        db4 = NA, db5 = NA, db6 = NA, db7 = NA,
                        SDb0 = 1, SDresid = 1, XWithin = TRUE,
                        ZWithin = TRUE, WWithin = TRUE){

  if (!type %in% c("raw", "effectSize", "slopeDifference") ){
    stop('Please specify argument "type" as either "raw" or "effectSize" or "slopeDifference".')
  }

  group <- rep(1:nCluster, each = nObs)
  obs <- rep(1:nObs, nCluster)

  if (type == "raw"){
    b0fixed <- b0
    b0rand <- rep(stats::rnorm(nCluster,b0fixed,SDb0),each = nObs)
    b1 <- b1
    b2 <- b2
    b3 <- b3
    b4 <- b4
    b5 <- b5
    b6 <- b6
    b7 <- b7

    db1 <- b1/sqrt(SDb0^2+SDresid^2)
    db2 <- b2/sqrt(SDb0^2+SDresid^2)
    db3 <- b3/sqrt(SDb0^2+SDresid^2)
    db4 <- b4/sqrt(SDb0^2+SDresid^2)
    db5 <- b5/sqrt(SDb0^2+SDresid^2)
    db6 <- b6/sqrt(SDb0^2+SDresid^2)
    db7 <- b7/sqrt(SDb0^2+SDresid^2)

    a_Dif <- 2*b5 + 2*b7
    b_Dif <- 2*b4 + 2*b7
    c_Dif <- 2*b4 - 2*b7
    d_Dif <- 2*b5 - 2*b7
    e_Dif <- 2*b4 + 2*b5
    f_Dif <- 2*b4 - 2*b5

  }

  if (type == "slopeDifference"){
    Xmat <- matrix(c(0,2,2,
                     2,0,2,
                     2,0,-2,
                     0,2,-2,
                     2,2,0,
                     2,-2,0), nrow = 6, ncol = 3, byrow = TRUE)
    Ymat <- c(a_Dif, b_Dif, c_Dif, d_Dif, e_Dif, f_Dif)
    b <- solve(t(Xmat)%*%Xmat)%*%t(Xmat)%*%Ymat

    b0fixed <- b0
    b0rand <- rep(stats::rnorm(nCluster,b0fixed,SDb0),each = nObs)
    b1 <- b1
    b2 <- b2
    b3 <- b3
    b4 <- b[1]
    b5 <- b[2]
    b6 <- b6
    b7 <- b[3]

    db1 <- b1/sqrt(SDb0^2+SDresid^2)
    db2 <- b2/sqrt(SDb0^2+SDresid^2)
    db3 <- b3/sqrt(SDb0^2+SDresid^2)
    db4 <- b4/sqrt(SDb0^2+SDresid^2)
    db5 <- b5/sqrt(SDb0^2+SDresid^2)
    db6 <- b6/sqrt(SDb0^2+SDresid^2)
    db7 <- b7/sqrt(SDb0^2+SDresid^2)
  }

  if (type == "effectSize"){
    b0fixed <- db0*sqrt(SDb0^2+SDresid^2)
    b0rand <- rep(stats::rnorm(nCluster,b0fixed,SDb0),each = nObs)
    b1 <- db1*sqrt(SDb0^2+SDresid^2)
    b2 <- db2*sqrt(SDb0^2+SDresid^2)
    b3 <- db3*sqrt(SDb0^2+SDresid^2)
    b4 <- db4*sqrt(SDb0^2+SDresid^2)
    b5 <- db5*sqrt(SDb0^2+SDresid^2)
    b6 <- db6*sqrt(SDb0^2+SDresid^2)
    b7 <- db7*sqrt(SDb0^2+SDresid^2)

    a_Dif <- 2*b5 + 2*b7
    b_Dif <- 2*b4 + 2*b7
    c_Dif <- 2*b4 - 2*b7
    d_Dif <- 2*b5 - 2*b7
    e_Dif <- 2*b4 + 2*b5
    f_Dif <- 2*b4 - 2*b5
  }


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

  coefMat <- data.frame("coefficients" = paste("b", 0:7, sep = ""),
                        "values" = c(b0fixed, b1, b2, b3, b4, b5, b6, b7),
                        "d" = c(NA, db1, db2, db3, db4, db5, db6, db7))

  slopeDifferences <- data.frame("test" = c("a", "b", "c", "d", "e", "f"),
                                 "values" = c(a_Dif, b_Dif, c_Dif, d_Dif, e_Dif, f_Dif))

  simDat <- data.frame("Cluster" = group, "Observation" = obs,Y,X,W,Z)
  resList <- list("simulatedData" = simDat, "coefficients" = coefMat, "slopeDifferences" = slopeDifferences)
  return(resList)
}




