#' Run a power analysis of slope differences.
#'
#' Creates a simulated data set based upon the supplied arguments.
#'
#' By specifying effect sizes (\code{db0}, ..., \code{db7}) of betas as seen in Dawson and Richter (2006),
#' number of clusters, and number of observations within cluster, this function creates a data set with the
#' appropriate relationships. Variables X, Y, and Z are each drawn from N~(0,1).
#'
#' @usage XYZpowerSim(R = 1000, clusters = NA, obs =  NA, sigLevel = .05 ...)
#' @param R Integer. Number of replications within each cell of simulation.
#' @param clusters A vector of integers. Number of clusters to be crossed with \code{obs} e.g., \code{c(5,10,15)}.
#' @param obs A vector of integers. Number of observations for each level of to be crossed with \code{clusters} e.g., \code{c(5,10,15)}.
#' @param sigLevel Numeric. $alpha$ value for determine significant p-values in each simulation.
#' @param ... Arguments passed to \link{simMixedXYZ}.
#' @return A data frame of simulated parameters.
#' @examples
#' library(mixedPowerSim)
#' nCluster <- c(5, 10, 15)
#' nObs <- c(5, 10, 15)
#' reliability <- .8
#' db0 <- .1
#' SDb0 <- 1
#' SDresid <- 5
#' XYZpowerSim(R = 3, clusters = nCluster, obs =  nObs, sigLevel = .05, type = "slopeDifference",
#' b0 = db0, a_Dif = 0, b_Dif = 0,c_Dif = 0,
#' d_Dif = 0, e_Dif = 0, f_Dif = 0,
#' SDb0 = SDb0, SDresid = SDresid, XWithin = TRUE, ZWithin = FALSE, WWithin = TRUE, reliability = reliability)
#' @export
#'

XYZpowerSim <- function(R = 1000, clusters = NA, obs =  NA, sigLevel = .05, ...){
  clustObs <- expand.grid(clusters, obs)
  simRes <- rep(NA,R)
  resMat <- data.frame(cbind(clustObs,NA,NA,NA,NA,NA,NA))
  colnames(resMat) <- c("numClusters", "numObs", "power a", "power b", "power c", "power d", "power e", "power f")

  pbCount <- 0
  pb <- utils::txtProgressBar(min = pbCount, max = length(clusters)*length(obs)*R, style = 3)

  for (i in 1:length(clusters)){
    for (j in 1:length(obs)){
      aSig <- 0
      bSig <- 0
      cSig <- 0
      dSig <- 0
      eSig <- 0
      fSig <- 0
      for (k in 1:R){

        simDat <- simMixedXYZ(nCluster = clusters[i], nObs = obs[j], ...)
        tempDif <- suppressWarnings(slopeDif(lme4::lmer(Y~X*Z*W + (1|Cluster), data = simDat$simulatedData)))

        aSig <- ifelse(tempDif$SlopeTests[1,6] < sigLevel, aSig + 1, aSig)
        bSig <- ifelse(tempDif$SlopeTests[2,6] < sigLevel, bSig + 1, bSig)
        cSig <- ifelse(tempDif$SlopeTests[3,6] < sigLevel, cSig + 1, cSig)
        dSig <- ifelse(tempDif$SlopeTests[4,6] < sigLevel, dSig + 1, dSig)
        eSig <- ifelse(tempDif$SlopeTests[5,6] < sigLevel, eSig + 1, eSig)
        fSig <- ifelse(tempDif$SlopeTests[6,6] < sigLevel, fSig + 1, fSig)

        pbCount <- pbCount + 1
        setTxtProgressBar(pb, pbCount)
      }
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 3] <- aSig/R
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 4] <- bSig/R
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 5] <- cSig/R
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 6] <- dSig/R
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 7] <- eSig/R
      resMat[which(resMat$numClusters == clusters[i] & resMat$numObs == obs[j]), 8] <- fSig/R
    }
  }

  close(pb)

  return(resMat)
}




