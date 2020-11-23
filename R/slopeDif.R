#' Conduct a slope difference test for three-way interaction mixed-effects models.
#'
#' Calculates the expected value of the outcome variable (Y) at all combinations of high and low values (+/- 1SD) of X, Z, and W.
#' Additionally, this function calculates the slope differences between all slope pairs.
#'
#'
#' @usage slopeDif(x)
#' @param x An object of class "lmerMod". The result of an \code{lmer()} function with a random intercept.
#' @return An object of class "slopetest" with the following components:
#' \itemize{
#' \item{Expected_Y: }{A data frame with the expected values of Y at high and low values of X, Z, and W.}
#' \item{SlopeTests: }{A data frame representing slope tests as defined in Dawson and Richter (2006).}
#' }
#' @examples
#' library(mixedPowerSim)
#' library(lme4)
#'
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
#'
#' simDat <- simMixedXYZ(nCluster = nCluster, nObs = nObs,
#'  db0 = db0, db1 = db1, db2 = db2,
#'  db3 = db3, db4 = db4, db5 = db5, db6 = db6, db7 = db7,
#'  SDb0 = SDb0, SDresid = SDresid, XWithin = TRUE, ZWithin = FALSE, WWithin = TRUE)
#'
#'  mod <- lmer(Y~X * Z * W + (1|Cluster), data = simDat)
#'
#'  slopeDif(mod)
#' @export
slopeDif <- function(x){
  b0 <- lme4::fixef(x)[1]
  b1 <- lme4::fixef(x)[2]
  b2 <- lme4::fixef(x)[3]
  b3 <- lme4::fixef(x)[4]
  b4 <- lme4::fixef(x)[5]
  b5 <- lme4::fixef(x)[6]
  b6 <- lme4::fixef(x)[7]
  b7 <- lme4::fixef(x)[8]

  Xhi <- mean(x@frame[,2]) + stats::sd(x@frame[,2])
  Xlo <- mean(x@frame[,2]) - stats::sd(x@frame[,2])
  Zhi <- mean(x@frame[,3]) + stats::sd(x@frame[,3])
  Zlo <- mean(x@frame[,3]) - stats::sd(x@frame[,3])
  Whi <- mean(x@frame[,4]) + stats::sd(x@frame[,4])
  Wlo <- mean(x@frame[,4]) - stats::sd(x@frame[,4])

  Y <- b0 +
    b1*c(Xhi, Xhi, Xhi, Xhi, Xlo, Xlo, Xlo, Xlo) +
    b2*c(Zhi, Zhi, Zlo, Zlo, Zhi, Zhi, Zlo, Zlo) +
    b3*c(Whi, Wlo, Whi, Wlo, Whi, Wlo, Whi, Wlo) +
    b4*c(Xhi, Xhi, Xhi, Xhi, Xlo, Xlo, Xlo, Xlo)*c(Zhi, Zhi, Zlo, Zlo, Zhi, Zhi, Zlo, Zlo) +
    b5*c(Xhi, Xhi, Xhi, Xhi, Xlo, Xlo, Xlo, Xlo)*c(Whi, Wlo, Whi, Wlo, Whi, Wlo, Whi, Wlo) +
    b6*c(Zhi, Zhi, Zlo, Zlo, Zhi, Zhi, Zlo, Zlo)*c(Whi, Wlo, Whi, Wlo, Whi, Wlo, Whi, Wlo) +
    b7*c(Xhi, Xhi, Xhi, Xhi, Xlo, Xlo, Xlo, Xlo)*c(Zhi, Zhi, Zlo, Zlo, Zhi, Zhi, Zlo, Zlo)*c(Whi, Wlo, Whi, Wlo, Whi, Wlo, Whi, Wlo)


  diff_a <- b5*(Whi - Wlo) + b7*Zhi*(Whi - Wlo)
  diff_b <- b4*(Zhi - Zlo) + b7*Whi*(Zhi - Zlo)
  diff_c <- b4*(Zhi - Zlo) + b7*Wlo*(Zhi - Zlo)
  diff_d <- b5*(Whi - Wlo) + b7*Zlo*(Whi - Wlo)
  diff_e <- b4*(Zhi - Zlo) + b5*(Whi - Wlo) + b7*(Zhi*Whi - Zlo*Wlo)
  diff_f <- b4*(Zhi - Zlo) + b5*(Wlo - Whi) + b7*(Zhi*Wlo - Zlo*Whi)

  sMat <- stats::vcov(x)[-1,-1]

  SE_a <- (Whi - Wlo) * sqrt(sMat[5,5] + (Zhi^2 * sMat[7,7]) + (2 * Zhi * sMat[5,7]))
  SE_b <- (Zhi - Zlo) * sqrt(sMat[4,4] + (Whi^2 * sMat[7,7]) + (2 * Whi * sMat[4,7]))
  SE_c <- (Zhi - Zlo) * sqrt(sMat[4,4] + (Wlo^2 * sMat[7,7]) + (2 * Wlo * sMat[4,7]))
  SE_d <- (Whi - Wlo) * sqrt(sMat[5,5] + (Zlo^2 * sMat[7,7]) + (2 * Zlo * sMat[5,7]))
  SE_e <- sqrt( (((Zhi - Zlo)^2) * sMat[4,4]) + (((Whi - Wlo)^2) * sMat[5,5]) +
                  (((Zhi*Whi - Zlo*Wlo)^2) * sMat[7,7]) +
                  2*(((Zhi-Zlo) * (Whi-Wlo) * sMat[4,5]) +
                       ((Zhi-Zlo) * (Zhi*Whi -Zlo*Wlo) * sMat[4,7]) +
                       ((Whi-Wlo) * (Zhi*Whi -Zlo*Wlo) * sMat[5,7])))
  SE_f <- sqrt( (((Zhi - Zlo)^2) * sMat[4,4]) + (((Wlo - Whi)^2) * sMat[5,5]) +
                  (((Zhi*Wlo - Zlo*Whi)^2) * sMat[7,7]) +
                  2*(((Zhi-Zlo) * (Wlo-Whi) * sMat[4,5]) +
                       ((Zhi-Zlo) * (Zhi*Wlo -Zlo*Whi) * sMat[4,7]) +
                       ((Whi-Wlo) * (Zhi*Wlo -Zlo*Whi) * sMat[5,7])))

  zWithin <- ifelse(sum(c(by(x@frame[,3], x@frame[,5], stats::var))) == 0, FALSE, TRUE)
  wWithin <- ifelse(sum(c(by(x@frame[,4], x@frame[,5], stats::var))) == 0, FALSE, TRUE)

  res <- data.frame("Y" = Y,
                    "X" = c("HIGH", "HIGH", "HIGH", "HIGH", "LOW", "LOW", "LOW", "LOW"),
                    "Z" = c("HIGH", "HIGH", "LOW", "LOW", "HIGH", "HIGH", "LOW", "LOW"),
                    "W" = c("HIGH", "LOW", "HIGH", "LOW", "HIGH", "LOW", "HIGH", "LOW"))

  dfs <- ifelse(zWithin == TRUE | wWithin == TRUE, nrow(x@frame)-8, length(unique(x@frame[,5]))-8)
  ts <- c(diff_a/SE_a, diff_b/SE_b, diff_c/SE_c, diff_d/SE_d, diff_e/SE_e, diff_f/SE_f)

  res2 <- data.frame("Label" = c("a", "b", "c", "d", "e", "f"),
                     "Value" = c(diff_a, diff_b, diff_c, diff_d, diff_e, diff_f),
                     "SE" = c(SE_a, SE_b, SE_c, SE_d, SE_e, SE_f),
                     "df" = rep(dfs, 6),
                     "t" = ts,
                     "p" = c(2*stats::pt(-abs(ts[1]),df=dfs),
                             2*stats::pt(-abs(ts[2]),df=dfs),
                             2*stats::pt(-abs(ts[3]),df=dfs),
                             2*stats::pt(-abs(ts[4]),df=dfs),
                             2*stats::pt(-abs(ts[5]),df=dfs),
                             2*stats::pt(-abs(ts[6]),df=dfs)))

  resList <- list("Expected_Y" = res, "SlopeTests" = res2)
  class(resList) <- "slopetest"
  if (zWithin == TRUE | wWithin == TRUE){
    warning("Normal calculation of df used. This may not be appropriate in some cases.")
  }
  return(resList)
}
