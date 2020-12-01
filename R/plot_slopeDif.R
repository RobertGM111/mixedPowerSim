#' Create an interaction plot from an object of class "slopetest".
#'
#' Creates an interaction plot comparing high and low values of X at high and low value combinations of Z and W.
#'
#'
#' @usage plot(x, legLoc = "topright", ...)
#' @aliases plot
#' @param x An object of class "slopetest". The result of an \code{slopeDif()} function.
#' @param legLoc A string. Where should the legend be located? Default = "topright".
#' @param ... Additional arguments passed to \code{plot()}.
#' @return A three-way interaction plot
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
#'  mod <- lmer(Y~X * Z * W + (1|Cluster), data = simDat$simulatedData)
#'
#'  difTest <- slopeDif(mod)
#'
#'  plot(difTest)
#'  ## OR
#'  plot(difTest)
#' @export
plot.slopetest <- function(x, legLoc = "topright", ...){
  graphics::plot(c(1,1,1,1,9,9,9,9),
                 c(x$Expected_Y[5,1],
                   x$Expected_Y[6,1],
                   x$Expected_Y[7,1],
                   x$Expected_Y[8,1],
                   x$Expected_Y[1,1],
                   x$Expected_Y[2,1],
                   x$Expected_Y[3,1],
                   x$Expected_Y[4,1]),
                 ylab = "Y",
                 xlab = NA,
                 xaxt = "n",
                 bty = "n",
                 xlim = c(0,10),
                 type = "n", ...)

  graphics::axis(1, at = c(1,9),
                 labels = c("Low X", "High X"))

  graphics::lines(c(1,9),
                  c(x$Expected_Y[5,1],
                    x$Expected_Y[1,1]),
                  lwd = 3)
  graphics::lines(c(1,9),
                  c(x$Expected_Y[6,1],
                    x$Expected_Y[2,1]),
                  lwd = 3)
  graphics::lines(c(1,9),
                  c(x$Expected_Y[7,1],
                    x$Expected_Y[3,1]),
                  lwd = 3)
  graphics::lines(c(1,9),
                  c(x$Expected_Y[8,1],
                    x$Expected_Y[4,1]),
                  lwd = 3)

  graphics::points(c(1,1,1,1,9,9,9,9),
                   c(x$Expected_Y[5,1],
                     x$Expected_Y[6,1],
                     x$Expected_Y[7,1],
                     x$Expected_Y[8,1],
                     x$Expected_Y[1,1],
                     x$Expected_Y[2,1],
                     x$Expected_Y[3,1],
                     x$Expected_Y[4,1]),
                   pch = c(22, 23, 22, 23, 22, 23, 22, 23),
                   bg = c("black","black","white","white","black","black","white","white"),
                   cex = 2)
  if (!is.na(legLoc)){
    graphics::legend(legLoc,
                     legend = c("Low W Low Z",
                                "Low W High Z",
                                "High W Low Z",
                                "High W High Z"),
                     lty = c(1,1,1,1),
                     lwd = c(2,2,2,2),
                     pch = c(23, 23, 22, 22),
                     pt.bg = c("white","black","white","black"))
  }
}
