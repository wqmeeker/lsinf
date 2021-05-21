#' Returns sqrt of the determinant of the FIM.
#' @param distribution The probability distribution.
#' @param prob.range range of probabilities for evaluation. 
#' @return Silently returns  a results vectors.
#' @examples
#' JeffreysPriorPlot()
#' @useDynLib lsinf
#' @export
###
"JeffreysPriorPlot" <- function(distribution="normal", prob.range = c(0.01, 0.99), number.points=100, cex.axis=1.5, cex.lab=1.5) {
###
### this first cut works only for the normal distribution
###
  z.range <- qnorm(prob.range)
  std.log.censor.time <- seq(z.range[1], z.range[2], length=number.points)
    big.ones <- std.log.censor.time > 10^10
###
### initalize
###
    det <- rep(NA, length(std.log.censor.time))
###
### handle the right-censored ones
###
    if (any(!big.ones)) {
        lsinf.out <- lsinf(std.log.censor.time[!big.ones], "right", distribution)
        det[!big.ones] <- lsinf.out$f11 * lsinf.out$f22 - lsinf.out$f12^2
     }
###
### handle the ones where the probability of failure is essentially 1
###
    if (any(big.ones)) {
        big.log.censor.time <- 10^100
        lsinf.out <- lsinf(big.log.censor.time, "uncensored", distribution)
        det[big.ones] <- lsinf.out$f11 * lsinf.out$f22 - lsinf.out$f12^2
    }
    result <- sqrt(det)
  par(mar=c(5.1, 5.1, 4.1, 2.1))
  plot(range(std.log.censor.time), range(result), xaxt = "n", type = "n", xlab="Expected fraction failing", ylab="Sqrt Det scaled FIM", las=1, cex.axis=cex.axis, cex.lab=cex.lab)
###
###  probability axis
###
  prob.labels=c("0.01", "0.05", "0.20", "0.50", "0.80", "0.95", "0.99")
  prob.locations= qnorm(as.numeric(prob.labels))
  axis(side = 1, at = prob.locations, labels = prob.labels, cex.axis=cex.axis)
  lines(std.log.censor.time, result, type="l", lwd=2)
    invisible(list(result, std.log.censor.time))
}
