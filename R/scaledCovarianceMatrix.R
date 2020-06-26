#' Returns elements of the scaled variance-covariance martrix for a location-scale distribution.
#' @param distribution The probability distribution.
#' @param std.log.censor.time Vector of standardized log censoring times. 
#' @param unlist If TRUE and std.log.censor.time is a scalar, return a vector instead of a list.
#' @return Silently returns  a results matrix containing the variance factors.
#' @examples
#' @useDynLib lsinf
#' @export
scaledCovarianceMatrix <- function(distribution, std.log.censor.time, unlist=TRUE) {
    big.ones <- std.log.censor.time > 10^10
    v11 <- rep(NA, length(std.log.censor.time))
    v12 <- rep(NA, length(std.log.censor.time))
    v22 <- rep(NA, length(std.log.censor.time))
###
### handle the right-censored ones
###
    if (any(!big.ones)) {
        lsinf.out <- lsinf(std.log.censor.time[!big.ones], "right", distribution)
        det <- lsinf.out$f11 * lsinf.out$f22 - lsinf.out$f12^2
        v11[!big.ones] <- as.vector(lsinf.out$f22/det)
        v12[!big.ones] <- as.vector(-lsinf.out$f12/det)
        v22[!big.ones] <- as.vector(lsinf.out$f11/det)
    }
###
### handle the ones where the probability of failure is essentially 1
###
    if (any(big.ones)) {
        big.log.censor.time <- 10^100
        lsinf.out <- lsinf(big.log.censor.time, "uncensored", distribution)
        det <- lsinf.out$f11 * lsinf.out$f22 - lsinf.out$f12^2
        v11[big.ones] <- as.vector(lsinf.out$f22/det)
        v12[big.ones] <- as.vector(-lsinf.out$f12/det)
        v22[big.ones] <- as.vector(lsinf.out$f11/det)
    }
    result <- list(v11 = v11, v12 = v12, v22 = v22)
###
### return a vector if there was only one std.log.censor.time (and unlist==TRUE)
###
    if(length(result$v11)==1 && unlist) result <- unlist(result)
    return(result)
}
