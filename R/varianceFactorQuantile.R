#' @useDynLib lsinf
#' @export
varianceFactorQuantile <- function(quant.of.interest, proportion.failing, distribution) {
###
###  length of quant.of.interest and proportion.failing must be the same
###  unless one of them has length 1
###
### do the check
###
  if(length(quant.of.interest) != length(proportion.failing)){
    if(!(length(quant.of.interest)==1||length(proportion.failing)==1))
    stop(paste("if length of quant.of.interest",length(quant.of.interest), "is not equal to the length of proportion.failing", length(proportion.failing), "the one of them should have length 1"))
  }
###
### do the work
###
    std.log.censor.time <- quant(as.numeric(proportion.failing), distribution)
    std.quant.of.int <- quant(quant.of.interest, distribution)
    var.out <- scaledCovarianceMatrix(distribution, std.log.censor.time, unlist=FALSE)
    scaled.asvar <- var.out$v11 + 2 * std.quant.of.int * var.out$v12 + (std.quant.of.int)^2 *  var.out$v22
    return(scaled.asvar)
}
