#' Returns elements of the scaled fisher information matrix for location-scale distributions
#' @param z standardized time
#' @param censor.type A characcter string: uncensored, left, or right
#' @param distribution A characcter string: sev, normal, lev, logistic.
#' @return A list containing vectors f11, f12 and f22 of Fisher information matrix elements.
#' @examples
#' lsinf(0.5, "right", "sev")
#' lsinf(0.5, "left", "normal"
#' lsinf(seq(-1, 1, by=.1), "right", "sev")
#' lsinf(seq(-2, 2, by=.2), "right", "normal")
#' @useDynLib lsinf
#' @export
lsinf <- function(z, censor.type, distribution) {
    if (!is.numeric(z)) 
        wqm.stop("z must be numeric")
    if (!is.character(censor.type)) 
        wqm.stop("censor.type must be character")
    if (!is.character(distribution)) 
        wqm.stop("distribution must be character string")
    switch(generic.distribution(distribution),
           weibull = , sev = idist <- 1,
           frechet = , lev = idist <- 2,
           lognormal = , normal = idist <- 3,
           loglogistic = , logistic = idist <- 4, 
        wqm.stop("Distribution must be sev, lev, normal, or logistic"))
    switch(censor.type,
           uncensored = icode <- 1,
           right = icode <- 2,
           left = icode <- 3, 
        wqm.stop("censor.type must be uncensored, left, or right"))
    nrows <- length(z)
    zout <- .Fortran("slsinf", as.integer(idist), as.integer(icode), as.double(z), as.double(z), 
        f11 = double(nrows), f12 = double(nrows), f22 = double(nrows), as.integer(nrows), 
        ifault = integer(1), irow = integer(1))
    if (zout$ifault > 1) 
        wqm.warning(paste(ier, "evaluation error in row", irow, sep = "", collapse = ""))
    return(list(f11 = zout$f11, f12 = zout$f12, f22 = zout$f22))
}
