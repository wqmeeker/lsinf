#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("lsinf version ", packageVersion("lsinf"), " loaded")
    packageStartupMessage(paste("Location:", system.file(package = "lsinf")))
    packageStartupMessage("Send bug reports to wqmeeker@iastate.edu")
}
