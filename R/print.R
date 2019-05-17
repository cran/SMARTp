#' Print for SMARTp class
#'
#' @param x SMARTp object to print
#' @param ... Other parameters for print
#'
#' @method print SMARTp
#'
#' @export

print.SMARTp <- function(x, ...){
  cat('Sample size: \n', sprintf("\t The sample size required is %s.", ceiling(x$N)))
}
