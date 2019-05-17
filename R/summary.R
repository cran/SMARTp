#' @title Summary for SMARTp class
#'
#' @param object SMARTp object to summarise
#' @param ... Other parameters for summary
#'
#' @method summary SMARTp
#'
#' @export

summary.SMARTp <- function(object, ...){
  cat('Sample size: \n', sprintf("\t The sample size required is %s.", ceiling(object$N)))
}
