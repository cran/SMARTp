#' An object of "SMARTp" class
#'
#' @slot N The estimated sample size
#' @slot sig.dd N*the variance or covariance matrix of the estimated regime means correspond to "regime"
#' @slot sig.e.sq N*the variance or covariance matrix of the difference between first and rest of estimated regime means correspond to \emph{regime}, sig.e.sq = sig.dd if the element number of \emph{regime} is one
#' @slot Del Effect size
#' @slot Del_std Standardized effect size
#' @slot ybar The estimated regime means corresponding to "regime"
#' @slot initr column matrix with dimension = the number of treatment paths, the elements are the corresponding row number of st1
#' @slot ga The response rates of initial treatments corresponding to each treatment path
#' @slot res A vector with binary indicators represent responders, or non-responders corresponding to a treatment path
#' @slot p_st1 The randomization probability of stage-1 for each treatment path
#' @slot p_st2 The randomization probability of stage-2 for each treatment path
#' @slot Sigma The CAR covariance matrix of the latent \eqn{Q_{it}}
#'
#' @importFrom methods new
#'
#' @export

`SMARTp-class` <- setClass(Class = 'SMARTp',
                           representation(N = 'numeric',
                                          sig.dd = 'numeric',
                                          sig.e.sq = 'numeric',
                                          Del = 'numeric',
                                          Del_std = 'numeric',
                                          ybar = 'numeric',
                                          initr = 'matrix',
                                          ga = 'matrix',
                                          res = 'matrix',
                                          p_st2 = 'matrix',
                                          p_st1 = 'matrix',
                                          Sigma = 'matrix'))
