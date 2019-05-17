#' @title The within-mouth covariance matrix with conditional autoregressive structure
#'
#' @description The covariance matrix of individual teeth measures for each subject follows a Conditional Autoregressive model (CAR) density
#' @usage CAR_cov_teeth(m, rho, tau)
#'
#' @param m Maximum number of units in each cluster, i.e., 28 teeth in each mouth (the 4 third-molars are usually ignored)
#' @param tau Variation parameter of the CAR model
#' @param rho Association parameter of the CAR model
#'
#' @details \emph{CAR_cov_teeth} gives the covariance matrix among the teeth within each mouth based on the CAR structure
#' (Besag \emph{et al.}, 1991), given the maximum number of teeth for each subject (\eqn{m}), the variance (\eqn{\tau}), and the
#' association (\eqn{\rho}) parameters.
#'
#' The CAR covariance matrix can be expressed as \eqn{\Sigma_{28\times 28} = \tau^2 (W - \rho D)^{-1}}, where \eqn{\tau^2 > 0}, and \eqn{\rho \in [0, 1]} are the
#' parameters that control the magnitude of variation and the degree of spatial association, respectively. For
#' matrix \eqn{D}, the element \eqn{D_{tt'}} is 1 if locations \eqn{t} and \eqn{t'} are adjacent and 0 otherwise. The matrix \eqn{W} is diagonal
#' with diagonal elements \eqn{W_{tt} = \sum_{t'} D_{tt'}}. Note, the argument \eqn{\tau} in \emph{CAR_cov_teeth} is the variance, and not the standard deviation. 
#'
#' @return The covariance matrix among the teeth in each mouth (assuming full dentition, i.e., 28 teeth) based on a CAR model. 
#'
#' @author Jing Xu, Dipankar Bandyopadhyay, Douglas Azevedo, Bibhas Chakraborty
#'
#' @references Besag, J., York, J. & Mollie, A. (1991), \emph{"Bayesian image restoration, with two applications in spatial statistics
#' (With Discussion)"}, Annals of the Institute of Statistical Mathematics 43, 159.
#'
#' Reich, B. & Bandyopadhyay, D. (2010), \emph{"A latent factor model for spatial data with informative missingness"},
#' The Annals of Applied Statistics 4, 439–459.
#'
#' @seealso MC_var_yibar_mis, SampleSize_SMARTp
#'
#' @examples m <- 28
#' rho <- 0.975
#' tau <- 0.85
#' Sigma <- CAR_cov_teeth(m = m, rho = rho, tau = tau)
#'
#' @export

CAR_cov_teeth <- function(m, rho, tau){
  ##-- Tests
  if(!is.numeric(m)) stop("m must be an integer number.")
  if(abs(rho) > 1) stop("rho must be a number in the (0, 1) interval.")
  if(!is.numeric(tau) | tau < 0) stop("tau must be a positive number.")

  if(m > 1){
    D_m <- abs(outer(1:m, 1:m, "-")) == 1
    M_m <- diag(colSums(D_m))
    Sigma <- tau^2*solve(M_m-rho*D_m)
  } else{
    Sigma <- tau^2
  }

  return(Sigma)
}
