#' @title Estimated mean and variance of the average change in CAL for each subject
#'
#' @description The estimated Monte Carlo mean and variance of the average change in clinical attachment level (CAL) for each subject
#' @usage MC_var_yibar_mis(mu, Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
#'
#' @param mu Mean matrix, where row represents each treatment path, and column represents each cluster unit
#' @param Sigma Within-mouth teeth covariance matrix
#' @param sigma1 Standard deviation of the residual for the continuous outcome \eqn{Y_{it}}
#' @param lambda The skewness parameter of the residual for the continuous outcome \eqn{Y_{it}}
#' @param nu The degree freedom, or kurtosis parameter of the residual for the continuous outcome \eqn{Y_{it}}
#' @param sigma0 Standard deviation of the residual for the binary outcome \eqn{M_{it}}
#' @param Num Number of samples to estimate mean or variance of \eqn{\bar{Y}_i}
#' @param a0 Intercept parameter in the probit model for the binary outcome \eqn{M_{it}}
#' @param b0 Slope parameter corresponding to the spatial random effect in the probit model for the binary outcome \eqn{M_{it}}
#' @param cutoff Cut-off value in the binary outcome regression
#'
#' @details \emph{MC_var_yibar_mis} computes the Monte-Carlo estimates of expectation and variance of the sample mean among the teeth within each mouth, i.e
#' \eqn{\bar{Y}_i = \sum Y_{it}(1 - M_{it})/\sum(1 - M_{it})}, where \eqn{Y_{it}} is the change in CAL (measured in mm) for patient \eqn{i} and tooth \eqn{t}, and \eqn{M_{it}}
#' is the misingness indicator, i.e., \eqn{M_{it} = 1} implies tooth \eqn{t} in subject \eqn{i} is mising. The joint regression models for \eqn{Y_{it}} and \eqn{M_{it}}
#' are available in Reich & Bandyopadhyay (2010, \emph{Annals of Applied Statistics}).
#'
#' @return The simulated dataset of CAL change "\eqn{Y_{it}}", missingness "\eqn{M_{it}}" and function inside the indicator of "\eqn{M_{it} I_{it}}" for
#' each tooth of each patient, with the corresponding estimated mean "\eqn{mY_i}", variance "\eqn{VarY_i}" and missing proportion "PM" for each patient
#'
#' @author Jing Xu, Dipankar Bandyopadhyay, Douglas Azevedo, Bibhas Chakraborty
#'
#' @references Besag, J., York, J. & Mollie, A. (1991), \emph{"Bayesian image restoration, with two applications in spatial statistics
#' (With Discussion)"}, Annals of the Institute of Statistical Mathematics 43, 159.
#'
#' Reich, B. & Bandyopadhyay, D. (2010), \emph{"A latent factor model for spatial data with informative missingness"},
#' The Annals of Applied Statistics 4, 439–459.
#'
#' @seealso CAR_cov_teeth, SampleSize_SMARTp
#'
#' @examples m <- 28
#' Num <- 1000
#' cutoff <- 0
#' sigma1 <- 0.95
#' sigma0 <- 1
#' lambda <- 0
#' nu <- Inf
#' b0 <- 0.5
#' a0 <- -1.0
#' rho <- 0.975
#' tau <- 0.85
#' del1 <- 0.5
#' del2 <- 2
#'
#' Sigma <- CAR_cov_teeth(m, rho, tau)
#' Sigma_comp <- array(Sigma, c(m, m, 4))
#' Sigma_sim <- array(Sigma, c(m, m, 10))
#'
#' mu_comp <- array(0, c(2, m, 2))
#' mu_comp[, , 1] <- rbind(rep(0, m), rep(del1, m))
#' mu_comp[, , 2] <- rbind(rep(0, m), rep(del2, m))
#'
#' VarYitd1R = MC_var_yibar_mis(mu = mu_comp[1, , 1], Sigma = Sigma,
#'                              sigma1 = sigma1,
#'                              lambda = lambda, nu = nu,
#'                              sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                              cutoff = cutoff)
#' PM <- VarYitd1R$PM
#' VarYid1R <- VarYitd1R$VarYi
#' mYid1R <- VarYitd1R$mYi
#' VarYitd1NR <- MC_var_yibar_mis(mu = mu_comp[2, , 1], Sigma = Sigma,
#'                                sigma1 = sigma1,
#'                                lambda = lambda, nu = nu,
#' sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0, cutoff = cutoff)
#'
#' PM <- VarYitd1NR$PM
#' VarYid1NR <- VarYitd1NR$VarYi
#' mYid1NR <- VarYitd1NR$mYi
#' VarYitd3R <- MC_var_yibar_mis(mu = mu_comp[1, , 2], Sigma = Sigma,
#'                               sigma1 = sigma1,
#'                               lambda = lambda, nu = nu,
#'                               sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                               cutoff = cutoff)
#'
#' PM <- VarYitd3R$PM
#' VarYid3R <- VarYitd3R$VarYi
#' mYid3R <- VarYitd3R$mYi
#' VarYitd3NR <- MC_var_yibar_mis(mu = mu_comp[2,,2], Sigma = Sigma,
#'                                sigma1 = sigma1,
#'                                lambda = lambda, nu = nu,
#' sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0, cutoff = cutoff)
#'
#' PM <- VarYitd3NR$PM
#' VarYid3NR <- VarYitd3NR$VarYi
#' mYid3NR <- VarYitd3NR$mYi
#'
#' @importFrom sn rst
#' @importFrom mvtnorm rmvnorm
#' @import stats
#'
#' @export

MC_var_yibar_mis <- function(mu, Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff){
  m <- dim(Sigma)[1]

  Qit <- rmvnorm(Num, rep(0, m), Sigma)
  Yit <- matrix(1, nrow = Num, ncol = 1)%*%matrix(mu, nrow = 1, ncol = m) + Qit + matrix(rst(n = Num*m, xi = 0, omega = sigma1, alpha = lambda, nu = nu), nrow = Num, ncol = m)

  Iit <- a0 + Qit*b0 + rmvnorm(Num, rep(0, m), diag(sigma0^2, m))
  Mit <- ifelse(Iit > cutoff, 1, 0)

  mYi <- mean(rowSums(Yit*(1 - Mit))/rowSums(1-Mit), na.rm = T)
  VarYi <- var(rowSums(Yit*(1 - Mit))/rowSums(1-Mit), na.rm = T)

  PM <- mean(rowSums(Mit)/m)

  return(list(Yit = Yit, Mit = Mit, Iit = Iit, mYi = mYi, VarYi = VarYi, PM = PM))
}
