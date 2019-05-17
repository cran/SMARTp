#' @title Sample size calculation under a clustered SMART design for non-surgical treatment of chronic periodontitis
#'
#' @description Sample size calculations to detect desired DTR effects, which includes (\eqn{i}) a single regime, (\eqn{ii}) difference between two regimes, and (\eqn{iii}) a specific regime is the best,
#' based on CAL changes under the proposed clustered, two-stage, SMART trial given type I and type II error rates
#' @usage SampleSize_SMARTp(mu, st1, dtr, regime, pow, a, rho, tau, sigma1,
#'                   lambda, nu, sigma0, Num, p_i, c_i, a0, b0, cutoff)
#'
#' @param mu Mean matrix, where row represents each treatment path from the SMART design diagram (see Xu \emph{et al.}, 2019), and
#' column represents each unit (i.e. tooth) within a cluster (i.e. mouth)
#' @param st1 Stage-1 treatment matrix, where rows represent the corresponding stage-1 treatments,
#' the 1st column includes the number of treatment options for the responder, the 2nd
#' column include the numbers of treatment options for the non-responder, the 3rd column
#' are the response rates, and the 4th column includes the row numbers
#' @param dtr Matrix of dimension (# of DTRs X 4), the 1st column represents the DTR numbers, the 2nd
#' column represents the treatment path number of responders for the corresponding DTRs in the
#' 1st column, the 3rd column represents the corresponding treatment path number of the non-responders
#' for the corresponding DTRs in the 1st column, while the 4th column represents
#' the corresponding initial treatment
#' @param regime Treatment regime vector. For detecting regime 1 as the best, use c(1, 2, 3, 4, 5, 6, 7, 8). Similarly, if regime 2 is the best, use c(2, 1, 3, 4, 5, 6, 7, 8), and so on
#' @param pow Power or 1 - Type II error rate, default is 0.8
#' @param a Type I error rate, default is 0.05
#' @param tau Variance parameter of the CAR model, default is 0.85
#' @param rho Association parameter of the CAR model, default is 0.975
#' @param sigma1 Standard deviation of the residual for the continuous outcome \eqn{Y_{it}}, default is 0.95
#' @param lambda Skewness parameter of the residual for the continuous outcome \eqn{Y_{it}}, default is 0
#' @param nu The degrees of freedom parameter of the residual for \eqn{Y_{it}}, default is Inf
#' @param sigma0 Standard deviation of the residual for the binary outcome \eqn{M_{it}}, default is 1
#' @param Num Iteration size to estimate variance of \eqn{\bar{Y}_i}, default is 100000
#' @param p_i The expected proportion of available teeth for subject \eqn{i}
#' @param c_i The average Pearson correlation coefficient between \eqn{Y_{it}} and \eqn{M_{it}} over the 28 teeth
#' @param a0 Intercept parameter in the probit model for the binary \eqn{M_{it}}, default is -1
#' @param b0 Slope parameter corresponding to the spatial random effect in the probit model for binary \eqn{M_{it}}, default is 0.5; note that \eqn{a_0} and \eqn{b_0} can be determined given \eqn{p_i} and \eqn{c_i}
#' @param cutoff Cut-off value of the binary outcome regression, default is 0
#'
#' @details SampleSize_SMARTp computes the sample size required to detect the dynamic treatment regime (DTR)
#' (Murphy, 2005, \emph{Statistics in Medicine}) effects in a study comparing non-surgical treatments of chronic periodontitis, via the sequential multiple
#' assignment randomized trial (SMART) design, with two-stages.
#'
#' Outcome measures (i.e. change in CAL) are continuous and clustered (i.e. tooth within a subject’s mouth,
#' where each subject/mouth is a cluster) with non-random missingness captured via a shared parameter setting, specified in
#' Reich and Bandyopadhyay (2010, \emph{Annals of Applied Statistics}). Each cluster sub-unit has a binary missingness indicator, which is
#' associated to its corresponding change of CAL through a joint model. The covariance structure within a cluster
#' is captured by the conditionally autoregressive (CAR) structure (Besag et al, 1991).
#'
#' The DTR effect can be detected based on either a single treatment regime, or the difference between two
#' treatment regimes (with or without sharing initial treatments), or when one regime is considered the best among others. The mean and variance of the CAL change for
#' each DTR can be estimated by the inverse probability weighting method via method of moments.
#'
#' Note that the first three inputs "mu", "st1" and "dtr" define the SMART design in term of matrices. From
#' Xu \emph{et al.} (2019+, Under Review), stage-1 includes two treatments, e.g., treatments "3" and "8". Participants who respond to the
#' stage-1 treatment will receive same treatment at stage-2, while non-responders will be randomly allocated to
#' other treatments, i.e. non-responders who received treatment "3" at stage-1 will be randomly allocated to
#' treatments "4"-"7" at stage-2, while non-responders receiving treatment "8" at stage-1 will be randomly
#' allocated to treatments "4"-"7" at stage-2.
#'
#' There are 8 treatment regimes for this design. They are 1 (treatment "3" at stage-1 and treatment "3" at stage-
#' 2 if responder, otherwise treatment "4"), 2 (treatment "3" at stage-1 and treatment "3" at stage-2 if responder,
#' otherwise treatment "5"), 3 (treatment "3" at stage-1 and treatment "3" at stage-2 if responder, otherwise
#' treatment "6"), 4 (treatment "3" at stage-1 and treatment "3" at stage-2 if responder, otherwise treatment "7"),
#' 5 (treatment "8" at stage-1 and treatment "8" at stage-2 if responder, otherwise treatment "4"), 6 (treatment "8" at
#' stage-1 and treatment "8" at stage-2 if responder, otherwise treatment "5"), 7 (treatment "8" at stage-1 and
#' treatment "8" at stage-2 if responder, otherwise treatment "6") and 8 (treatment "8" at stage-1 and treatment
#' "8" at stage-2 if responder, otherwise treatment "7"). See Figure 2 in Xu \emph{et al.} (2019+, Under Review)
#'
#' @return \item{N}{the estimated sample size}
#'\item{Del}{effect size}
#'\item{Del_std}{standardized effect size}
#'\item{ybar}{the estimated regime means corresponding to "regime"}
#'\item{Sigma}{the CAR covariance matrix corresponding to the latent \eqn{Q_{it}}; see Xu \emph{et al.} (2019+, Under Review)}
#'\item{sig.dd}{N*the variance or covariance matrix of the estimated regime means corresponding to "regime"}
#'\item{sig.e.sq}{N*the variance or covariance matrix of the difference between first and rest of estimated regime means corresponding to "regime", sig.e.sq = sig.dd if the element number of "regime" is one}
#'\item{p_st1}{the randomization probability of stage-1 for each treatment path}
#'\item{p_st2}{the randomization probability of stage-2 for each treatment path}
#'\item{res}{a vector with binary indicators represent responses or non-responses that corresponds to a treatment path}
#'\item{ga}{the response rates of initial treatments corresponding to each treatment path}
#'\item{initr}{column matrix with dimension = the number of treatment paths, the elements are the corresponding row number of st1}
#'
#' @author Jing Xu, Dipankar Bandyopadhyay, Douglas Azevedo, Bibhas Chakraborty
#'
#' @references Besag, J., York, J. & Mollie, A. (1991) \emph{"Bayesian image restoration, with two applications in spatial statistics
#' (with discussion)"}, Annals of the Institute of Statistical Mathematics 43, 159.
#'
#' Murphy, S. A. (2005), \emph{"An experimental design for the development of adaptive treatment strategies"},
#' Statistics in Medicine 24, 1455–1481.
#'
#' Reich, B. & Bandyopadhyay, D. (2010), \emph{A latent factor model for spatial data with informative missingness},
#' The Annals of Applied Statistics 4, 439–459.
#'
#' Xu, J., Bandyopadhyay, D., Mirzaei, S., Michalowicz, B and Bibhas Chakraborty. (2019+), \emph{"SMARTp: A SMART
#' design for non-surgical treatments of chronic periodontitis with spatially-referenced and non-randomly missing skewed outcomes"}, Under Review
#'
#' @seealso CAR_cov_teeth, MC_var_yibar_mis
#'
#' @examples m <- 28
#' pow <- 0.8
#' a <- 0.05
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
#'
#' Sigma <- CAR_cov_teeth(m = m, rho = rho, tau = tau)
#' p_i <- SMARTp:::pifun(cutoff = cutoff, a0 = a0, b0 = b0,
#'                      Sigma = Sigma, sigma0 = sigma0)
#' cit4 <- b0*diag(Sigma)/sqrt((diag(Sigma) +
#'             (sigma1^2 - 2/pi*sigma1^2*(0^2/(1+0^2))))*(b0^2*diag(Sigma) +
#'             sigma0^2))
#' c_i <- mean(cit4)
#'
#' del1 <- 5
#' del2 <- 0
#' del3 <- 0
#'
#' mu_sim <- matrix(0, 10, m)
#' mu_sim[2, ] <- rep(del1, m)
#' mu_sim[4, ] <- rep(del2, m)
#' mu_sim[7, ] <- rep(del3, m)
#'
#' st1 <- cbind(c(1, 1), c(4, 4), c(0.25, 0.5), 1:2)
#'
#' ##-- Stage-1 information
#' dtr <- cbind(1:8, c(rep(1, 4), rep(6, 4)),
#'              c(2, 3, 4, 5, 7, 8, 9, 10), c(rep(1, 4), rep(2, 4)))
#'
#' ##-- Detecting a single regime, e.g., Regime 1
#' regime <- 1
#' SampleSize <- SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr,
#'                                 regime = regime,
#'                                 pow = pow, a = a, rho = rho,
#'                                 tau = tau, sigma1 = sigma1, lambda = 0,
#'                                 nu = Inf, sigma0 = sigma0, Num = Num,
#'                                 p_i = p_i, c_i = c_i,
#'                                 cutoff = cutoff)
#'
#' N <- ceiling(SampleSize$N)
#' sig.e.sq <- SampleSize$sig.e.sq
#'
#' sqrt(diag(sig.e.sq)/N)
#' SampleSize$Del_std
#' SampleSize$Del
#' SampleSize$sig.dd
#' sqrt(diag(SampleSize$sig.dd)/N)
#' SampleSize$ybar
#'
#' ##-- Now using a0 and b0
#' SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr, regime = regime,
#'                   pow = pow, a = a, rho = rho,
#'                   tau = tau, sigma1 = sigma1, lambda = 0, nu = Inf,
#'                   sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                   cutoff = cutoff)
#' SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr, regime = regime,
#'                   p_i = p_i, c_i = c_i)
#'
#' ##-- Detecting the difference between two regimes that shares initial treatment,
#' ##-- e.g., Regimes 1 vs 3
#' regime <- c(1, 3)
#' SampleSize = SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr, regime = regime,
#'                                pow = pow, a = a, rho = rho,
#'                                tau = tau, sigma1 = sigma1, lambda = 0, nu = Inf,
#'                                sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                                cutoff = cutoff)
#'
#' N <- ceiling(SampleSize$N)
#' sig.e.sq <- SampleSize$sig.e.sq
#'
#' sqrt(diag(sig.e.sq)/N)
#' SampleSize$Del_std
#' SampleSize$Del
#' SampleSize$sig.dd
#'
#' ##-- Detecting the difference between two regimes that do not share initial treatment,
#' ##-- e.g., Regimes 1 vs 5
#' regime <- c(1, 5)
#'
#' SampleSize <- SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr, regime = regime,
#'                                 pow = pow, a = a, rho = rho,
#'                                 tau = tau, sigma1 = sigma1, lambda = 0, nu = Inf,
#'                                 sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                                 cutoff = cutoff)
#'
#' N <- ceiling(SampleSize$N)
#' sig.e.sq <- SampleSize$sig.e.sq
#'
#' sqrt(diag(sig.e.sq)/N)
#' SampleSize$Del_std
#' SampleSize$Del
#' SampleSize$sig.dd
#'
#' ##-- Detecting when Regime 1 is the best, e.g., comparing Regimes 1 vs 2, 3, 4, 5, 6, 7 and 8, i.e.
#' ##-- the alternative hypothesis is \mu_{d1}>\mu_{d2} & \mu_{d1}>\mu_{d3} ... & \mu_{d1}>\mu_{d8}
#' ##-- Note that this is a one-side test with Type-1 error rate of 0.025.
#' regime <- c(1, 2, 3, 4, 5, 6, 7, 8)

#' ##-- To detect Regime 2 is the best, just use regime = c(2, 1, 3, 4, 5, 6, 7, 8), and so on
#' SampleSize <- SampleSize_SMARTp(mu = mu_sim, st1 = st1, dtr = dtr, regime = regime,
#'                                 pow = pow, a = a, rho = rho,
#'                                 tau = tau, sigma1 = sigma1, lambda = 0, nu = Inf,
#'                                 sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0,
#'                                 cutoff = cutoff)
#'
#' N <- ceiling(SampleSize$N)
#' sig.e.sq <- SampleSize$sig.e.sq
#'
#' sqrt(diag(sig.e.sq)/N)
#' SampleSize$Del_std
#' SampleSize$Del
#' SampleSize$sig.dd
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @export

SampleSize_SMARTp = function(mu, st1, dtr,
                             regime, pow,
                             a,
                             rho, tau,
                             sigma1, lambda, nu, sigma0,
                             Num, p_i, c_i,
                             a0, b0,
                             cutoff){

  if(missing(regime)) stop('regime need to be specified')
  if(missing(pow)) pow <- 0.8
  if(missing(a)) a <- 0.05
  if(missing(rho)) rho <- 0.975
  if(missing(tau)) tau <- 0.85
  if(missing(sigma1)) sigma1 <- 0.95
  if(missing(lambda)) lambda <- 0
  if(missing(nu)) nu <- Inf
  if(missing(sigma0)) sigma0 <- 1
  if(missing(Num)) Num <- 100000
  if(missing(a0)) a0 <- -1
  if(missing(b0)) b0 <- 0.5
  if(missing(cutoff)) cutoff <- 0

  m <- dim(mu)[2]
  z.a.by.2 <- qnorm(1-a/2)
  Sigma <- CAR_cov_teeth(m = m, rho = rho, tau = tau)

  if(missing(p_i)){
    a0 <- a0
  } else{
    p_i <- pifun(cutoff, a0, b0, Sigma, sigma0)
    a0 <- uniroot(a0fun, c(-100, 100), tol = 0.0001, p_i = p_i, cutoff = cutoff, b0 = b0, Sigma = Sigma, sigma0 = sigma0)$root
  }

  if(missing(c_i)){
    b0 <- b0
  } else{
    b0 <- uniroot(b0fun, c(-1, 1), tol = 0.0001, c_i = c_i, Sigma = Sigma, sigma1 = sigma1, nu = nu, lambda = lambda, sigma0 = sigma0)$root
  }

  no.regime <- length(regime)

  initr <- as.matrix(rep(st1[, 4], st1[, 1] + st1[, 2]))
  ga <- as.matrix(rep(st1[, 3], st1[, 1] + st1[, 2]))
  res <- c()
  p_st2 <- c()

  for(i in 1:dim(st1)[1]){
    res <- c(res, rep(c(1, 0), st1[i,1:2]))
    p_st2 <- c(p_st2, rep(1/st1[i, 1:2], st1[i, 1:2]))
  }

  res <- as.matrix(res)
  p_st2 <- as.matrix(p_st2)

  ga_comp <- rep(0, no.regime)
  mu_comp <- array(0,c(2,m,no.regime))
  p2_comp <- matrix(0, 2, no.regime)
  p1_comp <- rep(0,no.regime)

  ga_comp <- ga[dtr[regime, 2],]

  for(j in 1:no.regime){
    mu_comp[, , j] <- rbind(mu[dtr[regime[j], 2],], mu[dtr[regime[j], 3],])
  }
  for(j in 1:no.regime){
    p2_comp[, j] <- rbind(p_st2[dtr[regime[j], 2], ], p_st2[dtr[regime[j], 3], ])
  }

  p_st1 <- as.matrix(rep((colSums(t(1/st1[, 1:2])*rbind(st1[, 3], (1-st1[, 3]))))^(-1)/sum((colSums(t(1/st1[, 1:2])*rbind(st1[, 3], (1-st1[, 3]))))^(-1)),
                         st1[, 1] + st1[, 2]))
  p1_comp <- p_st1[dtr[regime, 2], ]
  p1_comp <- as.vector(p1_comp)

  mud <- matrix(NA, no.regime, 2)
  Sigd <- matrix(NA, no.regime, 2)
  ybar <- rep(NA, no.regime)
  sig.dd <- matrix(NA, no.regime, no.regime)
  sig.e.sq <- matrix(NA, max(no.regime-1 , 1), max(no.regime-1, 1))

  for(j in 1:no.regime){
    YibardR <- MC_var_yibar_mis(mu = mu_comp[1, , j], Sigma = Sigma, sigma1 = sigma1,
                                lambda = lambda, nu = nu, sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0, cutoff = cutoff)
    YibardNR <- MC_var_yibar_mis(mu = mu_comp[2, , j], Sigma = Sigma, sigma1 = sigma1,
                                 lambda = lambda, nu = nu, sigma0 = sigma0, Num = Num, a0 = a0, b0 = b0, cutoff = cutoff)

    Sigd[j, 1] <- YibardR$VarYi
    mud[j, 1] <- YibardR$mYi

    Sigd[j, 2] <- YibardNR$VarYi
    mud[j, 2] <- YibardNR$mYi

    ybar[j] <- mud[j, 1]*ga_comp[j] + mud[j, 2]*(1-ga_comp[j])

    vr <- (ga_comp[j]/(p1_comp[j]*p2_comp[1, j]))*((Sigd[j, 1]) + (1-p1_comp[j]*p2_comp[1, j])*(mud[j, 1])^2)
    vnr <- ((1-ga_comp[j])/(p1_comp[j]*p2_comp[2, j]))*((Sigd[j, 2]) + (1-p1_comp[j]*p2_comp[2, j])*(mud[j, 2])^2)
    vrnrdiff <- ga_comp[j]*(1-ga_comp[j])*(mud[j, 1] - mud[j, 2])^2
    sig.dd[j,j] <- vr+vnr+vrnrdiff;
  }

  if(no.regime == 1){
    Del <- abs(ybar[1])
  } else{
    Del <- abs((rep(ybar[1], no.regime) - ybar)[-1])
  }

  for(i in 1:max(no.regime-1, 1)){
    for( j in min(i+1, no.regime):(no.regime-0)){
      if(no.regime == 1){
        covr <- sig.dd[1, 1]
      } else{
        if(dtr[regime[i], 4] == dtr[regime[j], 4]){
          covr <- (ga_comp[i]/(p1_comp[i]*p2_comp[1, i]))*((Sigd[i, 1]) + (mud[i, 1])^2) - (ga_comp[i]*ga_comp[j])*(mud[i, 1]*mud[j, 1]) - (ga_comp[i]*(1-ga_comp[j]))*(mud[i, 1]*mud[j, 2]) - ((1-ga_comp[i])*ga_comp[j])*(mud[i, 2]*mud[j, 1]) - ((1-ga_comp[i])*(1-ga_comp[j]))*(mud[i, 2]*mud[j, 2])
        } else{
          if(dtr[regime[i], 4] != dtr[regime[j], 4]){
            covr <- -(ga_comp[i]*ga_comp[j])*(mud[i, 1]*mud[j, 1]) - (ga_comp[i]*(1-ga_comp[j]))*(mud[i, 1]*mud[j, 2]) - ((1-ga_comp[i])*ga_comp[j])*(mud[i, 2]*mud[j, 1]) - ((1-ga_comp[i])*(1-ga_comp[j]))*(mud[i, 2]*mud[j, 2])
          }
        }
      }
      sig.dd[i,j]=covr
    }
  }

  for(j in 1:max(no.regime-1, 1)){
    sig.e.sq[j, j] <- sig.dd[1, 1] + ifelse(no.regime == 1, 0, sig.dd[j+1, j+1]) - 2*ifelse(no.regime == 1, 0, sig.dd[1, j+1])
  }

  for(i in 1:max(no.regime-2, 1)){
    for(j in ifelse(no.regime <= 2, 1, i+1):max(no.regime-1, 1)){
      sig.e.sq[i, j] <- sig.dd[1, 1] - ifelse(no.regime == 1, 0, sig.dd[1, (i+1)]) - ifelse(no.regime == 1, 0, sig.dd[1, (j+1)]) + ifelse(no.regime == 1, 0, sig.dd[(i+1), (j+1)] )
    }
  }

  sig.e.sq[lower.tri(sig.e.sq, diag = FALSE)] <- t(sig.e.sq)[lower.tri(sig.e.sq, diag = FALSE)]
  diag(sig.e.sq) <- diag(sig.e.sq) + 0.001

  sig.e.sq.cor <- cov2cor(sig.e.sq)
  Del_std <- Del/sqrt(diag(sig.e.sq)/2)
  ssize.fct <- function(N, Del_std, z.a.by.2, sig.e.sq.cor, pow){
    crit.vals <- z.a.by.2 - sqrt(N/2)*Del_std
    pmvnorm(upper = -as.vector(crit.vals), sigma = as.matrix(sig.e.sq.cor)) - pow
  }

  N <- uniroot(ssize.fct, interval = c(2, 1000), extendInt = "yes", Del_std = Del_std, z.a.by.2 = z.a.by.2, sig.e.sq.cor = sig.e.sq.cor, pow = pow)$root

  out <- list(N = N,
              sig.dd = sig.dd,
              sig.e.sq = sig.e.sq,
              Del = Del,
              Del_std = Del_std,
              ybar = ybar,
              initr = initr,
              ga = ga,
              res = res,
              p_st2 = p_st2,
              p_st1 = p_st1,
              Sigma = Sigma
  )

  class(out) <- "SMARTp"

  return(out)
}
