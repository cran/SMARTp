b0fun <- function(c_i, b0, Sigma, sigma1, nu, lambda, sigma0){
  if(nu < Inf){
    mean(b0*diag(Sigma)/sqrt((diag(Sigma) + (sigma1^2*nu/(nu - 2) - nu/pi*(gamma(0.5*(nu - 1))/gamma(0.5*nu))^2*sigma1^2*(lambda^2/(1 + lambda^2))))*(b0^2 * diag(Sigma) + sigma0^2))) - c_i
  } else{
    mean(b0*diag(Sigma)/sqrt((diag(Sigma) + (sigma1^2 - 2/pi*sigma1^2*(lambda^2/(1 + lambda^2))))*(b0^2 * diag(Sigma) + sigma0^2))) - c_i
  }
}

pifun <- function(cutoff, a0, b0, Sigma, sigma0){
  Epit <- rep(0, 28)
  for(j in 1:28){
    Epit[j] <- stats::pnorm((cutoff - a0)/sqrt(b0^2 * Sigma[j, j] + sigma0^2))
  }
  res <- mean(Epit)
  return(res)
}

a0fun <- function(p_i, cutoff, a0, b0, Sigma, sigma0){
  res <- pifun(cutoff, a0, b0, Sigma, sigma0) - p_i
  return(res)
}
