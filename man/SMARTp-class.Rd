% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMARTp-class.R
\docType{class}
\name{SMARTp-class}
\alias{SMARTp-class}
\title{An object of "SMARTp" class}
\description{
An object of "SMARTp" class
}
\section{Slots}{

\describe{
\item{\code{N}}{The estimated sample size}

\item{\code{sig.dd}}{N*the variance or covariance matrix of the estimated regime means correspond to "regime"}

\item{\code{sig.e.sq}}{N*the variance or covariance matrix of the difference between first and rest of estimated regime means correspond to \emph{regime}, sig.e.sq = sig.dd if the element number of \emph{regime} is one}

\item{\code{Del}}{Effect size}

\item{\code{Del_std}}{Standardized effect size}

\item{\code{ybar}}{The estimated regime means corresponding to "regime"}

\item{\code{initr}}{column matrix with dimension = the number of treatment paths, the elements are the corresponding row number of st1}

\item{\code{ga}}{The response rates of initial treatments corresponding to each treatment path}

\item{\code{res}}{A vector with binary indicators represent responders, or non-responders corresponding to a treatment path}

\item{\code{p_st1}}{The randomization probability of stage-1 for each treatment path}

\item{\code{p_st2}}{The randomization probability of stage-2 for each treatment path}

\item{\code{Sigma}}{The CAR covariance matrix of the latent \eqn{Q_{it}}}
}}

