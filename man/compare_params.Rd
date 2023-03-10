% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dabayes.R
\name{compare_params}
\alias{compare_params}
\title{Pairwise Comparisons of MCMC Output}
\usage{
compare_params(samps, hdi_prob = 0.95, ...)
}
\arguments{
\item{samps}{A matrix of MCMC output or an \code{mcmc.list}.}

\item{hdi_prob}{The probability to be covered by the HDI, defaults to 0.95.}

\item{...}{Other arguments to be passed down, currently unimplemented.}
}
\description{
Generates pairwise comparisons of MCMC output that can then be either
summarised or plotted. Creates differences by calculating the difference of
lhs_name - rhs_name.
}
