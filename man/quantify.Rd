% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sce_tools.R
\name{quantify}
\alias{quantify}
\alias{quantify.numeric}
\alias{quantify.forecasthistogram}
\title{Fit parametric distribution to histogram probabilities}
\usage{
\method{quantify}{numeric}(x, ub, fit_support = TRUE, support_limit = 38, ...)

\method{quantify}{forecasthistogram}(x, fit_support = TRUE, support_limit = 38, ...)
}
\arguments{
\item{x}{object of class `forecasthistogram'}

\item{ub}{vector of upper bounds}

\item{fit_support}{whether to choose support according to statistical criterium (default is \code{TRUE})}

\item{support_limit}{defaults to 38 (relevant only if \code{fit_support == FALSE})}

\item{...}{other inputs (currently not in use)}
}
\value{
object of class `forecasthistogram'. Use function \code{quantile.forecasthistogram} to compute quantiles, \code{plot.forecasthistogram} for plotting, and
\code{mean.forecasthistogram} for computing the implied mean.
}
\description{
Fit parametric distribution to histogram probabilities
}
\details{
The function fits a parametric distribution based on the histogram probabilities in \code{p}. The fitting procedure has been proposed by
Engelberg et al. (2009). Briefly, the method fits a generalized beta distribution if \code{p} features 3 or more bins,
and fits a simpler triangular distribution if \code{p} features at most two bins. See Krüger and Pavlova (2020, Appendix A) for details on the
present implementation.
}
\references{
Engelberg, J., Manski, C.F., and J. Williams (2009): `Comparing the Point Predictions
and Subjective Probability Distributions of Professional Forecasters', Journal of Business & Economic Statistics 27,
30-41, \doi{10.1198/jbes.2009.0003}

Krüger, F., and L. Pavlova (2020): `Quantifying Subjective
Uncertainty in Survey Expectations'. Working Paper, \doi{10.18452/21401}.
}
