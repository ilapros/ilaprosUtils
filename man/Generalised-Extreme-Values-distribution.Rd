% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gevglo.R
\name{Generalised Extreme Values distribution}
\alias{Generalised Extreme Values distribution}
\alias{dgev}
\alias{pgev}
\alias{rgev}
\alias{qgev}
\title{The Generalised Extreme Values distribution}
\usage{
dgev(x, loc, scale, sh, log = FALSE)

pgev(q, loc, scale, sh, lower.tail = TRUE, log.p = FALSE)

qgev(p, loc, scale, sh, lower.tail = TRUE, log.p = FALSE)

rgev(n, loc, scale, sh)
}
\arguments{
\item{x, q}{vector of quantiles}

\item{loc}{location parameter}

\item{scale}{scale parameter}

\item{sh}{shape parameter}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p)}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}}

\item{p}{vector of probabilities}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\value{
dgev gives the density, pgev gives the distribution function, qgev gives the quantile function, and rgev generates random deviates.
The length of the result is determined by n for rgev, and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments are recycled to the length of the result. 
Only the first elements of the logical arguments are used.
}
\description{
Density, distribution function, quantile function and random generation for 
the Generalised Extreme Values distribution with location parameter equal to \code{loc}, 
scale parameter equal to \code{scale} and shape parameter equal to \code{sh}. 
The functions use the Hosking and Wallis notation so that the domain of the distribution 
has an upper bound when the shape parameter is positive.
}
\details{
Extra care should be taken on the shape parameter of the distribution.
Different notations are used in the scientific literature: in one notation, presented 
for example in the Hosking and Wallis book and common in hydrology, 
the domain of the distribution has an upper bound when the shape parameter is positive. 
Conversely, in one notation, presented for example in Coles' book and Wikipedia and 
common in statistics, the domain of the distribution has a lower bound when the shape parameter is positive. 
The two notation only differ for the sign of the shape parameter. 
The functions in this package use the Hosking and Wallis notation for consistency with the \code{\link{pglo}} functions. 
Nevertheless the fitting in the \code{gev.fit} function is based on the \code{isemv::gev.fit} function written by Stuart Coles which uses the 
Coles' notation: be aware of these differences!
}
\examples{
plot(seq(-15,40,by=0.2),dgev(seq(-15,40,by=0.2),4,6,0.2),type="l")
plot(ecdf(rgev(100,4,6,0.2)))
lines(seq(-15,40,by=0.5),pgev(seq(-15,40,by=0.5),4,6,0.2),col=2)
qgev(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2) 
# notable quantiles
}
\references{
Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press.
}
\concept{gev distribution}
\keyword{Extreme}
\keyword{Generalised}
\keyword{Values}
\keyword{distribution}
