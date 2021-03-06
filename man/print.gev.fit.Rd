% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ismevExtension.R
\name{print.gev.fit}
\alias{print.gev.fit}
\title{Nicer print of gev.fit results}
\usage{
\method{print}{gev.fit}(x, ...)
}
\arguments{
\item{x}{a fitted object of the class gev.fit.}

\item{...}{further arguments passed to \code{print}.}
}
\description{
This function prints the MLE, se, convergence info and negative log-likelihood value.
}
\examples{
# library(ismev)
print(ismev::gev.fit(c(50,45,65,78,12,23),show=FALSE))
print(gevd.fit(c(50,45,65,78,12,23),show=FALSE))
}
\seealso{
gevd.fit, isemv::gev.fit
}
