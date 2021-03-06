% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ismevExtension.R
\name{print.gevcv.fit}
\alias{print.gevcv.fit}
\title{Nicer print of gevcv.fit results}
\usage{
\method{print}{gevcv.fit}(x, ...)
}
\arguments{
\item{x}{a fitted object of the class gevcv.fit.}

\item{...}{further arguments passed to \code{print}.}
}
\description{
This function prints the MLE, se, convergence info and negative log-likelihood value.
}
\examples{
print(gevcvd.fit(c(50,45,65,78,12,23),show=FALSE))
}
\seealso{
gevcvd.fit, isemv::gev.fit
}
