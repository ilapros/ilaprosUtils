% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snippets.R
\name{meanNOTna}
\alias{meanNOTna}
\title{Mean without NAs}
\usage{
meanNOTna(x)
}
\arguments{
\item{x}{a vector}
}
\value{
The mean of the vector when NAs are excluded, just \code{mean(x, na.rm=TRUE)}
}
\description{
Mean without NAs
}
\examples{
x <- c(5,12,NA,4)
meanNOTna(x)
}
