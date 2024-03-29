% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ismevExtension.R
\name{glod.fit}
\alias{glod.fit}
\title{Maximum Likelihood Fitting for the GLO distribution}
\usage{
glod.fit(
  xdat,
  ydat = NULL,
  mul = NULL,
  sigl = NULL,
  shl = NULL,
  mulink = identity,
  siglink = identity,
  shlink = identity,
  muinit = NULL,
  siginit = NULL,
  shinit = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  optimPars = NULL,
  maxit = 10000,
  fixedPars = list(mu = NULL, sig = NULL, sh = NULL),
  ...
)
}
\arguments{
\item{xdat}{A numeric vector of data to be fitted}

\item{ydat}{A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat}

\item{mul}{Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)}

\item{sigl}{As \code{mul} for the scale parameter}

\item{shl}{As \code{mul} for the shape parameter}

\item{mulink}{the link function for the location parameter - default to identity}

\item{siglink}{the link function for the scale parameter - default to identity}

\item{shlink}{the link function for the shape parameter - default to identity}

\item{muinit}{initial values for the location parameter}

\item{siginit}{initial values for the scale parameter}

\item{shinit}{initial values for the shape parameter}

\item{show}{Logical; if \code{TRUE} (the default), print details of the fit.}

\item{method}{The optimization method (see \code{optim} for details).}

\item{optimPars}{A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"}

\item{maxit}{The maximum number of iterations}

\item{fixedPars}{a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.}

\item{...}{Other control parameters for the optimization. These are passed to components of the control argument of optim.}
}
\value{
An object of the glo.fit class - with values which mirror the ones of the gev.fit class in ismev
}
\description{
Maximum-likelihood fitting for the generalized logistic distribution, 
including generalized linear modelling of each parameter. 
This function has the same structure as the \code{gevd.fit} and is inspired by \code{ismev::gev.fit}. 
The function allows any parameter to be kept fixed and to not be estimated.
}
\details{
The distribution is discussed in the Hosking and Wallis book and is used as the default distribution
for flood frequency estimation in the UK
}
\examples{
set.seed(12)
x <- runif(500)
y <- rglo(500,loc = 40+4*x,scale = 6,sh = 0.2)
fit1 <- glod.fit(y, show=FALSE)
fit1
## now add a regression model for the location
fit2 <- glod.fit(y, ydat = cbind(x), mul=1, show=FALSE)
fit2
## now a fit with a fixed shape parameter 
fitf <- glod.fit(y, show=FALSE, fixedPars = list(sh = 0.2))
fitf ## only two parameters are estimated 
}
\references{
Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 

Coles, S., 2001. An introduction to statistical modeling of extreme values. London: Springer.
}
\seealso{
\link{dglo}
}
