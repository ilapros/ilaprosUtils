% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kappa.R
\name{kappad.fit}
\alias{kappad.fit}
\title{Maximum Likelihood Fitting for the Kappa distribution}
\usage{
kappad.fit(
  xdat,
  ydat = NULL,
  mul = NULL,
  sigl = NULL,
  shl = NULL,
  sh2l = NULL,
  mulink = identity,
  siglink = identity,
  shlink = identity,
  sh2link = identity,
  muinit = NULL,
  siginit = NULL,
  shinit = NULL,
  sh2init = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  optimPars = NULL,
  maxit = 10000,
  fixedPars = list(mu = NULL, sig = NULL, sh = NULL, sh2 = NULL),
  ...
)
}
\arguments{
\item{xdat}{A numeric vector of data to be fitted}

\item{ydat}{A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat}

\item{mul}{Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)}

\item{sigl}{As \code{mul} for the scale parameter}

\item{shl}{As \code{mul} for the shape parameter}

\item{sh2l}{As \code{mul} for the second shape parameter}

\item{mulink}{the link function for the location parameter - default to identity}

\item{siglink}{the link function for the scale parameter - default to identity}

\item{shlink}{the link function for the shape parameter - default to identity}

\item{sh2link}{the link function for the second shape parameter - default to identity}

\item{muinit}{initial values for the location parameter}

\item{siginit}{initial values for the scale parameter}

\item{shinit}{initial values for the shape parameter}

\item{sh2init}{initial values for the second shape parameter}

\item{show}{Logical; if \code{TRUE} (the default), print details of the fit.}

\item{method}{The optimization method (see \code{optim} for details).}

\item{optimPars}{A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"}

\item{maxit}{The maximum number of iterations.}

\item{fixedPars}{a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.}

\item{...}{Other control parameters for the optimization. These are passed to components of the control argument of optim.}
}
\value{
An object of the kappa.fit class - with values which mirror the ones of the gev.fit class in \code{ismev}
}
\description{
Maximum-likelihood fitting for the Kappa distribution, 
including generalized linear modelling of each parameter. 
This function has the same structure as the \code{gevd.fit} and is inspired by \code{ismev::gev.fit}. 
The function allows any parameter to be kept fixed and to not be estimated.
}
\details{
The distribution is discussed in the Hosking and Wallis book and can be seen as a generelisation of several distributions used in extreme value modelling. 
Depending on the sh2 parameter value the Kappa distribution reduces to others commonly used distributions For example: 
when sh2 is equal to -1 the distribution reduces to a GLO, when sh2 is equal to 0 the distribution reduces to a GEV, 
when sh2 is equal to +1 the distribution reduces to a Generalised Pareto distribution (GPA).
}
\examples{
set.seed(12)
x <- runif(300)
y <- rkappa(300,loc = 40+4*x,scale = 6, sh = 0.2, sh2=-0.4)
fit1 <- kappad.fit(y, show=FALSE)
fit1
## now add a regression model for the location
fit2 <- kappad.fit(y, ydat = cbind(x), mul=1, show=FALSE)
fit2
## now a fit with a fixed shape parameter 
fitf <- kappad.fit(y, show=FALSE, fixedPars = list(sh = 0.2))
fitf ## only three parameters are estimated 
## could also fix the second shape parameter
#fitf2 <- kappad.fit(y, show=FALSE, fixedPars = list(sh2 = -0.4))
#fitf2 ## only three parameters are estimated 
}
\seealso{
\link{dkappa}
}
