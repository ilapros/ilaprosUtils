% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ismevExtension.R
\name{glo.rl}
\alias{glo.rl}
\title{Fitted return curve for a glo.fit object}
\usage{
glo.rl(
  a,
  mat,
  dat,
  nh = 0,
  nk = 0,
  X0 = NULL,
  f = c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
    0.99, 0.995, 0.999)
)
}
\arguments{
\item{a}{the mle estimates from a \code{glo.fit} object}

\item{mat}{the covariance function from a \code{glo.fit} object}

\item{dat}{data matrix from a \code{glo.fit} object}

\item{nh}{if historical data are used, the length of historical record - default is 0, no historical data}

\item{nk}{if historical data are used, the number of peaks above X0 - default is 0, no historical data}

\item{X0}{if historical data are used, the perception threshold - default is NULL, no historical data}

\item{f}{frequencies for which return level (and 95\%) confidence intervals are calculated}
}
\description{
This function mimics the \code{gev.rl} function for GLO models.
It plots the flood frequency curve (return curve) based on the data and the fitted parameters for a \code{glo.fit} model, including ones with historical data. 
Very ad-hoc and working under assumption that flood data are plotted (hence the y-axis lab).
Also outputs useful informations. 
Mostly written by Thomas Kjeldsen. 
This function is deprecated and has been replaced by the generic retPlot function.
}
\examples{
set.seed(7821567)
xx <- rglo(500, 40, 6, -0.2)
xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
h1 <- glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, show=FALSE)
s1 <- glod.fit(xxsist, show=FALSE) 
rls1 <- glo.rl(a=s1$mle,mat=s1$cov,dat=s1$data)
rlh1 <- glo.rl(h1$mle,h1$cov,h1$data,nh=h1$h,nk=h1$k,X0=h1$X0)
lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl+1.96*sqrt(rlh1$var), lty = 2)
lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl-1.96*sqrt(rlh1$var), lty = 2)
lines(log(rls1$f/(1-rls1$f)), rls1$rl, col = 2)
lines(log(rls1$f/(1-rls1$f)), rls1$rl-1.96*sqrt(rls1$var), lty = 2, col = 2)
lines(log(rls1$f/(1-rls1$f)), rls1$rl+1.96*sqrt(rls1$var), lty = 2, col = 2)
legend("topleft", col =c(1,2),legend = c("With historical","Systematic Only"), lty = 1)
## similar fitted curve - but obvious reduction in uncertainty 
}
\seealso{
\code{\link{retPlot}}
}
