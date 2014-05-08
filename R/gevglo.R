#' @title The Generalised Extreme Values distibution
#' @description Density, distribution function, quantile function and random generation for the Generalised Extreme Values distribution with location parameter equal to \code{loc}, scale paraemter equale to \code{scale} and shape parameter equal to \code{shape}
#' @param loc location parameter
#' @param scale scale parameter
#' @param sh shape parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @keywords Generalised Extreme Values distibution
#' @aliases pgev dgev rgev qgev
#' @name Generalised Extreme Values distribution 
#' @details Extra care should be taken on the shape parameter of the distribution as the HW and the Coles' version differ in the sign of the shape parameter: be careful!
#' @return These function mimic the standard output of distributions in R see \code{?pnorm}.  
#' @family gev distribution
#' @examples plot(seq(-15,40,by=0.2),dgev(seq(-15,40,by=0.2),4,6,0.2),type="l")
#' plot(ecdf(rgev(100,4,6,0.2)))
#' lines(seq(-15,40,by=0.5),pgev(seq(-15,40,by=0.5),4,6,0.2),col=2)
#' qgev(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2) # the 1-in-2years, 1-in-100years, 1-in-200years, 1-in-1000years events
#' 
dgev<-function(x,loc,scale,sh)
{
  # Density function
  if(sh> -1e-07 & sh< 1e-07) {tx<-exp(-(x-loc)/scale)}
  else {tx<-(1+(-(x-loc)/scale)*sh)^(1/sh)}
  dgev<-(1/scale)*(tx^(-sh+1))*exp(-tx)
  dgev
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
pgev<-function(q,loc,scale,sh, lower.tail=TRUE){
  # Cummulative distribution function
  if(sh> -1e-07 & sh< 1e-07) {pgev<-exp(-exp(-(q-loc)/scale))}
  else {pgev<-exp(-(1-sh*(q-loc)/scale)^(1/sh))}
  if(!lower.tail) pgev <- 1-pgev
  pgev
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
qgev<-function(p,loc,scale,sh, lower.tail=TRUE){
  # Quantile function
  if(!lower.tail) p <- 1-p
  # Send non-exceedance probability
  if(sh> -1e-07 & sh< 1e-07) {qgev<- loc-scale*log(-log(p))}
  else {qgev<-loc+(scale/sh)*(1-(-log(p))^(sh))}
  qgev
}


#' @name Generalised Extreme Values distribution 
#' @family gev distribution
rgev<-function(n,loc,scale,sh){
  # generate n random variates
  u<-runif(n)
  if(sh> -1e-07 & sh< 1e-07) {rgev<-loc-scale*log(-log(u))}
  else {rgev<-loc+(scale/sh)*(1-(-log(u))^sh)}
  rgev
}



#------glo--------------
#' @title The Generalised Logistic distibution
#' @description Density, distribution function, quantile function and random generation for the Generalised Logistic distribution (as in Hosking and Wallis' book) with location parameter equal to \code{loc}, scale paraemter equale to \code{scale} and shape parameter equal to \code{shape}
#' @param loc location parameter
#' @param scale scale parameter
#' @param sh shape parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @keywords Generalised Logistic distibution
#' @aliases pglo dglo rglo qglo
#' @return These function mimic the standard output of distributions in R see \code{?pnorm}.  
#' @name Generalised Logistic distribution 
#' @export 
#' @family glo distribution
#' @examples plot(seq(-15,40,by=0.2),dglo(seq(-15,40,by=0.2),4,6,0.2),type="l")
#' plot(ecdf(rglo(100,4,6,0.2)))
#' lines(seq(-15,40,by=0.5),pglo(seq(-15,40,by=0.5),4,6,0.2),col=2)
#' qglo(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2) # the 1-in-2years, 1-in-100years, 1-in-200years, 1-in-1000years events
#' 
dglo<-function(x, loc, scale, sh){
  # Denisty function f(x)
  if(sh> -1e-07 & sh< 1e-07) {tx<-(x-loc)/scale}
  else{tx<-(-1/sh)*log(1-sh*(x-loc)/scale)}
  dgev<-(1/scale)*exp(-(1-sh)*tx)/((1+exp(-tx))^2)
  dgev
}

pglo<-function(q, loc, scale, sh, lower.tail=TRUE){
  # Cummulative distribution function F(q)
  if(sh> -1e-07 & sh< 1e-07) {tx<-(x-loc)/scale}
  else {tx<-(-1/sh)*log(1-sh*(q-loc)/scale)}
  pglo<-1/(1+exp(-tx))
  if(!lower.tail) pglo <- 1-pglo
  pglo
}

qglo<-function (p, loc, scale, sh, lower.tail=TRUE){
  ## quantile function F^(-1) (p)
  if(!lower.tail) p <- 1-p
  if(sh> -1e-07 & sh< 1e-07) loc - scale * (log((1 - p)/p))
  else loc + (scale * (1-((1 - p)/p)^(sh)))/sh
}

rglo<-function(n, loc, scale, sh){
  # generate n random variates
  u<-runif(n)
  rglo<-qglo(u, loc, scale, sh)
  rglo
}


