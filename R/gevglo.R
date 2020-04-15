
#----- the GEV distributon -------
#' @title The Generalised Extreme Values distribution
#' @description Density, distribution function, quantile function and random generation for 
#' the Generalised Extreme Values distribution with location parameter equal to \code{loc}, 
#' scale parameter equal to \code{scale} and shape parameter equal to \code{sh}. 
#' The functions use the Hosking and Wallis notation so that the domain of the distribution 
#' has an upper bound when the shape parameter is positive. 
#' @param loc location parameter
#' @param scale scale parameter
#' @param sh shape parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @keywords Generalised Extreme Values distribution
#' @aliases pgev dgev rgev qgev
#' @name Generalised Extreme Values distribution 
#' @details Extra care should be taken on the shape parameter of the distribution.
#' Different notations are used in the scientific literature: in one notation, presented 
#' for example in the Hosking and Wallis book and common in hydrology, 
#' the domain of the distribution has an upper bound when the shape parameter is positive. 
#' Conversely, in one notation, presented for example in Coles' book and Wikipedia and 
#' common in statistics, the domain of the distribution has a lower bound when the shape parameter is positive. 
#' The two notation only differ for the sign of the shape parameter. 
#' The functions in this package use the Hosking and Wallis notation for consistency with the \code{\link{pglo}} functions. 
#' Nevertheless the fitting in the \code{gev.fit} function is based on the \code{isemv::gev.fit} function written by Stuart Coles which uses the 
#' Coles' notation: be aware of these differences! 
#' @return dgev gives the density, pgev gives the distribution function, qgev gives the quantile function, and rgev generates random deviates.
#' The length of the result is determined by n for rgev, and is the maximum of the lengths of the numerical arguments for the other functions.
#' The numerical arguments are recycled to the length of the result. 
#' Only the first elements of the logical arguments are used.
#' @family gev distribution
#' @export
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' @examples 
#' plot(seq(-15,40,by=0.2),dgev(seq(-15,40,by=0.2),4,6,0.2),type="l")
#' plot(ecdf(rgev(100,4,6,0.2)))
#' lines(seq(-15,40,by=0.5),pgev(seq(-15,40,by=0.5),4,6,0.2),col=2)
#' qgev(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2) 
#' # notable quantiles
dgev<-function(x, loc, scale, sh){
  # Density function
  allTogether <- cbind(x,loc,scale,sh)
  apply(allTogether,1, dgev_int)
}

dgev_int<-function(x){
  loc = x[2]; scale= x[3]; sh= x[4]; x = x[1]
  # Density function
  if(sh> -1e-07 & sh< 1e-07) {tx<-exp(-(x-loc)/scale)}
  else {tx<-(1+(-(x-loc)/scale)*sh)^(1/sh)}
  dgev<-(1/scale)*(tx^(-sh+1))*exp(-tx)
  dgev
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
#' @export
pgev<-function(q,loc,scale,sh, lower.tail=TRUE){
    allTogether <- cbind(q,loc,scale,sh)
    apply(allTogether,1, pgev_int, lower.tail=lower.tail)
}
  
pgev_int<-function(x,lower.tail=lower.tail){
  q = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Cummulative distribution function
  if(sh> -1e-07 & sh< 1e-07) {pgev<-exp(-exp(-(q-loc)/scale))}
  else {pgev<-exp(-(1-sh*(q-loc)/scale)^(1/sh))}
  if(!lower.tail) pgev <- 1-pgev
  pgev
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
#' @export
qgev<-function(p,loc,scale,sh, lower.tail=TRUE){
  allTogether <- cbind(p,loc,scale,sh)
  apply(allTogether,1, qgev_int, lower.tail=lower.tail)
}
  
qgev_int<-function(x,lower.tail=lower.tail){
  p = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Quantile function
  if(!lower.tail) p <- 1-p
  # Send non-exceedance probability
  if(sh> -1e-07 & sh< 1e-07) {qgev<- loc-scale*log(-log(p))}
  else {qgev<-loc+(scale/sh)*(1-(-log(p))^(sh))}
  qgev
}


#' @name Generalised Extreme Values distribution 
#' @family gev distribution
#' @export
rgev<-function(n,loc,scale,sh){
  # generate n random variates
  allTogether <- cbind(rep(loc, length.out=n), 
                       rep(scale, length.out=n),
                       rep(sh, length.out=n))
  ## add random percentiles
  allTogether <- cbind(runif(n), allTogether)
  apply(allTogether, 1, qgev_int, lower.tail=FALSE)
#  apply(allTogether, 1, rgev_int)
}

# rgev_int<-function(x){
#   loc = x[1]; scale= x[2]; sh= x[3]
#   qgev(p = runif(1), loc = loc, scale = scale, sh = sh)
# }



#------glo--------------
#' @title The Generalised Logistic distribution
#' @description Density, distribution function, quantile function and random generation for the
#' Generalised Logistic distribution (as in Hosking and Wallis' book) with location parameter 
#' equal to \code{loc}, scale parameter equal to \code{scale} and shape parameter equal to \code{sh}
#' @param loc location parameter
#' @param scale scale parameter
#' @param sh shape parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @return dglo gives the density, pglo gives the distribution function, qglo gives the quantile function, 
#' and rglo generates random deviates. The length of the result is determined by n for rglo, 
#' and is the maximum of the lengths of the numerical arguments for the other functions.
#' The numerical arguments are recycled to the length of the result. 
#' Only the first elements of the logical arguments are used.
#' @export
#' @name Generalised Logistic distribution
#' @family glo distribution
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' @examples 
#' plot(seq(-26,80,by=0.2),dglo(seq(-26,80,by=0.2),4,6,-0.2),type="l")
#' plot(ecdf(rglo(100,4,6,-0.2)))
#' lines(seq(-26,80,by=0.2),pglo(seq(-26,80,by=0.2),4,6,-0.2),col=2)
#' qglo(c(0.5,0.99,0.995,0.995,0.999),4,6,-0.2) 
#' # notable quantiles

dglo<-function(x, loc, scale, sh){
  # Density function
  allTogether <- cbind(x,loc,scale,sh)
  apply(allTogether,1, dglo_int)
}

dglo_int<-function(x){
  loc = x[2]; scale= x[3]; sh= x[4]; x = x[1]
  # Density function
  if(sh> -1e-07 & sh< 1e-07) {tx <- (x-loc)/scale}
  else{tx <- (-1/sh)*log(1-sh*(x-loc)/scale)}
  dglo <- (1/scale)*exp(-(1-sh)*tx)/((1+exp(-tx))^2)
  dglo
}



#' @name Generalised Logistic distribution
#' @family glo distribution
#' @export
pglo<-function(q,loc,scale,sh, lower.tail=TRUE){
  allTogether <- cbind(q,loc,scale,sh)
  apply(allTogether,1, pglo_int, lower.tail=lower.tail)
}

pglo_int<-function(x,lower.tail=lower.tail){
  q = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Cummulative distribution function
  if(sh> -1e-07 & sh< 1e-07) {tx<-(x-loc)/scale}
  else {tx<-(-1/sh)*log(1-sh*(q-loc)/scale)}
  pglo<-1/(1+exp(-tx))
  if(!lower.tail) pglo <- 1-pglo
  pglo
}



#' @name Generalised Logistic distribution
#' @family glo distribution
#' @export
qglo<-function(p,loc,scale,sh, lower.tail=TRUE){
  allTogether <- cbind(p,loc,scale,sh)
  apply(allTogether,1, qglo_int, lower.tail=lower.tail)
}

qglo_int<-function(x,lower.tail=lower.tail){
  p = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Quantile function
  if(!lower.tail) p <- 1-p
  if(sh> -1e-07 & sh< 1e-07) loc - scale * (log((1 - p)/p))
  else loc + (scale * (1-((1 - p)/p)^(sh)))/sh
}


#' @name Generalised Logistic distribution
#' @family glo distribution
#' @export
rglo<-function(n, loc, scale, sh){
    # generate n random variates
    ### the parameters
    allTogether <- cbind(rep(loc, length.out=n), 
                     rep(scale, length.out=n),
                     rep(sh, length.out=n))
    ## add random percentiles
    allTogether <- cbind(runif(n), allTogether)
    apply(allTogether, 1, qglo_int, lower.tail=FALSE)
}

# rglo_int<-function(x){
#   loc = x[1]; scale= x[2]; sh= x[3]
#   qglo(p = runif(1), loc = loc, scale = scale, sh = sh)
# }




