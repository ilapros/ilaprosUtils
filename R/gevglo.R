
#----- the Gumbel distribution -------
#' @title The Gumbel distribution
#' @description Density, distribution function, quantile function and random generation for 
#' the Gumbel distribution with location parameter equal to \code{loc}, 
#' scale parameter equal to \code{scale}. This corresponds to a GEV 
#' distribution with shape parameter equal to 0. 
#' @param loc location parameter
#' @param scale scale parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @keywords Gumbel distribution
#' @aliases pgum dgum rgum qggum
#' @name Gumbel distribution 
#' @return dgum gives the density, pgum gives the distribution function, 
#' qgum gives the quantile function, and rgum generates random deviates.
#' The length of the result is determined by n for rgev, 
#' and is the maximum of the lengths of the numerical arguments for the other functions.
#' The numerical arguments are recycled to the length of the result. 
#' Only the first elements of the logical arguments are used.
#' @family gumbel distribution
#' @export
#' @examples 
#' curve(dgum(x,4,6), from=-15, to = 40,type="l")
#' plot(ecdf(rgum(100,4,6)))
#' lines(seq(-15,40,by=0.5),pgum(seq(-15,40,by=0.5),4,6),col=2)
#' qgum(c(0.5,0.99,0.995,0.995,0.999),4,6) 
#' # notable quantiles
dgum<-function(x, loc, scale, log = FALSE){
  # Density function
  allTogether <- cbind(x,loc,scale)
  allTogether <- cbind(allTogether,rep(0,nrow(allTogether)))
  apply(allTogether,1, dgev_int, log = log)
}

#' @name Gumbel distribution 
#' @family gumbel distribution
#' @export
pgum<-function(q,loc,scale,lower.tail=TRUE, log.p = FALSE){
  allTogether <- cbind(q,loc,scale)
  allTogether <- cbind(allTogether,rep(0,nrow(allTogether)))
  apply(allTogether,1, pgev_int, lower.tail=lower.tail,log.p=log.p)
}

#' @name Gumbel distribution 
#' @family gumbel distribution
#' @export
qgum<-function(p,loc,scale,lower.tail=TRUE,log.p =FALSE){
  allTogether <- cbind(p,loc,scale)
  allTogether <- cbind(allTogether,rep(0,nrow(allTogether)))
  apply(allTogether,1, qgev_int, lower.tail=lower.tail,log.p =log.p)
}

#' @name Gumbel distribution 
#' @family gumbel distribution
#' @export
rgum<-function(n,loc,scale){
  # generate n random variates
  allTogether <- cbind(rep(loc, length.out=n), 
                       rep(scale, length.out=n),
                       rep(0, length.out=n))
  ## add random percentiles
  allTogether <- cbind(runif(n), allTogether)
  apply(allTogether, 1, qgev_int, lower.tail=FALSE)
}



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
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p)
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
dgev<-function(x, loc, scale, sh, log = FALSE){
  # Density function
  allTogether <- cbind(x,loc,scale,sh)
  apply(allTogether,1, dgev_int, log = log)
}

dgev_int<-function(x, log = FALSE){
  loc = x[2]; scale= x[3]; sh= x[4]; x = x[1]
  # Density function
  if(sh> -1e-07 & sh< 1e-07) {tx<-exp(-(x-loc)/scale)}
  else {tx<-(1+(-(x-loc)/scale)*sh)^(1/sh)}
  dgev<-pmax((1/scale)*(tx^(-sh+1))*exp(-tx),0, na.rm = TRUE)
  if(log) dgev <- log(dgev)
  dgev
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
#' @export
pgev<-function(q,loc,scale,sh, lower.tail=TRUE, log.p = FALSE){
    allTogether <- cbind(q,loc,scale,sh)
    apply(allTogether,1, pgev_int, lower.tail=lower.tail,log.p=log.p)
}
  
pgev_int<-function(x,lower.tail=lower.tail,log.p=log.p){
  q = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Cumulative distribution function
  if(sh> -1e-07 & sh< 1e-07) {pgev<-exp(-exp(-(q-loc)/scale))}
  else {
    ## if sh < 0 upper bound is loc+scale/sh
    ## if sh > 0 lower bound is loc+scale/sh
    q <- ifelse(sh < -10^-7, 
                max(q, loc+scale/sh), 
                min(q, loc+scale/sh))
    pgev<-exp(-(1-sh*(q-loc)/scale)^(1/sh))
    }
  if(lower.tail) return(ifelse(log.p, log(pgev), pgev))
  else return(ifelse(log.p, log(1-pgev), 1-pgev))
}

#' @name Generalised Extreme Values distribution 
#' @family gev distribution
#' @export
qgev<-function(p,loc,scale,sh, lower.tail=TRUE,log.p =FALSE){
  allTogether <- cbind(p,loc,scale,sh)
  apply(allTogether,1, qgev_int, lower.tail=lower.tail,log.p =log.p)
}
  
qgev_int<-function(x,lower.tail=TRUE,log.p =FALSE){
  p = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  p <- ifelse(log.p,exp(p),p)
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
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p)
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

dglo<-function(x, loc, scale, sh,log=FALSE){
  # Density function
  allTogether <- cbind(x,loc,scale,sh)
  apply(allTogether,1, dglo_int,log=log)
}

dglo_int<-function(x,log){
  loc = x[2]; scale= x[3]; sh= x[4]; x = x[1]
  # Density function
  if(sh> -1e-07 & sh< 1e-07) {tx <- (x-loc)/scale}
  else{tx <- (-1/sh)*log(pmax(1-sh*(x-loc)/scale,0))}
  dglo <- pmax((1/scale)*exp(-(1-sh)*tx)/((1+exp(-tx))^2), 0, na.rm = TRUE)
  dglo <- ifelse(log,log(dglo),dglo)
  dglo
}



#' @name Generalised Logistic distribution
#' @family glo distribution
#' @export
pglo<-function(q,loc,scale,sh, lower.tail=TRUE,log.p=FALSE){
  allTogether <- cbind(q,loc,scale,sh)
  apply(allTogether,1, pglo_int, lower.tail=lower.tail,log.p=log.p)
}

pglo_int<-function(x,lower.tail=lower.tail,log.p=FALSE){
  q = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  # Cumulative distribution function
  if(sh> -1e-07 & sh< 1e-07) {tx<-(x-loc)/scale}
  else {
    q <- ifelse(sh < -10^-7, 
                max(q, loc+scale/sh), 
                min(q, loc+scale/sh))
    tx<-(-1/sh)*log(1-sh*(q-loc)/scale)
    }
  pglo<-1/(1+exp(-tx))
  if(lower.tail) return(ifelse(log.p, log(pglo), pglo))
  else return(ifelse(log.p, log(1-pglo), 1-pglo))
}



#' @name Generalised Logistic distribution
#' @family glo distribution
#' @export
qglo<-function(p,loc,scale,sh, lower.tail=TRUE,log.p=FALSE){
  allTogether <- cbind(p,loc,scale,sh)
  apply(allTogether,1, qglo_int, lower.tail=lower.tail,log.p=log.p)
}

qglo_int<-function(x,lower.tail=TRUE,log.p=FALSE){
  p = x[1];loc = x[2]; scale= x[3]; sh= x[4] 
  p <- ifelse(log.p,exp(p),p)
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




