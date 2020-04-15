
######### pretty printing ----- 

#' Nicer print of gev.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class gev.fit.
#' @param ... further arguments passed to \code{print}.
#' @seealso gevd.fit, isemv::gev.fit 
#' @examples
#' # library(ismev)
#' print(ismev::gev.fit(c(50,45,65,78,12,23),show=FALSE))
#' print(gevd.fit(c(50,45,65,78,12,23),show=FALSE))
#' @export 
#' @import ismev
#' @import stats
print.gev.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}


#' Nicer print of gevcv.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class gevcv.fit.
#' @param ... further arguments passed to \code{print}.
#' @examples
#' print(gevcvd.fit(c(50,45,65,78,12,23),show=FALSE))
#' @seealso gevcvd.fit, isemv::gev.fit 
#' @export
print.gevcv.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}

#' Nicer print of glo.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class glo.fit. 
#' @param ... further arguments passed to \code{print}.
#' @seealso glod.fit
#' @examples
#' library(ismev)
#' print(glod.fit(c(50,45,65,78,12,23),show=FALSE))
#' @export
print.glo.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}




#' Nicer print of glocv.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class gevcv.fit. 
#' @param ... further arguments passed to \code{print}.
#' @examples
#' print(glocvd.fit(c(50,45,65,78,12,23),show=FALSE))
#' @seealso glocvd.fit, glod.fit 
#' @export
print.glocv.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}


#' Nicer print of gpd.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class pp.fit. 
#' @param ... further arguments passed to \code{print}.
#' @keywords gpd.fit
#' @examples
#' y <- c(53, 52, 49, 58, 50, 48, 47, 50, 46, 46, 49, 51, 47, 49, 50) 
#' a <- ismev::gpd.fit(y, threshold = 46, show=FALSE)
#' a
#' @export
print.gpd.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}


#' nicer print of pp.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class pp.fit. 
#' @param ... further arguments passed to \code{print}.
#' @keywords pp.fit
#' @examples
#' data(rain, package = "ismev")
#' a <- ismev::pp.fit(rain, 10)
#' a
#' @export
print.pp.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}


######### GEV fitting ----- 

#' @title Maximum Likelihood Fitting for the GEV distribution
#' @description Maximum-likelihood fitting for the generalized extreme value distribution, 
#' including generalized linear modelling of each parameter. 
#' This function has the same structure as the \code{gev.fit} from \code{ismev}. 
#' Additionally it allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param sigl As \code{mul} for the scale parameter
#' @param shl As \code{mul} for the shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param show  Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details).
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations 
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the gev.fit class as for objects obtained using ismev::gev.fit, so that \code{ismev} functions
#' for gev.fit work on these objects when no parameter is kept fixed (i.e. \code{fixedPars} is a list of NULL)   
#' @seealso \link{dgev}
#' @details The form of the GEV used is that of Coles (2001) Eq (3.2) - so it is not the same as the one used in \link{dgev}. 
#' Specifically, positive values of the shape parameter imply a heavy tail, and negative values imply a bounded upper tail.
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' 
#' Coles, S., 2001. An introduction to statistical modeling of extreme values. London: Springer.
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(500)
#' y <- rgev(500,loc = 40+4*x,scale = 6,sh = 0.2)
#' fit1 <- gevd.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- gevd.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed shape parameter (notice the sign)
#' fitf <- gevd.fit(y, show=FALSE, fixedPars = list(sh = -0.2))
#' fitf ## only two parameters are estimated 
gevd.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 
                      mulink = identity, siglink = identity, shlink = identity, 
                      muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
                      method = "Nelder-Mead", optimPars = NULL, maxit = 10000, 
                      fixedPars = list(mu = NULL, sig = NULL, sh = NULL), ...) {
  z <- list()
  npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
  npsc <- ifelse(is.null(fixedPars[["sig"]]), length(sigl) + 1, 0)
  npsh <- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  
  ### include in inital values only the parameters for whcih estiamtion is required 
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  if(is.null(fixedPars[["mu"]])) muinit <- muinit
  if(!is.null(fixedPars[["mu"]])) muinit <- NULL
  if(is.null(fixedPars[["sig"]])) siginit <- siginit
  if(!is.null(fixedPars[["sig"]])) siginit <-NULL
  if(is.null(fixedPars[["sh"]])) shinit <- shinit
  if(!is.null(fixedPars[["sh"]])) shinit <- NULL
  init <- c(muinit, siginit, shinit)
  gev.lik <- function(a) {
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu])) 
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["sig"]])) sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
    if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
  }
  x <- eval(parse(text = paste0(
    "optim(init, gev.lik, ",
    optimPars,
    ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
    # x <- optim(init, gev.lik, hessian = TRUE, method = method, 
    #          control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu])) 
  if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
  if(is.null(fixedPars[["sig"]])) sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
  if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
  # mu <- mulink(mumat %*% (x$par[1:npmu]))
  # sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  # xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  vals <- cbind(mu, sc, xi)
  colnames(vals) <- c("location","scale","shape")
  z$vals <- vals 
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "gev.fit"
  invisible(z)
}



#' @title Maximum Likelihood Fitting for the GEV distribution - CV model
#' @description Maximum-likelihood fitting for the generalized extreme value distribution, 
#' including generalized linear modelling of each parameter. The function differs from 
#' \code{gevd.fit} because it uses a different parametrisation of the distribution based on 
#' the \eqn{\tau}, the ratio of the scale parameter and the location parameter, which is a monotonic 
#' function of the coefficient of variation. This means that when regression models are applied for the location,
#' these also affect the scale. 
#' The function allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param taul As \code{mul} for the tau parameter
#' @param shl As \code{mul} for the shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param taulink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param tauinit initial values for the tau parameter
#' @param shinit initial values for the shape parameter
#' @param show  Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details)
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations 
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the \code{gevcd.fit} class, similar to \code{gev.fit} objects obtained using ismev::gev.fit. 
#' The \code{ismev} functions for gev.fit will not work on these objects. 
#' 
#' In the output the \code{vals} matrix gives the location and scale values obtained as scale = \eqn{\tau} *  location.   
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(500)
#' y <- rgev(500,loc = 40+4*x, scale = 0.2*(40+4*x), sh = 0.15)
#' fit1 <- gevcvd.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- gevcvd.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed tau parameter
#' fitf <- gevcvd.fit(y, ydat = cbind(x), mul=1, show=FALSE, fixedPars = list(tau = 0.2))
#' fitf ## only two parameters are estimated (location and shape)
gevcvd.fit <- function (xdat, ydat = NULL, mul = NULL, taul = NULL, shl = NULL, 
                         mulink = identity, taulink = identity, shlink = identity, 
                         muinit = NULL, tauinit = NULL, shinit = NULL, show = TRUE, 
                         method = "Nelder-Mead", optimPars = NULL, maxit = 10000, 
                         fixedPars = list(mu = NULL, tau = NULL, sh = NULL),...) {
  z <- list()
  npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
  npcv <- ifelse(is.null(fixedPars[["tau"]]), length(taul) + 1, 0)
  npsh <- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
  in2 <- in2/in1
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(taul)) {
    taumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(tauinit)) 
      tauinit <- in2
  }
  else {
    z$trans <- TRUE
    taumat <- cbind(rep(1, length(xdat)), ydat[, taul])
    if (is.null(tauinit)) 
      tauinit <- c(in2, rep(0, length(taul)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  z$model <- list(mul, taul, shl)
  z$link <- deparse(substitute(c(mulink, taulink, shlink)))
  if(is.null(fixedPars[["mu"]])) muinit <- muinit
  if(!is.null(fixedPars[["mu"]])) muinit <- NULL
  if(is.null(fixedPars[["tau"]])) tauinit <- tauinit
  if(!is.null(fixedPars[["tau"]])) tauinit <-NULL
  if(is.null(fixedPars[["sh"]])) shinit <- shinit
  if(!is.null(fixedPars[["sh"]])) shinit <- NULL
  init <- c(muinit, tauinit, shinit)
  gev.lik <- function(a) {
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu])) 
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (a[seq(npmu + 1, length = npcv)]))
    if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
    if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (a[seq(npmu + npcv + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
    sc <- cv*mu
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 
                                                    1))
  }
  x <- eval(parse(text = paste0(
    "optim(init, gev.lik, ",
    optimPars,
    ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
  # x <- optim(init, gev.lik, hessian = TRUE, method = method, 
  #            control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu])) 
  if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
  if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (x$par[seq(npmu + 1, length = npcv)]))
  if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
  if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (x$par[seq(npmu + npcv + 1, length = npsh)]))
  if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
  sc <- mu * cv
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  vals <- cbind(mu, cv*mu, xi)
  colnames(vals) <- c("location","scale","shape")
  z$vals <- vals 
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "gevcv.fit"
  invisible(z)
}



######### GLO fitting ----- 

#' @title Maximum Likelihood Fitting for the GLO distribution
#' @description Maximum-likelihood fitting for the generalized logistic distribution, 
#' including generalized linear modelling of each parameter. 
#' This function has the same structure as the \code{gevd.fit} and is inspired by \code{ismev::gev.fit}. 
#' The function allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param sigl As \code{mul} for the scale parameter
#' @param shl As \code{mul} for the shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param show  Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details).
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations 
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the glo.fit class - with values which mirror the ones of the gev.fit class in ismev
#' @details The distribution is discussed in the Hosking and Wallis book and is used as the default distribution
#' for flood frequency estimation in the UK
#' @seealso \link{dglo}
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' 
#' Coles, S., 2001. An introduction to statistical modeling of extreme values. London: Springer.
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(500)
#' y <- rglo(500,loc = 40+4*x,scale = 6,sh = 0.2)
#' fit1 <- glod.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- glod.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed shape parameter 
#' fitf <- glod.fit(y, show=FALSE, fixedPars = list(sh = 0.2))
#' fitf ## only two parameters are estimated 
glod.fit<-function(xdat, ydat = NULL, 
                  mul = NULL, sigl = NULL, shl = NULL, 
                  mulink = identity, siglink = identity, shlink = identity, 
                  muinit = NULL, siginit = NULL, shinit = NULL, 
                  show = TRUE, method = "Nelder-Mead", optimPars = NULL, maxit = 10000, 
                  fixedPars = list(mu = NULL, sig = NULL, sh = NULL), ...) {
    #
    # obtains mles etc for glo distn
    #
    z <- list()
    npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
    npsc <- ifelse(is.null(fixedPars[["sig"]]), length(sigl) + 1, 0)
    npsh <- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
    z$trans <- FALSE
    in2 <- sqrt(6 * var(xdat))/pi
    in1 <- mean(xdat) - 0.57722 * in2
    if (is.null(mul)) {
      mumat <- as.matrix(rep(1, length(xdat)))
      if (is.null(muinit)) 
        muinit <- in1
    }
    else {
      z$trans <- TRUE
      mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
      if (is.null(muinit)) 
        muinit <- c(in1, rep(0, length(mul)))
    }
    if (is.null(sigl)) {
      sigmat <- as.matrix(rep(1, length(xdat)))
      if (is.null(siginit)) 
        siginit <- in2
    }
    else {
      z$trans <- TRUE
      sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
      if (is.null(siginit)) 
        siginit <- c(in2, rep(0, length(sigl)))
    }
    if (is.null(shl)) {
      shmat <- as.matrix(rep(1, length(xdat)))
      if (is.null(shinit)) 
        shinit <- 0.1
    }
    else {
      z$trans <- TRUE
      shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
      if (is.null(shinit)) 
        shinit <- c(0.1, rep(0, length(shl)))
    }
    
    ### include in inital values only the parameters for which estiamtion is required 
    z$model <- list(mul, sigl, shl)
    z$link <- deparse(substitute(c(mulink, siglink, shlink)))
    if(is.null(fixedPars[["mu"]])) muinit <- muinit
    if(!is.null(fixedPars[["mu"]])) muinit <- NULL
    if(is.null(fixedPars[["sig"]])) siginit <- siginit
    if(!is.null(fixedPars[["sig"]])) siginit <-NULL
    if(is.null(fixedPars[["sh"]])) shinit <- shinit
    if(!is.null(fixedPars[["sh"]])) shinit <- NULL
    init <- c(muinit, siginit, shinit)
    glo.lik <- function(a) {
      # computes neg log lik of glo model
      if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu])) 
      if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
      if(is.null(fixedPars[["sig"]])) sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
      if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
      if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
      if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
      y <- ((1-xi*(xdat - mu)/sc))
      if(any(!y>0) || any(sc <= 0)) return(10^6)
      y<- -log(y)/xi
      if(any(abs(xi) <10^-5)) {y <- (xdat - mu[1])/sc[1]}
      (sum(log(sc)) + sum(y*(1-xi)) + 2*sum(log(1+exp(-y) )))
    }
    x <- eval(parse(text = paste0(
      "optim(init, glo.lik, ",
      optimPars,
      ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
    # x <- optim(init, glo.lik, hessian = TRUE, method = method, 
    #            control = list(maxit = maxit, ...))
    z$conv <- x$convergence
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu])) 
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["sig"]])) sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
    if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
    z$nllh <- x$value
    z$data <- xdat
    if(z$trans) {
      z$data <- - log(as.vector((1 - (xi * (xdat - mu))/sc)^(1/xi)))
    }
    z$mle <- x$par
    z$cov <- solve(x$hessian)
    z$se <- sqrt(diag(z$cov))
    vals <- cbind(mu, sc, xi)
    colnames(vals) <- c("location","scale","shape")
    z$vals <- vals 
    if(show) {
      if(z$trans)
        print(z[c(2, 3, 4)])
      else print(z[4])
      if(!z$conv)
        print(z[c(5, 7, 9)])
    }
    class(z) <- "glo.fit"
    invisible(z)
}




#' @title Maximum Likelihood Fitting for the GLO distribution - CV model
#' @description Maximum-likelihood fitting for the generalized logistic distribution, 
#' including generalized linear modelling of each parameter. The function differs from 
#' \code{glod.fit} because it uses a different parametrisation of the distribution based on 
#' the \eqn{\tau}, the ratio of the scale parameter and the location parameter, which is a monotonic 
#' function of the coefficient of variation. This means that when regression models are applied for the location,
#' these also affect the scale. 
#' The function allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param taul As \code{mul} for the tau parameter
#' @param shl As \code{mul} for the shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param taulink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param tauinit initial values for the tau parameter
#' @param shinit initial values for the shape parameter
#' @param show  Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details).
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations 
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the \code{glocd.fit} class, similar to \code{glo.fit}. 
#' #' 
#' In the output the \code{vals} matrix gives the location and scale values obtained as scale = \eqn{\tau} *  location. 
#' @details The distribution is discussed in the Hosking and Wallis book and is used as the default distribution
#' for flood frequency estimation in the UK
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(500)
#' y <- rglo(500,loc = 40+4*x, scale = 0.2*(40+4*x), sh = 0.15)
#' fit1 <- glocvd.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- glocvd.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed tau parameter
#' fitf <- glocvd.fit(y, ydat = cbind(x), mul=1, show=FALSE, fixedPars = list(tau = 0.2))
#' fitf ## only two parameters are estimated (location and shape)
glocvd.fit <- function (xdat, ydat = NULL, mul = NULL, taul = NULL, shl = NULL, 
                        mulink = identity, taulink = identity, shlink = identity, 
                        muinit = NULL, tauinit = NULL, shinit = NULL, show = TRUE, 
                        method = "Nelder-Mead", optimPars = NULL, maxit = 10000, 
                        fixedPars = list(mu = NULL, tau = NULL, sh = NULL),...) {
  z <- list()
  npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
  npcv <- ifelse(is.null(fixedPars[["tau"]]), length(taul) + 1, 0)
  npsh <- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
  in2 <- in2/in1
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(taul)) {
    taumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(tauinit)) 
      tauinit <- in2
  }
  else {
    z$trans <- TRUE
    taumat <- cbind(rep(1, length(xdat)), ydat[, taul])
    if (is.null(tauinit)) 
      tauinit <- c(in2, rep(0, length(taul)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  z$model <- list(mul, taul, shl)
  z$link <- deparse(substitute(c(mulink, taulink, shlink)))
  if(is.null(fixedPars[["mu"]])) muinit <- muinit
  if(!is.null(fixedPars[["mu"]])) muinit <- NULL
  if(is.null(fixedPars[["tau"]])) tauinit <- tauinit
  if(!is.null(fixedPars[["tau"]])) tauinit <-NULL
  if(is.null(fixedPars[["sh"]])) shinit <- shinit
  if(!is.null(fixedPars[["sh"]])) shinit <- NULL
  init <- c(muinit, tauinit, shinit)
  glo.lik <- function(a) {
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu])) 
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (a[seq(npmu + 1, length = npcv)]))
    if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
    if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (a[seq(npmu + npcv + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
    sc <- cv*mu
    y <- ((1-xi*(xdat - mu)/sc))
    if(any(!y>0) || any(sc <= 0)) return(10^6)
    y<- -log(y)/xi
    if(any(abs(xi) <10^-5)) {y <- (xdat - mu[1])/sc[1]}
    (sum(log(sc)) + sum(y*(1-xi)) + 2*sum(log(1+exp(-y) )))
  }
  x <- eval(parse(text = paste0(
    "optim(init, glo.lik, ",
    optimPars,
    ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
  # x <- optim(init, glo.lik, hessian = TRUE, method = method, 
  #            control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu])) 
  if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
  if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (x$par[seq(npmu + 1, length = npcv)]))
  if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
  if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (x$par[seq(npmu + npcv + 1, length = npsh)]))
  if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
  sc <- mu * cv
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- - log(as.vector((1 - (xi * (xdat - mu))/sc)^(1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  vals <- cbind(mu, cv*mu, xi)
  colnames(vals) <- c("location","scale","shape")
  z$vals <- vals 
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "glocv.fit"
  invisible(z)
}


####### corollary functions ---- 
#' @import graphics 

gumbX <- function(x, emp=FALSE) {
  if(emp) x <- (1:length(x)-0.44)/(length(x) + 1 - 0.88)
  -log(-log(x))
}

redVarX <- function(x, emp=FALSE) {
  if(emp) x <- (1:length(x)-0.44)/(length(x) + 1 - 0.88)
  log((x)/(1-x))
}


#' @title Plots of return curves and confidence intervals estimated by the delta method
#' @description This function is based on the \code{ismev::gev.rl} function but also allows the case in which historical data are added to the systematic record.
#' It plots the return curve (flood frequency curve in hydrology) based on the data and the fitted parameters for a \code{glo.fit} or \code{gev.fit} object, including ones with historical data. 
#' The data points are plotted using the Gringorten plotting positions. 
#' The function also outputs useful the estimated return levels and the corresponding standard error as estimated via the delta method (see Coles, 2001). 
#' @param obj a \code{glo.fit} or \code{gev.fit} object
#' @param p non-exceedance probabilities for which return level and standard errors are calculated. If the given non-exceedance probability don't cover the whole range of empirical non-exceedance probabilities of the data the figure will automatically draw a line covering the whole range, but the output dataframe will only contain the specified non-exceedance probabilities.
#' @param vtype a character specifying the reduced variate type to be used in the figure. The types allowed are "Gumbel", corresponding to -log(-log(p)), and "redVar", corresponding to log(p/(1-p)). The default is set to "redVar".
#' @param sign.alpha significance level required for the confidence intervals - default to 0.05
#' @param pch pch parameter to be used for the (systematic) data - default to the current setting in \code{par}
#' @param pchHist pch parameter to be used for the historical data (if present) - default to 15
#' @param plot.out logical, indicating whether the plot should actually be displayed; set to FALSE to only compute the return levels and standard errors
#' @param ...	Arguments to be passed to methods, such as graphical parameters (see par)
#' @return a return levels figure and a data.frame containing the estimated return levels and corresponding standard errors for the specified exceedance probabilities
#' @export
#' @examples 
#' set.seed(7821567)
#' xx <- rglo(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' s1 <- glod.fit(xxsist, show=FALSE) 
#' rls1 <- retPlot(s1, sign.alpha = 0.1, col = 4, 
#'    p = c(seq(0.01,0.991,length=45),seq(0.992,0.9992,length=120)))
#' h1 <- glo.hist.fit(c(xxhist,xxsist), 
#'    k = length(xxhist), h = 470, X0 = 80, show=FALSE)
#' rlh1 <- retPlot(h1, vtype = "Gumbel", col = 1, 
#'    sign.alpha = 0.05, p = rls1$p, 
#'    xlab = "Gumbel reduced variate (-log(-log(1-1/T)))")
#' lines(-log(-log(rls1$p)), rls1$retLev, col = 2)
#' lines(-log(-log(rls1$p)), 
#'    rls1$retLev+qnorm(0.025)*rls1$se, lty = 2, col = 2)
#' lines(-log(-log(rls1$p)), 
#'    rls1$retLev-qnorm(0.025)*rls1$se, lty = 2, col = 2)
#' legend("topleft", col =c(1,2),
#'    legend = c("With historical","Systematic Only"), lty = 1)
#' ## similar fitted curve - but large reduction in uncertainty for rare events
retPlot <- function(obj, p=NULL, 
                    sign.alpha = 0.05,
                    plot.out = TRUE,
                    vtype = "redVar", 
                    pch = par()$pch, pchHist = 15, ...) { # , gridProb = NULL, gridLab,
  if(!(class(obj) %in% c("glo.fit","gev.fit"))) stop("Function works only for gev.fit and glo.fit objects")
  if(!(length(as.vector(obj$mle)) == 3)) stop("The function assumes a model with three estimated parameters (location, scale, shape)")
  a <- as.vector(obj$mle); mat <- as.matrix(obj$cov); dat <- as.vector(obj$data)
  nh <- 0; nk <- 0; X0 <- 0
  if(exists("X0",obj)){
    nh <- obj$h; nk <- obj$k; X0 <- obj$X0
  }
  a <- a + 1e-006
  if(is.null(p)) p <- c(0.01, 0.02, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 0.95, 0.99, 0.991, seq(0.992, 0.995, by = 0.001))
  f <- 1-p
  fline <- f
  ### function to calculate the right type of reduced variate
  xcalc <- function(x, vtype, emp = FALSE) switch(vtype,"Gumbel" = gumbX(x, emp = emp),"redVar" = redVarX(x, emp = emp))
  
  xemp <- xcalc(c(dat,rep(mean(dat),nh)), vtype = vtype, emp = TRUE)
  if(max(xcalc(1-f, vtype = vtype)) < max(xemp)) f <- sort(c(f, seq(1/((length(dat)+nh)), 1/(2*(length(dat)+nh)),l=120)))
  xv <- xcalc(1-f, vtype = vtype)

  ## quantiles
  y <- switch(class(obj),
    "glo.fit" = gloq(f,loc = a[1],scale = a[2],sh = a[3]),
    "gev.fit" = ismev::gevq(p = f, a = a))

  # Confidence intervals by means of delta method
  v <- switch (class(obj),
                  "glo.fit" = apply(t(glo.rl.gradient(a=a, p=f)), 1, ismev::q.form, m = mat),
                  "gev.fit" = apply(t(ismev::gev.rl.gradient(a=a, p=f)), 1, ismev::q.form, m = mat))
  # v <- v[seq(length(v),1)]
  if(plot.out){
    mydots <- list(...)  
    print(!any(names(mydots) == "xlim"))
    print(!any(names(mydots) == "ylim"))
    if(!any(names(mydots) == "ylim")) ylim <- range(c(y, dat))
    if(!any(names(mydots) == "xlim")) xlim <- range(c(xv, xemp))
    ### no limits given
    if((!any(names(mydots) == "xlim")) & (!any(names(mydots) == "ylim"))) plot(xv, y, type = "n", ylim = ylim, xlim = xlim, ...)
    ## only ylim given
    if(!any(names(mydots) == "ylim") & any(names(mydots) == "xlim"))  plot(xv, y, type = "n", ylim = ylim, ...)
    ### only xlim given 
    if(!any(names(mydots) == "xlim") & any(names(mydots) == "ylim"))  plot(xv, y, type = "n", xlim = xlim, ...)
    ### all limits given 
    if(any(names(mydots) == "xlim") & any(names(mydots) == "ylim"))  plot(xv, y, type = "n", ...)
    lines(xv, y, ...)
    lines(xv, y + qnorm(1-sign.alpha/2) * sqrt(v), ...)
    lines(xv, y - qnorm(1-sign.alpha/2) * sqrt(v), ...)
    if(!(nh > 0)){
      points(xemp, sort(dat), pch = pch, ...)
  }  
  if(nh > 0){
    ### following Bayliss and Reed
    nex <- length(dat[dat>X0])   # Number of observations > X0 (both in systematic and hist data when present)
    n<-length(dat)-nk+nh 
    fplot1<-NULL 
    # fplot1 is the nex points that exceed the perception threshold 
    # (both in systematic and hist data)
    if(nex!=0) fplot1<-((1:nex-0.44)/(nex + 1 - 0.88))*(nex/n)
    # flpot2 is the points that do not exceed the perception threshold in systematic data
    fplot2 <- 
      (nex/n)+((n-nex)/n)*((seq((nex+1),length(dat))-nex-0.44)/(n-nh-length(dat[-seq(0,nk)][dat[-seq(0,nk)]>X0])+1-0.88))
    points(xcalc(1-fplot2, vtype = vtype),sort(dat[dat<X0], decreasing=TRUE), pch = pch, ...)
    ### we use a different symbol for the historical data 
    ### and it gets messy
    dat[1:nk]<-(dat[1:nk]+0.00001) ### to make sure we identify the right ones!!!
    histpoints<-seq(1,length(dat[dat>X0]))[(sort(dat[dat>X0]) %in% dat[1:nk])]
    fplot3 <- sort(fplot1[length(dat[dat>X0])+1-histpoints])
    points(xcalc(1-fplot3, vtype = vtype),
           sort((sort(dat[dat>X0]))[histpoints],decreasing = TRUE), pch=pchHist, ... )
    fplot4 <- sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])
    points(xcalc(1-fplot4, vtype = vtype),
           sort((sort(dat[dat>X0]))[-histpoints],decreasing = TRUE),pch = pch, ...)
  }
  }
  invisible(data.frame(p = 1-f, retLev = y, se = sqrt(v))[f %in% fline, ])
}


# xl <- ifelse(any(names(list(...)) == "xlab"), 
#              list(...)["xlab"],  
#              switch(vtype,
#                     "Gumbel" = "Gumbel variate: -log(-log(p))",
#                     "redVar" = "Reduced variate: log(p/(1-p))"))


#' @title Fitted return curve for a glo.fit object
#' @description This function mimics the \code{gev.rl} function for GLO models.
#' It plots the flood frequency curve (return curve) based on the data and the fitted parameters for a \code{glo.fit} model, including ones with historical data. 
#' Very ad-hoc and working under assumption that flood data are plotted (hence the y-axis lab).
#' Also outputs useful informations. 
#' Mostly written by Thomas Kjeldsen. 
#' This function is deprecated and has been replaced by the generic retPlot function.
#' @param a the mle estimates from a \code{glo.fit} object
#' @param mat the covariance function from a \code{glo.fit} object
#' @param dat data matrix from a \code{glo.fit} object
#' @param nh if historical data are used, the length of historical record - default is 0, no historical data
#' @param nk if historical data are used, the number of peaks above X0 - default is 0, no historical data
#' @param X0 if historical data are used, the perception threshold - default is NULL, no historical data
#' @param f frequencies for which return level (and 95\%) confidence intervals are calculated
#' @export
#' @seealso \code{\link{retPlot}}
#' @examples 
#' set.seed(7821567)
#' xx <- rglo(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' h1 <- glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, show=FALSE)
#' s1 <- glod.fit(xxsist, show=FALSE) 
#' rls1 <- glo.rl(a=s1$mle,mat=s1$cov,dat=s1$data)
#' rlh1 <- glo.rl(h1$mle,h1$cov,h1$data,nh=h1$h,nk=h1$k,X0=h1$X0)
#' lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl+1.96*sqrt(rlh1$var), lty = 2)
#' lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl-1.96*sqrt(rlh1$var), lty = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl, col = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl-1.96*sqrt(rls1$var), lty = 2, col = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl+1.96*sqrt(rls1$var), lty = 2, col = 2)
#' legend("topleft", col =c(1,2),legend = c("With historical","Systematic Only"), lty = 1)
#' ## similar fitted curve - but obvious reduction in uncertainty 
glo.rl <- function(a, mat,dat,nh=0,nk=0,X0=NULL, f=c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 0.95, 0.99, 0.995, 0.999)) {
    #
    # produces return level curve and 95 % confidence intervals 
    # on usual scale
    #
    # TRKJ, 18 December
    # a: maximum likelihood estimates of 3 parameters
    # mat: covariance matrix
    # dat: amax data time series
    # nh:  period covered by historical data
    # nk: number of historic data, 0 = no historic data
    # X0: NULL no perception threshold !NULL no perception threshold (actual values)
    #
    # f non-exceedance probablity P(X<x), thus log(f/(1-f)) = log(T-1), where T is return period
    #
    #
    if(is.null(X0)  & nh>0) stop("historic data with no perception threshold")
    if(!is.null(X0) & nh<1) stop("perception threshold with no historic data")
    eps <- 1e-006
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    ## quantile function
    gloqq<-
      function (a, p) 
      {
        if (a[3] != 0) 
          a[1] + (a[2] * (1-((1 - p)/p)^(-a[3])))/a[3]
        else a[1] - a[2] * (log((1 - p)/p))
      }
    # Points on the quantile function
    q <- gloqq(a, 1 - f)
    d <- t(glo.rl.gradient( a=a, p=1-f))
    v <- apply(d, 1, q.form, m = mat); print(v)
    
    # set-up FEH-style plotting surface diagram
    par(mgp=c(2.75,1,0))
    plot(log(f/(1-f)), q,type = "n", 
         xlab = paste(paste("Logistic reduced variate","\n"),"Ln(T-1)"), 
         ylab = expression(paste("Peak flow ",(m^3/s))),
         bty="l",xaxt = "n")
    axis(1)
    # Plot return period axis
    # Axis is plotted at the level of the 1.11-year event (10th entry in the q-array)
    # Upper limit currently fixed at 200 years. Maybe should be a plotting option?
    #
    yaxis<-q[10]
    
    Taxis<-c(2,5,10,25,50,100,200)
    et_axis<-log(Taxis-1)
    paxis<-rep(yaxis,length(Taxis))
    
    #points(gLaxis,paxis,
    points(et_axis,paxis,
           pch="|",cex=0.75)
    lines(x=c(log(2-1),log(200-1)),y=c(yaxis,yaxis),
          lwd=1)
    
    rtp<-c("2","5","10","25","50","100","200")
    #text(gLaxis,paxis-5,rtp,cex=0.5)
    text(et_axis,paxis-0.20*yaxis,rtp,cex=0.75)
    text(4.0,1.5*yaxis,"Return period [years]",cex=0.75)
    
    
    # Plot FFC curves
    #
    lines(log(f/(1-f)), q)
    #   x<-log(f/(1-f))
    #   curve((a1+(a2/a3)*(1-exp(x)^(-a3))),from=1,to=499,n=500,add=TRUE)
    #     points(-1/(log((1:length(dat))/(length(dat) + 1))), sort(dat))
    #     lines(1/(1-f), q,col=2)
    #     lines(-1/log(f), q)
    #     lines(1/(1-f), q + 1.96 * sqrt(v), col = 4)
    #     lines(1/(1-f), q - 1.96 * sqrt(v), col = 4)
    #### change from ismev: use Gringorten plotting position
    # points(-1/(log((1:length(dat)-0.44)/(length(dat) + 1 - 0.88))), sort(dat),col=1)
    nex<-ifelse(is.null(X0),0,length(dat[dat>X0]))   # Number of observations > X0
    n<-length(dat)-nk+nh 
    fplot1<-NULL 
    # fplot1 is the nex points that exceed the perception threshold 
    # (both in systematic and hist data)
    if(nex!=0) fplot1<-((1:nex-0.44)/(nex + 1 - 0.88))*(nex/n)
    # flpot2 is the points that do not exceed the perception threshold in systematic data
    fplot2<- 
      (nex/n)+((n-nex)/n)*((seq((nex+1),length(dat))-nex-0.44)/(n-nh-length(dat[-seq(0,nk)][dat[-seq(0,nk)]>X0])+1-0.88))
    #    points(-1/log(1-c(fplot1,fplot2)),sort(dat,decreasing = TRUE),col=2,pch=19)
    if(is.null(X0)) X0<-max(dat)+1
    ##    points((-1/log(1-fplot2)),sort(dat[dat<X0],decreasing=TRUE),col=1,pch=1,cex=0.9)
    points(log((1-fplot2)/(fplot2)),sort(dat[dat<X0],decreasing=TRUE),col=1,pch=1,cex=0.9)
    ### we use a different symbol for the historical data 
    ### and it gets messy
    dat[1:nk]<-(dat[1:nk]+0.00001) ### to make sure we identify the right ones!!!
    histpoints<-seq(1,length(dat[dat>X0]))[(sort(dat[dat>X0]) %in% dat[1:nk])]
    if(nh>0){  ## if historical data exists
      #### this is not very pretty, but it is correct,it's all because of sorting issues
      ##      points(-1/log(1-sort(fplot1[length(dat[dat>X0])+1-histpoints])),
      fplot3<-(1-sort(fplot1[length(dat[dat>X0])+1-histpoints]))
      points(log(fplot3/(1-fplot3)),
             sort((sort(dat[dat>X0]))[histpoints],decreasing = TRUE),
             col=1,pch=17,cex=0.9)
      ##      points(-1/log(1-sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])),
      fplot4<-1-sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])
      cat(log(fplot4/(1-fplot4)))
      points(log(fplot4/(1-fplot4)),
             sort((sort(dat[dat>X0]))[-histpoints],decreasing = TRUE),
             col=1,pch=1,cex=0.9)
    }
    #    points(-1/(log((1:length(dat)-0.44)/(length(dat) + 1 - 0.88))), sort(dat),col=1)   
    ##    title("Return Level Plot")        
    ## add an output, so we can use the rl again
    res<-list(f=f,rl=q,var=v,dat=dat,n=length(dat),mle=a)
    res
}




#' @title GLO quantile function
#' @description This function gives quantiles the GLO distribution. It's built mirroring the \code{gev.fit} function of \code{ismev}. This is the same as qglo, just used to mimic ismev
#' @param loc,scale,sh the locatio, scale, shape parameters
#' @param p the non-exceedance probabiliy
gloq<-function (p, loc, scale, sh) {
  if (sh != 0) loc + (scale * (1-((1 - p)/p)^(-sh)))/sh
  else loc - scale * (log((1 - p)/p))
}


glo.rl.gradient<-function (a, p) {
  scale <- a[2]
  shape <- a[3]
  out <- matrix(NA, nrow = 3, ncol = length(p))
  out[1, ] <- 1
  
  yp <- p/(1 - p)
  out[2, ] <- shape^(-1) * (1 - yp^(shape))
  # out[3, ] <- scale * (shape^(-2)) * (1 - yp^(-shape)) -
  #   scale * shape^(-1) * yp^(-shape) * log(yp)
  out[3, ] <- - scale * (shape^(-2)) * (1 - yp^(shape)) -
    scale * shape^(-1) * yp^(shape) * log(yp)
  return(out)
}




