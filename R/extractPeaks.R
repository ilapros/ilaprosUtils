if(! existsFunction("turnpoints") )
# this corresponds to pastecs::turnpoints 
# copied here since no other function from the package is used
turnpoints<-function (x) {
  data <- deparse(substitute(x))
  if (is.null(ncol(x)) == FALSE) 
    stop("Only one series can be treated at a time")
  x <- as.vector(x)
  n <- length(x)
  diffs <- c(x[1] - 1, x[1:(n - 1)]) != x
  uniques <- x[diffs]
  n2 <- length(uniques)
  poss <- (1:n)[diffs]
  exaequos <- c(poss[2:n2], n + 1) - poss - 1
  if (n2 < 3) {
    warning("Less than 3 unique values, no calculation!")
    nturns <- NA
    firstispeak <- FALSE
    peaks <- rep(FALSE, n2)
    pits <- rep(FALSE, n2)
    tppos <- NA
    proba <- NA
    info <- NA
  }
  else {
    m <- n2 - 2
    ex <- matrix(uniques[1:m + rep(3:1, rep(m, 3)) - 1], 
                 m)
    peaks <- c(FALSE, apply(ex, 1, max, na.rm = TRUE) == 
                 ex[, 2], FALSE)
    pits <- c(FALSE, apply(ex, 1, min, na.rm = TRUE) == ex[, 
                                                           2], FALSE)
    tpts <- peaks | pits
    if (sum(tpts) == 0) {
      nturns <- 0
      firstispeak <- FALSE
      peaks <- rep(FALSE, n2)
      pits <- rep(FALSE, n2)
      tppos <- NA
      proba <- NA
      info <- NA
    }
    else {
      tppos <- (poss + exaequos)[tpts]
      tptspos <- (1:n2)[tpts]
      firstispeak <- tptspos[1] == (1:n2)[peaks][1]
      nturns <- length(tptspos)
      if (nturns < 2) {
        inter <- n2 + 1
        posinter1 <- tptspos[1]
      }
      else {
        inter <- c(tptspos[2:nturns], n2) - c(1, tptspos[1:(nturns - 
                                                              1)]) + 1
        posinter1 <- tptspos - c(1, tptspos[1:(nturns - 
                                                 1)])
      }
      posinter2 <- inter - posinter1
      posinter <- pmax(posinter1, posinter2)
      proba <- 2/(inter * gamma(posinter) * gamma(inter - 
                                                    posinter + 1))
      info <- -log(proba, base = 2)
    }
  }
  res <- list(data = data, n = n, points = uniques, pos = (poss + 
                                                             exaequos), exaequos = exaequos, nturns = nturns, firstispeak = firstispeak, 
              peaks = peaks, pits = pits, tppos = tppos, proba = proba, 
              info = info)
  class(res) <- "turnpoints"
  res
}





#' @title Extracting peaks
#' @description Extracts independent peaks
#' @param vecTime an optional numeric vector containing information on the time variable. 
#' The vector is used to order the vecObs variable, so it should be provided if the data are not in chronological order
#' If \code{vecTime = NULL} (the default) it is assumed that the data are correctly ordered. 
#' @param vecObs the vector of observations, the variable from which peaks should be extracted 
#' @param mintimeDiff the time after which peaks should be considered independent
#' @param thrConst the constant defining how much the through should be smaller than the peaks to identify independent peaks 
#' @name extractPeaks
#' @details This functions has been (sort of) tested in extracting peaks from 15 minutes data. 
#' The function assumes that two peaks are independent if they are at least \code{mintimeDiff} steps apart and
#' if the minimumn amount of water in the though is less than \code{thrConst} the value of the maximum peak.
#' \code{mintimeDiff} can not be given as an actual time measuraments, it is just an indicator of how many postions apart 
#' in the \code{vecTIme} vector two peaks need to be to be considered independent. 
#' Therefore \code{vecTIme} needs to be equally spaced, 
#' e.g. if missing values are present in the orignal file they should be infilled. 
#' @return a list of 0-1 values which indicate whether the \code{vecObs}value is a peak 
#' @examples 
#' zz <- data.frame(hour=seq(1, 240),
#'                  flow=sample(c(rnorm(200, 50, 5),rgev(40, 70,3,0.2))))
#' zz$isPeak <- extractPeaks(vecTime=zz$hour,vecObs=zz$flow,mintimeDiff=4)
#' plot(zz$hour, zz$flow, type = "l")
#' points(zz$hour[zz$isPeak==1], zz$flow[zz$isPeak==1], col=2)
#' ## more hours needed to be assumed independent
#' zz$isPeak <- extractPeaks(vecTime=zz$hour,vecObs=zz$flow,mintimeDiff=8)
#' points(zz$hour[zz$isPeak==1], zz$flow[zz$isPeak==1], col=4)
#' @export
extractPeaks<-function(vecTime = NULL, vecObs, mintimeDiff=73, thrConst=(2/3)){ 
  #### selects peaks at least mintimeDiff apart
  ##  relies on turnpoints from pastecs
  if(is.null(vecTime)) vecTime <- seq_along(vecObs)
  if(!is.numeric(vecTime)) warning("Time should be numeric")
  tt<-as.data.frame(cbind(vecTime,vecObs))
  names(tt)<-c("time","obs")
  tt<-tt[order(tt$time),]
  cc<-turnpoints(tt$obs)  ### from pastecs
  sub<-as.data.frame(cbind(tt$time[cc$pos],tt$obs[cc$pos],cc$peaks,cc$pits))
  names(sub)<-c("time","flow","peak","pit")
  keep<-rep(0,nrow(sub))  ### indicator on whether the obs is a peak
  ObsEA<-seq(1,nrow(tt))[-cc$pos] ### position of exaeqos - needed later
  obsPeaks<-(seq(1,nrow(sub))[sub$peak==1])
  obsPits<-(seq(1,nrow(sub))[sub$pit==1])
  keep[obsPeaks[1]]<-1
  for(i in 2:length(obsPeaks)){
    now<-obsPeaks[i]
    prev<-max(obsPeaks[1:(i-1)][obsPeaks[1:(i-1)] %in% seq(1,length(keep))[keep==1]])
    keep[now]<-1
    tp<-sub$time[prev];fp<-sub$flow[prev] ### time and flow at previous peak
    tn<-sub$time[now]; fn<-sub$flow[now]  ### time and flow at current peak
    ## issue if the two peaks happen to be the same values, results could be changing at different runs
    ## as max.col will identify different columns as being the max
    if(fn==fp) smallest<-prev 
    else smallest <- c(prev,now)[!c(prev,now) %in% seq(prev,now)[max.col(t(sub[prev:now,2]))]]
    if((tn-tp) < mintimeDiff) keep[smallest]<-0
    else{
      minThrough<-min(sub$flow[obsPits[obsPits<now & obsPits>prev]])
      if(minThrough > max(fp,fn)*(thrConst)) keep[smallest]<-0
    }
  }  
  isPeak<-rep(0,nrow(tt))
  isPeak[cc$pos]<-keep
  isPeak[ObsEA]<- 0
  isPeak
}

