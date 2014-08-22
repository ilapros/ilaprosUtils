
#' @title Extracting maximas
#' @description Given a data.frame or matrix, keeps rows which correspond to the maximum for each level of indFac. 
#' @param zz the dataset of matrix whcih should be processed
#' @param indFac the column number which identifies the grouping column. 
#' The column doesn't need to be a factor (\code{as.factor} enforced). 
#' Use ordere factors if the data need to be in the same order. 
#' @param indVec the column number which identifies the variable of interest column.
#' @name inst2daily
#' @details As the name betrays this was built to process 15-minutes data and identify daily maxima. 
#' It still work to extract monthy or annual maxima. 
#' 
#' 
#' 
#' Could probably be achived via dplyr as well.
#' @return an object of the same class and size as zz, 
#' in which only values wich are the maximum value for the indVec column in the indFac groups are kept.
#' @examples 
#' set.seed(145)
#' zz <- data.frame(month = c(rep("Jan",5),rep("Feb",5),rep("Mar",5),rep("Apr",5)), 
#'                 flow = rnorm(20,50,6),
#'                 day=rep(c(2,9,16,23,28),4))
#' inst2daily(zz, indVec=2, indFac = 1)
#' zz$month <- factor(zz$month, levels = c("Jan","Feb","Mar","Apr"), ordered = TRUE)
#' inst2daily(zz, indVec=2, indFac = 1)
#' @export
inst2daily<-function(zz, indFac, indVec){
  ## extract from the instantenuos record a daily record
  ## only implemented for max
  if(!is.ordered(zz[,indFac])) zz[,indFac] <- as.factor(as.character(zz[,indFac]))
  zz <- zz[order(zz[,indFac]), ]
  ind <- as.numeric(tapply(cbind(zz[,indVec]),
                         INDEX=zz[,indFac], FUN=which.max))
  ll <- as.numeric(tapply(cbind(zz[,indVec]),
                        INDEX=zz[,indFac], FUN=length))
  ll <-c(0, cumsum(ll[-length(ll)]))
  ind <- ind + ll
  out <- zz[ind,]
  return(out)
}





