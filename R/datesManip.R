#' Add together two numbers
#'
#' @param x a vector
#' @return The length of \code{x} without the NA
#' @examples
#' x <- c(5,12,NA,4)
#' lengthNOTna(x)
lengthNOTna <- function(x) length(x[!is.na(x)])


#' Functions to transform dates into other useful information
#' @name datas transformation 
#' @family manipulate dates
#' @param datevec a vetor of dates (\code{class(x)} must be Date)
#' @return The length of \code{x} without the NA
#' @export
date2year <- function(datevec) as.numeric(substring(datevec,first=1,last=4)) 
#' @name datas transformation
#' @family manipulate dates
#' @export
date2wy <- function(datevec, endMonth = 9) as.numeric(substring(datevec,first=1,last=4)) - 1*as.numeric(as.numeric(substring(datevec,first=6,last=7) ) <= endMonth)
#' @name datas transformation
#' @family manipulate dates
#' @export
date2month <- function(datevec)  as.numeric(substring(datevec,first=6,last=7))
#' @name datas transformation
#' @family manipulate dates
###@aliases date2year date2month date2summer date2wy
#' @export
date2summer <- function(datevec, begMonth = 4, endMonth = 9) as.numeric(date2month(datevec) > (begMonth-1) & date2month(datevec) <= endMonth)

