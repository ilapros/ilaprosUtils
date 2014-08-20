#' @title Functions to transform dates into other useful information
#' @name dates transformation 
#' @param datevec a vetor of dates (\code{class(x)} must be Date)
#' @return The year, water year or month from the input date. Or a 0-1 indicator of summer.
#' @export
date2year <- function(datevec) as.numeric(substring(datevec,first=1,last=4)) 
#' @name dates transformation
#' @family manipulate dates
#' @export
date2wy <- function(datevec, endMonth = 9) as.numeric(substring(datevec,first=1,last=4)) - 1*as.numeric(as.numeric(substring(datevec,first=6,last=7) ) <= endMonth)
#' @name dates transformation
#' @family manipulate dates
#' @export
date2month <- function(datevec)  as.numeric(substring(datevec,first=6,last=7))
#' @name dates transformation
#' @family manipulate dates
#' @export
date2summer <- function(datevec, begMonth = 4, endMonth = 9) as.numeric(date2month(datevec) > (begMonth-1) & date2month(datevec) <= endMonth)



