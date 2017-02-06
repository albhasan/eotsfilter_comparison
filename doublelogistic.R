# Estimation of NDVI using a double logistic function. Taken from Global land surface phenology trends from GIMMS database by Julien & Sobrino 2009. Equation 1 
#
# @param doy    A numeric. Day of the year (0:364)
# @param w.ndvi A numeric. Winter NDVI value
# @param m.ndvi A numeric. Maximum NDVI value
# @param S      A numeric. Increasing inflection point (spring date)
# @param A      A numeric. Drecreasing inflection point (autumn date)
# @param mS     A numeric. Rate of increase at S
# @param mA     A numeric. Rate of decrease at A
# @return       A numeric. An estimation od NDVI for a given day of the year
dlog1 <- function(doy, w.ndvi, m.ndvi, S, A, mS, mA){
  # Global land surface phenology trends from GIMMS database by Julien & Sobrino. Equation 1 
  # TEST: plot(dlog1(doy = 0:364, w.ndvi = 0.07, m.ndvi = 0.68, S = 119, A = 282, mS = 0.19, mA = 0.13), type = "l")
  return(w.ndvi + (m.ndvi - w.ndvi) * ((1/(1 + exp(-mS * (doy - S)))) + (1/(1 + exp( mA * (doy - A)))) - 1))
}



# Estimation of NDVI using a double logistic function for arid and semi-arid areas. Taken from Global land surface phenology trends from GIMMS database by Julien & Sobrino 2009. Equation 2
#
# @param doy    A numeric. Day of the year (0:364)
# @param w.ndvi A numeric. Winter NDVI value
# @param m.ndvi A numeric. Maximum NDVI value
# @param S      A numeric. Increasing inflection point (spring date)
# @param A      A numeric. Drecreasing inflection point (autumn date)
# @param mS     A numeric. Rate of increase at S
# @param mA     A numeric. Rate of decrease at A
# @return       A numeric. An estimation od NDVI for a given day of the year
dlog2 <- function(doy, w.ndvi, m.ndvi, S, A, mS, mA){
  # Global land surface phenology trends from GIMMS database by Julien & Sobrino. Equation 1 
  return(m.ndvi - (m.ndvi - w.ndvi) * ((1/(1 + exp(-mS * (doy - S)))) + (1/(1 + exp( mA * (doy - A)))) - 1))
}






# Transforms a date into the year-day_of_the_year date (YYYYDOY)
#
# @param date.vec   A Date vector
# @return           A character representing a date as day-of-the-year (YYYYDOY)
date2ydoy <- function(date.vec){
  yearOriginText <- paste(format(date.vec, "%Y"), "/01/01", sep="")
  yearOrigin <- as.POSIXlt(yearOriginText)
  doy <- as.numeric(as.Date(date.vec) - as.Date(yearOrigin)) + 1
  res <- paste(format(date.vec, "%Y"), sprintf("%03d", doy), sep="")
  return(res)
}
