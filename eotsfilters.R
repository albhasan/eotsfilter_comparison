################################################################################
# preliminars
################################################################################
# install required packages
# install.packages(c("devtools", "knitr", "ptw", "tidyr", "dplyr", "optimx"))
# library(devtools)
# install_github("luizassis/wtss.R")
library(wtss.R)
library(knitr)
library(ptw)
library(zoo)
library(tidyr)
library(dplyr)
library(ggplot2)
library(optimx)

# colors for the color blinded
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
################################################################################
# get data from WTSS
################################################################################
veg_index <- "evi"
latitude <- -14.919100049
longitude <- -59.11781088

chronos <- WTSS("http://www.dpi.inpe.br/tws/wtss")
coverages <- listCoverages(chronos)
cv <- describeCoverage(chronos,coverages[2])
ts <- timeSeries(chronos,
                 names(cv),
                 attributes = c(veg_index, "reliability", "quality"),
                 latitude = latitude,
                 longitude = longitude,
                 start="2000-02-18",
                 end="2016-01-01"
)
data.zoo <- ts[[1]][[2]]
data.df <- as.data.frame(coredata(data.zoo), stringsAsFactors = FALSE)
data.df["date"] <- index(data.zoo)
data.df[, veg_index] <- data.df[, veg_index]/10000
################################################################################
# whitaker
################################################################################
lambda.vec <- c(1e0, 1e1, 1e2, 1e3)                                             # whitaker's lambdas
plot.df <- data.frame()                                                         # data to plot
data.wt1.df <- data.wt2.df <- data.df[, c("date", "reliability", "quality", "evi")]  # re-order columns
#-----------------------------------
# plot first order whitaker
#-----------------------------------
whit1.df <- lapply(lambda.vec, function(lambda, veg_index, data.df){            # compute whitaker
  return(whit1(y = as.vector(unlist(data.df[veg_index])), lambda = lambda))
},
veg_index = veg_index,
data.df = data.df
)
whit1.df <- do.call("cbind", whit1.df)
colnames(whit1.df) <- paste("whit1_", lambda.vec, sep = "")
data.wt1.df <- cbind(data.wt1.df, whit1.df)                                     # formated data frame
plot.df <- gather(data.wt1.df, smoothing, value, evi:whit1_1000)                # plot data.frame
ggplot(plot.df, aes(date, value, group = smoothing, color = smoothing)) +
  geom_line() +
  labs(
    title = paste(toupper(veg_index), " using Weighted Whittaker smoothing", sep = " "),
    subtitle="First order finite difference penalty") +
  scale_colour_manual(values=cbPalette)
#-----------------------------------
# plot second order whitaker
#-----------------------------------
whit2.df <- lapply(lambda.vec, function(lambda, veg_index, data.df){            # compute whitaker
  return(whit2(y = as.vector(unlist(data.df[veg_index])), lambda = lambda))
},
veg_index = veg_index,
data.df = data.df
)
whit2.df <- do.call("cbind", whit2.df)
colnames(whit2.df) <- paste("whit2_", lambda.vec, sep = "")
data.wt2.df <- cbind(data.wt2.df, whit2.df)                                     # formated data frame
plot.df <- gather(data.wt2.df, smoothing, value, evi:whit2_1000)                # plot data.frame
ggplot(plot.df, aes(date, value, group = smoothing, color = smoothing)) +
  geom_line() +
  labs(
    title = paste(toupper(veg_index), " using Weighted Whittaker smoothing", sep = " "),
    subtitle = "Second order finite difference penalty") #+
  scale_colour_manual(values=cbPalette)
################################################################################
# Apply the Fourier transformation as a filter to the time-series
################################################################################
source("/home/alber/Documents/Dropbox/alberLocal/inpe/projects/tsFiltering/fourier.R")
#-----------------------------------
# FFT parameters
#-----------------------------------
trajectory <- as.vector(unlist(data.df[veg_index]))
time <- 1
acq.freq <- (length(trajectory) - 1)/time
ts <- seq(0, time, 1/acq.freq)
#-----------------------------------
# de-trend the time-series
#-----------------------------------
trend <- lm(trajectory ~ ts)
detrended.trajectory <- trend$residuals
obs.pred <- predict(trend)
#-----------------------------------
# plot raw & detrended data
#-----------------------------------
#plot.df <- data.frame(obs = c(trajectory, detrended.trajectory, obs.pred), x = data.df$date)
#plot.df["type"] <- rep(c("original", "de-trended", "trend"), each = length(trajectory))
#ggplot(plot.df, aes(x = x, y = obs, color = type)) + geom_line() + labs(title = "Original versus de-trended series")
#-----------------------------------
# plot harmonics
#-----------------------------------
X.k <- fft(detrended.trajectory)
#plot.frequency.spectrum(X.k, xlimits = c(0, acq.freq / 2))
#-----------------------------------
# build the filtered signal
#-----------------------------------
keep <- c(0.10, 0.20, 0.50, 0.75) # ratio of frequencies to keep
plot.df <- data.frame(
  x = data.df$date,
  obs = trajectory,
  type = rep(veg_index, length(trajectory))
)
#-----------------------------------
# build the plot
#-----------------------------------
for(r in keep){
  X.k.tmp <- X.k
  X.k.tmp[round(r * length(X.k)):length(X.k)] <- 0                              # unwanted higher frequencies become 0
  x.n <- get.trajectory(X.k.tmp, ts, acq.freq) / acq.freq                       # TODO: why the scaling?
  x.n <- x.n + (trajectory - trend$residuals)                                   # add back the trend
  plot.df <- rbind(plot.df, data.frame(x = data.df$date, obs = Re(x.n), type = paste("keep", r, sep = " ")))
}
ggplot(plot.df, aes(x = x, y = obs, color = type)) +
  geom_line() +
  labs(
    title = paste(toupper(veg_index), " using Fourier filtering", sep = " "),
    subtitle="Keeping a ratio of the lower frequencies"
  ) +
  scale_colour_manual(values=cbPalette)
################################################################################
# double logistic
################################################################################
source("/home/alber/Documents/Dropbox/alberLocal/inpe/projects/tsFiltering/doublelogistic.R")
#-----------------------------------
# add the Day Of the Year as a column
#-----------------------------------
data.df["doy"] <- as.integer(
  substr(date2ydoy(date.vec = data.df$date), start = 5, stop = 7)
)
#-----------------------------------
# move the first day of the year 6 months when in the south hemisphere
#-----------------------------------
if(latitude < 0){
  halfyear <- as.numeric(
    substr(date2ydoy(date.vec = c(as.Date("2010-07-01"))), start = 5, stop = 7)
  )
  data.df["doy"] <- data.df[, "doy"] + halfyear
}
#-----------------------------------
# get the double-logistic function parameters
#-----------------------------------
dlog.df <- data.frame(x = data.df$doy, y = data.df[, veg_index])                # data to fit
start1 <- c(w.ndvi = 0.07, m.ndvi = 0.68, S = 119, A = 282, mS = 0.19, mA = 0.13) # initial parameters for optimization algorithm
dlog1.f <- function(b, mydata){                                                 # function to optimize
  sum((mydata$y - dlog1(mydata$x, w.ndvi = b[1], m.ndvi = b[2], S = b[3], A = b[4], mS = b[5], mA = b[6]))^2)
}
dlog.optx <- optimx(                                                            # optimization
  par = start1, fn = dlog1.f, mydata = dlog.df,
  control = list(
    method = "BFGS",
    #all.methods =  TRUE,
    save.failures = TRUE,
    maxit = 2500
  )
)
# dlog.optx # review results of optimization - all.methods =  TRUE
params <- as.vector(unlist(dlog.optx["BFGS", 1:6]))
data.df["dlog"] <- dlog1(doy = data.df$doy, w.ndvi = params[1], m.ndvi = params[2], S = params[3], A = params[4], mS = params[5], mA = params[6])
#-----------------------------------
# plot
#-----------------------------------
plot.df <- rbind(
  data.frame(date = data.df$date, val = data.df[, veg_index], type = veg_index),
  data.frame(date = data.df$date, val = data.df[, "dlog"], type = "D-Log")
)
ggplot(plot.df, aes(x = date, y = val, group = type, color = type)) +
  geom_line() +
  labs(
    title = paste("Double logistic function", sep = " "),
    subtitle="Adjusted using Day-Of-Year 0:365"
  ) +
  scale_colour_manual(values=cbPalette)
################################################################################
# KFAS - structural model
################################################################################
library(KFAS)
fit <- StructTS(data.df[, veg_index], type = "level")
plot.df <- rbind(
  data.frame(date = data.df$date, val = data.df[, veg_index], type = veg_index),
  data.frame(date = data.df$date, val = as.numeric(fit$fitted), type = "SSM")
)
ggplot(plot.df, aes(x = date, y = val, group = type, color = type)) +
  geom_line() +
  labs(
    title = paste("Fitted structural model", sep = " "),
    subtitle="Using maximum likehood"
  ) +
  scale_colour_manual(values=cbPalette)
################################################################################
# 
################################################################################





