# Analysis with Benchmark training and test data 
# The fit and predict function are applied to each forcast horizon and 
# each station separately (local estimation)
 

# Install R package ensAR from github, as some functions of the package
# are used in this code
library(devtools)
install_github("JuGross/ensAR")
library(ensAR)


# Load preprocessed R data sets
load(file="benchmark_t2m_train.Rdata")
load(file="benchmark_t2m_test.Rdata")

# Source fit and predict functions
source(file="FitArEmos.R")
source(file="PredictArEmos.R")


# Extract all existing leadtimes and stations in both data sets
# (here extracted from training data, they are the same in the test data)
leadtimes <- unique(benchmark_t2m_train$leadtime)
stations <- unique(benchmark_t2m_train$id) 

# Number of different leadtimes
n.leadtimes <- length(leadtimes)

# For leadtimes up to 72h the argument skip in ar_emos_fit is set to 0,
# for leatimes larger 72h it is set to 1
skip.vals <- c(rep(0,13), rep(1,8))


# Number of ensemble members in test data
m <- 51
# Quantile levels at which a sample from the predictive distribution is obtained
qlevels <- (1:m)/(m+1)

# Output Array with 51 quantiles (ensemble forecasts) from the predictive distribution 
# on all test dates (730), for each station (229) and each leadtime (21)
ArEmos_quantiles <- array(0, dim = c(m, n.leadtimes, 730, 229))


# Loop through all leadtimes in the data
for(h in 1:n.leadtimes)
{

print(h)

# Obtain subset of data corresponding to the leadtime of current loop step
dat <- subset(benchmark_t2m_train, subset=(leadtime==leadtimes[h]))
dat.test <- subset(benchmark_t2m_test, subset=(leadtime==leadtimes[h]))

# Delete data from 2017 in training set, as the test set also contains 2017
# data, this way, the data is not used twice for training and testing
ind.2017 <- which(substr(dat$init_date, start=1, stop=4) == "2017")
dat <- dat[-ind.2017,]

# Total number of dates in training set after removing 2017
init.dates <- unique(dat$init_date)
n.dates.train <- length(init.dates)


# Count proportion of NaN/NA for each station
freq.na <- numeric(length(stations))
for(i in 1:length(stations))
{
print(i)
dat.i <- subset(dat, subset=(id==stations[i])) 
n <- dim(dat.i)[1]
freq.na[i] <- sum(is.na(dat.i$observation))/n
}

# Stations with more than 50% NA will not be used for fitting at all
# set predictions in test data for those to NA later
ind.station <- which(freq.na <= 0.5)
ind.station.out <- which(freq.na > 0.5)


# Aditionally, the station numbers 89, 100, 113 are not used for fitting as they
# contain further long blocks of NA which deteriorate the parameter estimates 
# Results in the following index vectors for stations where the model is 
# fitted and stations where the model is not fitted
ind.station <- ind.station[-which(ind.station %in% c(89,100,113))]
ind.station.out <- sort(c(ind.station.out, 89, 100, 113))

# Loop through all stations for which <= 50% of data non-NA
for(s in ind.station)
{

print(s)

# Obtain training data for the station of current loop step
dat.s <- subset(dat, subset=(id==stations[s])) 
# Due to large blocks of sucessive NA in the beginning of the data for
# most stations, which cannot be reasonably handeled by AR-EMOS
# the first 2200 observations are discarded for the estimation to make sure that  
# no long blocks of NA are present at any station. The number of data points
# to delete in order to have no problems at any station was determined 
# in a preliminary analysis
dat.s <- dat.s[-(1:2200),]
n.s <- dim(dat.s)[1]
# For station number 152 (id "bla") further NA are present
# therefore here some more data is removed
if(s=="152"){dat.s <- dat.s[-((n.s-625):n.s),]}

# Obtain test data for the station of current loop step
# test data only contains few single NA, so no additional data needs to be
# discarded
dat.s.test <- subset(dat.test, subset=(id==stations[s])) 
n.s.test <- dim(dat.s.test)[1]

# Fit the model on training data for current station and leadtime
train.s <- ar_emos_fit(ens=dat.s, obs_col=14, mem_col=15:25 , skip = skip.vals[h]) 

# Predict on test data
predict.s <- ar_emos_predict(fit=train.s, ens=dat.s.test, obs_col=14, mem_col=15:65) 

# Obtain equidistant quantiles from predictive normal distribution with
# estimated mean and sd from AR-EMOS approach
quants <- function(x){qnorm(x, mean=predict.s$pred.mean, sd=predict.s$pred.sd)}
quant.pred.s <- sapply(qlevels, quants)

# When making predictions for a time point t on test data, the data at 
# past time points t-1,..., t-p, where p is the order of the estimated AR process, 
# are needed to compute the predition at time point t.
# Thus, the very first p observations in the data have to be deleted (set to NA/NaN), 
# as a prediction with AR-EMOS is only possible from time point t+1 on.  
n.NA <- n.s.test - length(predict.s$pred.mean)
NA.mat <- matrix(NaN, nrow=n.NA, ncol=m)
quant.pred.s <- rbind(NA.mat, quant.pred.s)

# Assign the predictive quantiles for current leadtime h and station s
# to appropriate position in the output array
ArEmos_quantiles[ ,h, ,s] <- t(quant.pred.s)
# In order to make the end result more homogeneous 
# set the first 15 (maximal order of AR process set in the fit function)
# time points to NA/NaN although the actual order might be lower
ArEmos_quantiles[ ,h, 1:15 ,s] <- NaN

# end for s
}

# Those stations for which >50% of the data is NA, are not fitted and predicted.
# Aditionally, set the results for station numbers 89, 100, 113, as they
# contain further long blocks of NA which deteriorate the parameter estimates  
for(l in ind.station.out)
{
ArEmos_quantiles[ ,h, ,l] <- NaN 
}

# end for h
}
