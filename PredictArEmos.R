# Function to predict with fitted AR-EMOS on a fixed static test set
# Modification of the implementation in the R package ensAR 
# that is adapted to the specific needs of the training (and test) data used for 
# the benchmark comparison

# Arguments: 
# fit - An object from the ar_emos_fit function containing the fitted model 
#       from a static training set
# ens - A data frame (possibly new test data) with one observation column, 
#       at least one forecast column, and at least one additional column (e.g. date)
# obs_col - The column number containing the observations
# mem_col - The column number(s) containing the ensemble forecast(s)	


ar_emos_predict <- function(fit, ens, obs_col, mem_col)
{

require(ensAR)

# Observation and ensemble forecasts, in the benchmark test data it is 50 members
y <- ens[, obs_col]
ens.dat <- ens[, mem_col]
n <- length(y)

# Error series of mean of ensemble and error series of ctrl forecast
z.mean <- (y - apply(ens.dat[,-1], 1, mean))
z.ctrl <- (y - ens$ctrl)


# Compute AR-corrected ensemble members and empirical variance of
# AR-corrected ensemble on (test) data, based on the parameter estimates from the fit
# object obtained on static training data

# The maximum of the 2 estimated orders 
p.max <- max(fit$order.ctrl, fit$order.mean)


# AR-corrected control member
ens.pred.l <- length(ens$ctrl[-(1:p.max)])  
xtilde.ctrl <- NULL

for(i in 1:ens.pred.l)
{
xtilde.ctrl[i] <- ens$ctrl[i+p.max] + fit$intercept.ctrl + 
sum(fit$coefs.ctrl * (z.ctrl[i:(i+fit$order.ctrl-1)] - fit$intercept.ctrl), na.rm=TRUE)
}

# Variance of the AR process of the error series of ctrl member
varAR.ctrl <- var_ar(ar = fit$coefs.ctrl, i_var = fit$var.ar.ctrl)


# AR-corrected 50 ensemble members 
# Apply the AR-correction based on parameters estimated from the ensemble mean
# to each of the members individually 
ens.dat.ex <- ens.dat[,-1]
n.ens.pred <- length(mem_col[-1])
xtilde.ens <- matrix(NA, ncol=n.ens.pred, nrow=ens.pred.l)
for(j in 1:n.ens.pred)
{
for(i in 1:ens.pred.l)
{
xtilde.ens[i,j] <- ens.dat.ex[i+p.max, j] + fit$intercept.mean + 
sum(fit$coefs.mean * (z.mean[i:(i+fit$order.mean-1)] - fit$intercept.mean), na.rm=TRUE)
}
}

# Variance of the AR process of the error series of mean of the 50 members
varAR.ens <- var_ar(ar = fit$coefs.mean, i_var = fit$var.ar.mean)


# Resulting AR-corrected ensemble members on (test) data and AR variance
# for the full ensemble containing control member and the 50 other members
xtilde <- data.frame(xtilde.ctrl, xtilde.ens)
varAR <- cbind(varAR.ctrl, varAR.ens)

# Resulting mean of AR-corrected ensemble and empirical sd of 
# AR-corrected ensemble on (test) data
sample_sd <- function(x) sqrt(mean(scale(x, scale = FALSE)^2))
mu_pred <- apply(xtilde, 1, mean)
sd_pred_1 <- sqrt(apply(varAR, 1, mean))
sd_pred_2 <- apply(xtilde, 1, sample_sd) 
sd_pred <- fit$w *sd_pred_1 + (1-fit$w) * sd_pred_2


# Output: predictive mean and sd 
return(list(pred.mean=mu_pred, pred.sd=sd_pred))

 }