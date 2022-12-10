# Function to fit AR-EMOS for a fixed static training period
# Modification of the implementation in the R package ensAR 
# that is adapted to the specific needs of the training (and test) data used for 
# the benchmark comparison

# Arguments: 
# ens - A data frame with one observation column, at least one forecast column, 
#       and at least one additional column (e.g. date)
# obs_col - The column number containing the observations
# mem_col - The column number(s) containing the ensemble forecast(s)	
# skip - A numerical value corresponding to the forecast horizon 
# (0 for horizons not greater than 24 hours, 1 for horizons greater than 
# 24 hours and not greater than 48 hours, and so on).

ar_emos_fit <- function (ens, obs_col, mem_col, skip = 0) 
{
require(ensAR)

    if (!is.data.frame(ens)) 
        stop("ensemble input must be a data frame")
        
# In case individual/small numbers of successive missing values occur 
# in the training data time series these are imputed by spline interpolation
# A longer series of missing values cannot be reasonably handeled
   
    if (any(!complete.cases(ens[, c(obs_col, mem_col)]))) {
        warning("occurring NA's were replaced by spline interpolation")
        ens_approx <- apply(ens[, c(obs_col, mem_col)], 2, zoo::na.spline, 
            na.rm = FALSE)
        ens[, obs_col] <- ens_approx[, 1]
        ens[, mem_col] <- ens_approx[, -1]
    }

# Observations
    y <- ens[, obs_col]


# AR parameters are estimated oce on static training data,
# one set of pars for mean of 10 ensemble members 
# and one set of pars for ctrl member.
# This is a modification compared to the orginal implementation in ensAR, 
# where the AR process is estimated for each member indiviudally.
# This modification was necessary as the number of members in training and test
# data differs, however, both data sets contain a control member

# Error series of mean of the 10 members   
z.mean <- (y - apply(ens[,mem_col[-1]], 1, mean))
# Error series of control member
z.ctrl <- (y - ens$ctrl)

             # For higher forecast horizons the observations und thus the errors
             # at one or more past time points (depending on the horizon) 
             # are not yet observed, and are thus predicted by an AR process,
             # see paper on AR-Emos for details
            if (skip > 0) 
            {
                zs.mean <- z.mean[1:(length(z.mean) - skip)]
                zs.mean_mod <- ar(x = zs.mean, aic = TRUE, order.max = NULL)
                zs.mean_pred <- predict(object = zs.mean_mod, n.ahead = skip, 
                  se.fit = FALSE)
                z.mean <- c(zs.mean, zs.mean_pred)
                
                zs.ctrl <- z.ctrl[1:(length(z.ctrl) - skip)]
                zs.ctrl_mod <- ar(x = zs.ctrl, aic = TRUE, order.max = NULL)
                zs.ctrl_pred <- predict(object = zs.ctrl_mod, n.ahead = skip, 
                  se.fit = FALSE)
                z.ctrl <- c(zs.ctrl, zs.ctrl_pred)
            }
            
            # Estimate AR process of error series of mean of 10 members, 
            # here order p of the process is restricted to 15.
            # In the original ensAR implementation the order is not restricted,
            # this is a modifcation made to have more homogeneous results
            # for the benchmark comparison
            z.mean_mod <- ar(x = z.mean, aic = TRUE, order.max = 15)
            p.mean <- z.mean_mod$order
            mu.mean <- z.mean_mod$x.mean
            a.mean <- z.mean_mod$ar
            ar_var.mean <- z.mean_mod$var.pred
            
            # Estimate AR process of error series of ctrl member, 
            # here order p of the process is restricted to 15 (see above)
            z.ctrl_mod <- ar(x = z.ctrl, aic = TRUE, order.max = 15)
            p.ctrl <- z.ctrl_mod$order
            mu.ctrl <- z.ctrl_mod$x.mean
            a.ctrl <- z.ctrl_mod$ar
            ar_var.ctrl <- z.ctrl_mod$var.pred



# Compute AR-corrected ensemble members on training data and empirical variance 
# of AR-corrected ensemble on training data

# The maximum of the 2 estimated orders 
p.max <- max(p.ctrl, p.mean)

# AR-corrected control member
ens.new.l <- length(ens$ctrl[-(1:p.max)])  
xtilde.ctrl <- NULL

for(i in 1:ens.new.l)
{
xtilde.ctrl[i] <- ens$ctrl[i+p.max] + mu.ctrl + sum(a.ctrl * (z.ctrl[i:(i+p.ctrl-1)] - mu.ctrl))
}

# Variance of the AR process of the error series of ctrl member
varAR.ctrl <- var_ar(ar = a.ctrl, i_var = z.ctrl_mod$var.pred)

# AR-corrected ensemble members 
# Apply the AR-correction based on parameters estimated from the ensemble mean
# to each of the 10 members individually 
ens.dat <- ens[,mem_col[-1]]
n.ens.train <- length(mem_col[-1])
xtilde.ens <- matrix(NA, ncol=n.ens.train, nrow=ens.new.l)
for(j in 1:n.ens.train)
{
for(i in 1:ens.new.l)
{
xtilde.ens[i,j] <- ens.dat[i+p.max, j] + mu.mean + sum(a.mean * (z.mean[i:(i+p.mean-1)] - mu.mean))
}
}

# Variance of the AR process of the error series of mean of the 10 members
varAR.ens <- var_ar(ar = a.mean, i_var = z.mean_mod$var.pred)

# Resulting AR-corrected ensemble members on training data and AR variance
# for the full ensemble containing control member and the 10 other members
xtilde <- data.frame(xtilde.ctrl, xtilde.ens)
varAR <- cbind(varAR.ctrl, varAR.ens)

# Resulting mean of AR-corrected ensemble and empirical sd of 
# AR-corrected ensemble on training data
sample_sd <- function(x) sqrt(mean(scale(x, scale = FALSE)^2))
mu_out <- apply(xtilde, 1, mean)
sd_out_1 <- sqrt(apply(varAR, 1, mean))
sd_out_2 <- apply(xtilde, 1, sample_sd) 


# crps objective function to estimate weights for combining the 2 variance terms
   m_crps <- function(par) 
    {
                m_obs <- y[-(1:p.max)]
                m_mu <- mu_out
                m_sd <- par[1] * sd_out_1 + (1 - par[1]) * sd_out_2
                z <- (m_obs - m_mu)/m_sd
                vals <- m_sd * (z * (2 * pnorm(z) - 1) + 2 * 
                  dnorm(z) - 1/sqrt(pi))
                mean(vals)
    }

# Numerical optimization    
res <- optim(par = c(w = 0.5), fn = m_crps, method = "L-BFGS-B", 
                lower = 0, upper = 1)

res.w <- as.vector(res$par["w"])


# Output: estimated orders and coefficients of the AR-processes, and estimated weight 
#         when combining the 2 variance terms
return(list(order.mean=p.mean, intercept.mean=mu.mean, coefs.mean=a.mean, var.ar.mean=ar_var.mean,
       order.ctrl=p.ctrl, intercept.ctrl=mu.ctrl, coefs.ctrl=a.ctrl, var.ar.ctrl=ar_var.ctrl, w=res.w))


}