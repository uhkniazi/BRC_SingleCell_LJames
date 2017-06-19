# File: utilities.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global functions used in some scripts
# Date: 12/06/2017


### calculate model fits
## first write the log predictive density function
lpd = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## likelihood function with posterior theta
  return(sum(dbinom(resp, 1, iFitted, log=T)))
}


## calculate one data point at a time
## log pointwise predictive density
lppd = function(theta, data){
  betas = t(theta) # matrix of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = matrix(data$mModMatrix, nrow = 1, byrow = T)
  # calculate fitted value
  iFitted = as.vector(mModMatrix %*% betas)
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## likelihood function with posterior theta
  return(mean(dbinom(resp, 1, iFitted, log=F)))
}

getStanSD = function(obj){
  return(apply(extract(obj)$betas, 2, sd))
}
getStanMean = function(obj){
  return(apply(extract(obj)$betas, 2, mean))
}
getStanPValue = function(obj){
  pnorm(-abs(getStanMean(obj)/getStanSD(obj)))*2
}

## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  return(iFitted)
}