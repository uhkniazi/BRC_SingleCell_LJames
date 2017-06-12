# File: 16_signatureGenesCellTypes.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: identify signature genes to classify cells
# Date: 05/06/2017


## set variables and source libraries
source('header.R')

## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did2
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 13) AND (MetaFile.comment like "%NanoString%")')
dfSample = dbGetQuery(db, q)
dfSample

# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location[2], dfSample$name[2])

load(n)

## make count matrix
names(lNanoString)
lNanoString$desc

mCounts.norm = lNanoString$normalizedData$normalized.data
mCounts.norm = mCounts.norm[mCounts.norm$Code.Class == 'Endogenous', ]
rownames(mCounts.norm) = mCounts.norm$Name
mCounts = as.matrix(mCounts.norm[,-c(1:3)])

dim(mCounts)
i = which(rowSums(mCounts) == 0)
mCounts = mCounts[-i,]
dim(mCounts)
mCounts = scale(t(mCounts))
head(apply(mCounts, 2, sd))

############### follow the regression based approach
## remove correlated variables first
## find correlated variables
mCor = cor(mCounts, use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.6, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})
i = which(colnames(mCounts) %in% n)
cvTopGenes = colnames(mCounts)[-i]

mCounts = mCounts[,cvTopGenes]
rm(mCor)
gc()

########## perform binomial regression with shrinkage 
###### on each category 
dfData = data.frame(mCounts)
colnames(dfData) = gsub('\\.', '_', colnames(dfData))


## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'CD19'] = 1
fGroups[dfData$fCellID != 'CD19'] = 0

dfData$fCellID = factor(fGroups)


## setup functions to fit model and calculate log posterior
## and log likelihood
logit.inv = function(p) {exp(p)/(exp(p)+1) }

## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

library(LearnBayes)
source('utilities.R')
library(parallel)
tryCombinations = function(iCombinationIndex){
  ## setup the data
  f = paste('fCellID ~ ', paste(mCombinations[,iCombinationIndex], collapse='+'), collapse = ' ')
  lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1),
               mModMatrix=model.matrix(as.formula(f), data=dfData))
  
  # set starting values for optimiser
  start = c(rep(0, times=ncol(lData$mModMatrix)))
  names(start) = colnames(lData$mModMatrix)
  # fit model
  fit.lap = laplace(mylogpost, start, lData)
  ### lets take a sample from this 
  ## parameters for the multivariate t density
  tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=5)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 5000, lData)
  colnames(s) = colnames(lData$mModMatrix)
  #fit.lap$sir = s
  
  ## averages of posterior from sir sample
  post = apply(s, 2, mean)
  
  # calculate AIC
  iAIC = (lpd(post, lData) - 4) * -2
  
  # calculate E(lpd(theta))
  eLPD = mean(sapply(1:nrow(s), function(x) lpd(s[x,], lData)))
  
  # calclate ilppd
  ilppd = sum(log(sapply(seq_along(lData$resp), function(x) {
    d = list(resp=lData$resp[x], mModMatrix = lData$mModMatrix[x,])
    lppd(s, d)
  })))
  
  ## effective numbers of parameters pWAIC1
  pWAIC1 = 2 * (ilppd - eLPD)
  
  iWAIC = -2 * (ilppd - pWAIC1)
  
  fit.lap$modelCheck = list('AIC'=iAIC, 'pWAIC1'=pWAIC1, 'WAIC'=iWAIC)
  return(fit.lap)
}

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(cvTopGenes, 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(cvTopGenes, 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(cvTopGenes, 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(cvTopGenes, 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

lCD19 = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

names(lFits.1var) = 1:length(lFits.1var)
mOne = t(do.call(cbind, lapply(lFits.1var, function(x) unlist(x$modelCheck))))

names(lFits.2var) = 1:length(lFits.2var)
mTwo = t(do.call(cbind, lapply(lFits.2var, function(x) unlist(x$modelCheck))))

names(lFits.3var) = 1:length(lFits.3var)
mThree = t(do.call(cbind, lapply(lFits.3var, function(x) unlist(x$modelCheck))))

names(lFits.4var) = 1:length(lFits.4var)
mFour = t(do.call(cbind, lapply(lFits.4var, function(x) unlist(x$modelCheck))))

iAIC = c(min(mOne[,'AIC']), min(mTwo[,'AIC']), min(mThree[,'AIC']), min(mFour[,'AIC']))
pWAIC = c(min(mOne[,'pWAIC1']), min(mTwo[,'pWAIC1']), min(mThree[,'pWAIC1']), min(mFour[,'pWAIC1']))
WAIC = c(min(mOne[,'WAIC']), min(mTwo[,'WAIC']), min(mThree[,'WAIC']), min(mFour[,'WAIC']))




lFits.3var = lapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

table(sapply(lFits.3var, is.null))

lFits.3var$mCombinations = mCombinations
names(lFits.3var) = 1:(length(lFits.3var)-1)
lFits.3var.sub = lFits.3var[!sapply(lFits.3var, is.null)]
length(lFits.3var.sub)

iWAIC3v = sapply(1:(length(lFits.3var.sub)-1), function(x){
  lFits.3var.sub[[x]]$modelCheck$WAIC
})

iWAIC3v.param = sapply(1:(length(lFits.3var.sub)-1), function(x){
  lFits.3var.sub[[x]]$modelCheck$pWAIC
})

iAIC3v = sapply(1:(length(lFits.3var.sub)-1), function(x){
  lFits.3var.sub[[x]]$modelCheck$AIC
})

mFitParam.3v = cbind(iWAIC3v, iWAIC3v.param, iAIC3v)
rownames(mFitParam.3v) = names(lFits.3var.sub)[1:(length(lFits.3var.sub)-1)]

########### try 2 variables
mCombinations = combn(cvTopGenes, 2)

lFits.2var = lapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

table(sapply(lFits.2var, is.null))

lFits.2var$mCombinations = mCombinations
names(lFits.2var) = 1:(length(lFits.2var)-1)
lFits.2var.sub = lFits.2var[!sapply(lFits.2var, is.null)]
length(lFits.2var.sub)

iWAIC2v = sapply(1:(length(lFits.2var.sub)-1), function(x){
  lFits.2var.sub[[x]]$modelCheck$WAIC
})

iWAIC2v.param = sapply(1:(length(lFits.2var.sub)-1), function(x){
  lFits.2var.sub[[x]]$modelCheck$pWAIC
})

iAIC2v = sapply(1:(length(lFits.2var.sub)-1), function(x){
  lFits.2var.sub[[x]]$modelCheck$AIC
})

mFitParam.2v = cbind(iWAIC2v, iWAIC2v.param, iAIC2v)
rownames(mFitParam.2v) = names(lFits.2var.sub)[1:(length(lFits.2var.sub)-1)]


########### try 4 variables
mCombinations = combn(cvTopGenes, 4)

library(parallel)

lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

table(sapply(lFits.4var, is.null))

lFits.4var$mCombinations = mCombinations
names(lFits.4var) = 1:(length(lFits.4var)-1)
lFits.4var.sub = lFits.4var[!sapply(lFits.4var, is.null)]
length(lFits.4var.sub)

iWAIC4v = sapply(1:(length(lFits.4var.sub)-1), function(x){
  lFits.4var.sub[[x]]$modelCheck$WAIC
})

iWAIC4v.param = sapply(1:(length(lFits.4var.sub)-1), function(x){
  lFits.4var.sub[[x]]$modelCheck$pWAIC
})

iAIC4v = sapply(1:(length(lFits.4var.sub)-1), function(x){
  lFits.4var.sub[[x]]$modelCheck$AIC
})

mFitParam.4v = cbind(iWAIC4v, iWAIC4v.param, iAIC4v)
rownames(mFitParam.4v) = names(lFits.4var.sub)[1:(length(lFits.4var.sub)-1)]









f = paste('fCellID ~ ', paste(mCombinations[,7703], collapse='+'), collapse = ' ')
  
fit.bin = glm(as.formula(f), data=dfData, family = binomial(link='logit'))
summary(fit.bin)



# library(car)
# ## utility function to calculate hyper-prior variance parameters
# ## gamma shape function
# ## this function is from the book: Bayesian data analysis
# gammaShRaFromModeSD = function( mode , sd ) {
#   # function changed a little to return jeffery non-informative prior
#   if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
#   rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
#   shape = 1 + mode * rate
#   return( list( shape=shape , rate=rate ) )
# }
# unlist(gammaShRaFromModeSD(mode = logit(sd(lData$resp+0.5)/2), 
#                            sd = logit(2*sd(lData$resp+0.5))))


## try with stan first
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=3, 
                    pars=c('betas'), cores=3)
# some diagnostics for stan
print(fit.stan, digits=3)
plot(fit.stan)
m = extract(fit.stan)
b = m$betas
s = m$betaSigma



traceplot(fit.stan, ncol=1, nrow=6, inc_warmup=F)
#stan_diag(fit.stan)
## some sample diagnostic plots
library(coda)
oCoda = As.mcmc.list(fit.stan)
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])







## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood 
  lp = dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

## set starting values and data
f = paste('fCellID ~ ', paste(colnames(dfData)[1:2], collapse='+'), collapse = ' ')
lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1), mModMatrix=model.matrix(fCellID ~  ., data=dfData))

start = c(betas=rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.lap = laplace(mylogpost, start, lData)
op = optimx(start, mylogpost, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)
summary(op) ##rvmmin seems to converge better

library(rstan)
stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'betaSigma'))
# some diagnostics for stan
print(fit.stan, digits=3)
plot(fit.stan)
traceplot(fit.stan, ncol=1, nrow=6, inc_warmup=F)
#stan_diag(fit.stan)
## some sample diagnostic plots
library(coda)
oCoda = As.mcmc.list(fit.stan)
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])























################################ old


dfData = data.frame(t(mCounts))
## add the cell type id
dfData$fCellID = factor(lNanoString$metaData$group2)

library(randomForest)
fit.rf = randomForest(fCellID ~ ., data=dfData)

# get variables importance
varImpPlot(fit.rf)
dfRF = data.frame(importance(fit.rf))
head(dfRF)
ivScore = dfRF$MeanDecreaseGini
names(ivScore) = rownames(dfRF)
ivScore = sort(ivScore, decreasing = T)
ivScore
summary(ivScore)
f = quantile(ivScore, 0.90)
table(ivScore >= f)
ivScore = ivScore[ivScore >= f]
tail(ivScore)
## keep the top names
## for some reason the - in the name is replaced by a .
## replace those first
cvTopGenes = gsub('\\.', '-', names(ivScore))
table(cvTopGenes %in% rownames(mCounts))

## find correlated variables
mCounts = mCounts[cvTopGenes,]
mCor = cor(t(mCounts), use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.6, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})
i = which(cvTopGenes %in% n)
cvTopGenes.cor = cvTopGenes[-i]

mCounts = mCounts[cvTopGenes.cor,]

dfData = data.frame(scale(t(mCounts)))
## add the cell type id
dfData$fCellID = factor(lNanoString$metaData$group2)

library(MASS)

fit.lda = lda(fCellID ~ ., data=dfData)
plot(fit.lda)
#############################################################################################

## try with binomial regression
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'CD19'] = 1
fGroups[dfData$fCellID != 'CD19'] = 0
dfData$fCellID = factor(fGroups )
dfData$patient = factor(lNanoString$metaData$group1)

library(lme4)
# define a modification of the laplace function from learnbayes
library(LearnBayes)
library(car)
library(numDeriv)
library(optimx)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

f = paste(c('fCellID ~ ', paste(colnames(dfData)[1:2], collapse='+'), '+ (1|patient)'), collapse = ' ')
fit.glm = glmer(f, data=dfData, family = binomial)

summary(fit.glm)


myloglike.random = function(theta, data){
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
  ## parameters to track/estimate
  sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[grep('ran', names(theta))]# random effects jitters for the group deflections
  
  ## random effect jitter for the population intercept
  # each group contributes a jitter centered on 0
  # population slope + random jitter
  ivBetaRand = betas[1] + iGroupsJitter
  # create a matrix of betas with the new interceptr/unique intercept for each random effect
  ivIntercept = ivBetaRand[groupIndex] # expand this intercept
  iFitted = mModMatrix[,2:ncol(mModMatrix)] %*% betas[2:ncol(mModMatrix)]
  iFitted = ivIntercept + iFitted
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## first level, log priors
  lhp = dunif(sigmaRan, 0, 2, log=T)
  lran = sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T))
  lp = sum(dcauchy(betas, 0, 10, log=T))
  # write the likelihood function
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lhp + lran + lp + lik
  return(val)
}



mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=5000), method='Nelder-Mead', data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}

## set starting values and data
f = paste('fCellID ~ ', paste(colnames(dfData)[1:2], collapse='+'), collapse = ' ')
lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1), mModMatrix=model.matrix(fCellID ~  MME+CXCR3, data=dfData))
lData$groupIndex = as.numeric(dfData$patient)

start = c(sigmaRan=log(1), betas=rep(0, times=ncol(lData$mModMatrix)), 
          ran=rep(0, times=length(unique(lData$groupIndex))))

myloglike.random(start, lData)

fit.lap = mylaplace(myloglike.random, start, lData)

op = optimx(start, myloglike.random, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)
summary(op) ##rvmmin seems to converge better

mylaplace2 = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optimx(mode, logpost, control = list(maximize=T, usenumDeriv=T), method='Rvmmin', data=data)
  # calculate hessian
  hes = hessian(logpost, coef(fit), data=data)
  colnames(hes) = names(mode)
  rownames(hes) = names(mode)
  options(warn = 0)
  mode = coef(fit)
  h = -solve(hes)
  stuff = list(mode = mode, var = h, converge = fit$convcode == 
                 0)
  return(stuff)
}

####################################################### try without random effects 
## it doesnt seem to be required

f = paste('fCellID ~ ', paste(colnames(dfData)[1:2], collapse='+'), collapse = ' ')
fit.bin = glm(f, data=dfData, family = binomial(link='logit'))

## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood 
  lp = dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

## set starting values and data
f = paste('fCellID ~ ', paste(colnames(dfData)[1:2], collapse='+'), collapse = ' ')
lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1), mModMatrix=model.matrix(fCellID ~  ., data=dfData))

start = c(betas=rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.lap = laplace(mylogpost, start, lData)
op = optimx(start, mylogpost, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)
summary(op) ##rvmmin seems to converge better

library(rstan)
stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'betaSigma'))
# some diagnostics for stan
print(fit.stan, digits=3)
plot(fit.stan)
traceplot(fit.stan, ncol=1, nrow=6, inc_warmup=F)
#stan_diag(fit.stan)
## some sample diagnostic plots
library(coda)
oCoda = As.mcmc.list(fit.stan)
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])

