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

