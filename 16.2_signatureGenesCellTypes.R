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

####### load single cell data
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12) AND (MetaFile.comment like "%scater%")')
dfSample = dbGetQuery(db, q)
dfSample
# close connection after getting data
dbDisconnect(db)
library(scater)
library(org.Hs.eg.db)
n = paste0(dfSample$location, dfSample$name)
load(n)

# load the single cell normalised count matrix
mCounts.SC = exprs(oSce.F)

# load the nano string count matrix
names(lNanoString)
lNanoString$desc

mCounts.norm = lNanoString$normalizedData$normalized.data
mCounts.norm = mCounts.norm[mCounts.norm$Code.Class == 'Endogenous', ]
rownames(mCounts.norm) = mCounts.norm$Name

## match the names between the 2 data sets 
## i.e. symbols in nanostring and enterez id in single cell
dfSymbols = AnnotationDbi::select(org.Hs.eg.db, rownames(mCounts.norm), 
                                  columns = 'ENTREZID', keytype = 'SYMBOL')
dfSymbols = na.omit(dfSymbols)
# how many match b/w the 2 data sets
table(dfSymbols$ENTREZID %in% rownames(mCounts.SC))
# 
# FALSE  TRUE 
# 293   262 
dfSymbols = dfSymbols[(dfSymbols$ENTREZID %in% rownames(mCounts.SC)),]
# match these symbols with nanostring data
i = match(dfSymbols$SYMBOL, rownames(mCounts.norm))
identical(dfSymbols$SYMBOL, rownames(mCounts.norm)[i])
mCounts.norm = mCounts.norm[i,]
# replace these names by enterez ids
rownames(mCounts.norm) = dfSymbols$ENTREZID

## create the count matrix from nano string data for variable selection
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
# data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
i = which(colnames(mCounts) %in% n)
cvTopGenes = colnames(mCounts)[-i]

mCounts = mCounts[,cvTopGenes]
rm(mCor)
gc()

## cvTopGenes = 16 genes

########## perform binomial regression on each category 
dfData = data.frame(mCounts)
# this conversion to data.frame tends to put an X before variable
# names so use that when making formulas

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
  tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=4)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 5000, lData)
  colnames(s) = colnames(lData$mModMatrix)
  #fit.lap$sir = s
  
  ## averages of posterior from sir sample
  post = apply(s, 2, mean)
  
  # calculate AIC
  iAIC = (lpd(post, lData) - length(start)) * -2
  
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
# variable names need an X
mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lCD19 = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

########## repeat the variable selection for other groups, GC
dfData = data.frame(mCounts)

## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'GC'] = 1
fGroups[dfData$fCellID != 'GC'] = 0
table(fGroups, dfData$fCellID)

dfData$fCellID = factor(fGroups)

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lGC = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

########## repeat the variable selection for other groups, Mem
dfData = data.frame(mCounts)


## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'Mem'] = 1
fGroups[dfData$fCellID != 'Mem'] = 0
table(fGroups, dfData$fCellID)

dfData$fCellID = factor(fGroups)

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lMem = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

########## repeat the variable selection for other groups, Naive
dfData = data.frame(mCounts)


## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'Naive'] = 1
fGroups[dfData$fCellID != 'Naive'] = 0
table(fGroups, dfData$fCellID)

dfData$fCellID = factor(fGroups)

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lNaive = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

########## repeat the variable selection for other groups, PB
dfData = data.frame(mCounts)


## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'PB'] = 1
fGroups[dfData$fCellID != 'PB'] = 0
table(fGroups, dfData$fCellID)

dfData$fCellID = factor(fGroups)

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lPB = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

########## repeat the variable selection for other groups, PreGC
dfData = data.frame(mCounts)


## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'PreGC'] = 1
fGroups[dfData$fCellID != 'PreGC'] = 0
table(fGroups, dfData$fCellID)

dfData$fCellID = factor(fGroups)

## generate the combination matrix
## using 3-4 variables at the most i.e. log(18)

mCombinations = combn(paste('X', cvTopGenes, sep=''), 1)
lFits.1var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


mCombinations = combn(paste('X', cvTopGenes, sep=''), 2)
lFits.2var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 3)
lFits.3var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})

mCombinations = combn(paste('X', cvTopGenes, sep=''), 4)
lFits.4var = mclapply(1:ncol(mCombinations), function(iIndexSub) {
  tryCatch(tryCombinations(iIndexSub), error=function(e) NULL)
})


names(lFits.1var) = 1:length(lFits.1var)
names(lFits.2var) = 1:length(lFits.2var)
names(lFits.3var) = 1:length(lFits.3var)
names(lFits.4var) = 1:length(lFits.4var)
# save the object
lPreGC = list(one=lFits.1var, two=lFits.2var, three=lFits.3var, four=lFits.4var)

rm(list = c('lFits.1var', 'lFits.2var', 'lFits.3var', 'lFits.4var'))
gc()

## save the objects
lVarSelection = list(lCD19, lGC, lMem, lNaive, lPB, lPreGC)

n = make.names(paste('Binomial Variable Selection for NanoString Data from louisa single cell project rds'))
lVarSelection$desc = paste('Binomial Variable Selection for NanoString Data from louisa single cell project', date())
names(lVarSelection) = c('cd19', 'gc', 'mem', 'naive', 'pb', 'pregc', 'desc')
n2 = paste0('~/Data/MetaData/', n)

save(lVarSelection, file=n2)

# comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did2, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='Binomial Variable Selection for NanoString Data from louisa single cell project')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

## check each result
lFits.1var = lVarSelection$pb$one
lFits.2var = lVarSelection$pb$two
lFits.3var = lVarSelection$pb$three
lFits.4var = lVarSelection$pb$four

mOne = t(do.call(cbind, lapply(lFits.1var, function(x) unlist(x$modelCheck))))

mTwo = t(do.call(cbind, lapply(lFits.2var, function(x) unlist(x$modelCheck))))

mThree = t(do.call(cbind, lapply(lFits.3var, function(x) unlist(x$modelCheck))))

mFour = t(do.call(cbind, lapply(lFits.4var, function(x) unlist(x$modelCheck))))

# iAIC = c(min(mOne[,'AIC']), min(mTwo[,'AIC']), min(mThree[,'AIC']), min(mFour[,'AIC']))
# pWAIC = c(min(mOne[,'pWAIC1']), min(mTwo[,'pWAIC1']), min(mThree[,'pWAIC1']), min(mFour[,'pWAIC1']))
# WAIC = c(min(mOne[,'WAIC']), min(mTwo[,'WAIC']), min(mThree[,'WAIC']), min(mFour[,'WAIC']))

### make some plots
boxplot(mOne[,'AIC'], mTwo[,'AIC'], mThree[,'AIC'], mFour[,'AIC'])
boxplot(mOne[,'pWAIC1'], mTwo[,'pWAIC1'], mThree[,'pWAIC1'], mFour[,'pWAIC1'])
boxplot(mOne[,'WAIC'], mTwo[,'WAIC'], mThree[,'WAIC'], mFour[,'WAIC'])

boxplot(mOne[,'AIC'], mOne[,'WAIC'], mTwo[,'AIC'], mTwo[,'WAIC'], mThree[,'AIC'], mThree[,'WAIC'],
        mFour[,'AIC'], mFour[,'WAIC'], col=rep(grey.colors(2), times=4), main='Scores for model vs model size',
        xlab='No. of variables', ylab='Score', xaxt='n', pch=20, cex=0.5)
axis(1, at = 1:8, labels = c(1, 1, 2, 2, 3, 3, 4, 4))
legend('bottomleft', legend = c('AIC', 'WAIC'), fill=grey.colors(2))

### select the variables with lowest scores in each model size
getAICVar = function(m, l){
  iA = which.min(m[,'AIC'])
  names(l[[names(iA)]]$mode)[-1]  
}

getWAICVar = function(m, l){
  iW = which.min(m[,'WAIC'])
  names(l[[names(iW)]]$mode)[-1]
}

## get the top variables based on WAIC in each comparison
lVarSelection$desc = NULL

lTopVariables = lapply(lVarSelection, function(lx){
  ## get the matrix 
  lmats = lapply(lx, function(lx2){
    t(do.call(cbind, lapply(lx2, function(x) unlist(x$modelCheck))))
  })
  lwaic = lapply(seq_along(lmats), function(lx2){
    return(getWAICVar(lmats[[lx2]], lx[[lx2]]))
  })
})

## predict on these variables and test in single cell data
dfData = data.frame(mCounts)
# this conversion to data.frame tends to put an X before variable
# names so use that when making formulas

## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'PreGC'] = 1
fGroups[dfData$fCellID != 'PreGC'] = 0

dfData$fCellID = factor(fGroups)

## calculate coefficients for the top variables and model size
lCoef = lapply(lTopVariables$pregc, function(lx){
  f = paste('fCellID ~ ', paste(lx, collapse='+'), collapse = ' ')
  lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1),
               mModMatrix=model.matrix(as.formula(f), data=dfData))
  # set starting values for optimiser
  start = c(rep(0, times=ncol(lData$mModMatrix)))
  names(start) = colnames(lData$mModMatrix)
  # fit model
  fit.lap = laplace(mylogpost, start, lData)
  ### lets take a sample from this 
  ## parameters for the multivariate t density
  tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=4)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 5000, lData)
  colnames(s) = colnames(lData$mModMatrix)
  ## averages of posterior from sir sample
  post = apply(s, 2, mean)
  return(post)
})

## new data, the same training data
dfData.new = data.frame(mCounts)
## perform prediction on this
lPred = lapply(lCoef, function(lx){
  X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,names(lx)[-1]]))
  colnames(X) = names(lx)
  lData = list(mModMatrix=X)
  return(mypred(lx, lData))
})

dfPred.train = do.call(cbind, lPred)
dfPred.train = round(dfPred.train, 3)
dfPred.train = data.frame(dfPred.train)
dfPred.train$actual = dfData$fCellID

###### looking at the pWAIC, WAIC and prediction errors on the training data
## a 3 variable model appears to be generally good enough 
## use this 3 variable model to make predictions in single cell data
dfData.new = data.frame(scale(t(mCounts.SC)))
dfData.new = dfData.new[,names(lCoef[[3]])[-1]]
dim(dfData.new)
## perform prediction on this
lPred = lapply(lCoef[3], function(lx){
  X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,names(lx)[-1]]))
  colnames(X) = names(lx)
  lData = list(mModMatrix=X)
  return(mypred(lx, lData))
})

#dfSingleCellPred = data.frame(cd19=lPred[[1]])
dfSingleCellPred$pregc = lPred[[1]]

n = make.names(paste('Predictions for single cell classes using binomial classification louisa rds'))
n2 = paste0('~/Data/MetaData/', n)

save(dfSingleCellPred, file=n2)

## comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='Predictions for single cell classes using binomial classification louisa j project')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

############################################# 
##### second classifier using random forest and lda
dfData = data.frame(mCounts)
# this conversion to data.frame tends to put an X before variable
# names so use that when making formulas

## add the cell type id
dfData$fCellID = factor(lNanoString$metaData$group2)
table(dfData$fCellID)
dim(dfData)

set.seed(123)
library(randomForest)
fit.rf = randomForest(fCellID ~ ., data=dfData)

# get variables importance
varImpPlot(fit.rf)
dfRF = data.frame(importance(fit.rf))
head(dfRF)
ivScore = dfRF$MeanDecreaseGini
names(ivScore) = rownames(dfRF)
ivScore = sort(ivScore, decreasing = T)
head(ivScore)
summary(ivScore)
f = quantile(ivScore, 0.75)
table(ivScore >= f)
ivScore = ivScore[ivScore >= f]
tail(ivScore)
## keep the top names
cvTopGenes = gsub('X', '', names(ivScore))
table(cvTopGenes %in% colnames(mCounts))
## find correlated variables
mCor = cor(mCounts[,cvTopGenes], use="na.or.complete")
# library(caret)
# ### find the columns that are correlated and should be removed
# n = findCorrelation((mCor), cutoff = 0.6, names=T)
# data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# i = which(cvTopGenes %in% n)
# cvTopGenes.cor = cvTopGenes[-i]
# 
# mCounts = mCounts[cvTopGenes.cor,]

## train the model
dfData = data.frame(mCounts[,cvTopGenes])
# this conversion to data.frame tends to put an X before variable
# names so use that when making formulas
## add the cell type id
dfData$fCellID = factor(lNanoString$metaData$group2)
table(dfData$fCellID)
dim(dfData)

library(MASS)
fit.lda = lda(fCellID ~ ., data=dfData)
plot(fit.lda)

# prediction error on training data
p = predict(fit.lda)

dfPred.train.lda = data.frame(round(p$posterior, 3), actual=dfData$fCellID)

### predict on the new data from single cells
dfData.new = data.frame(scale(t(mCounts.SC)))
dfData.new = dfData.new[,paste('X', cvTopGenes, sep='')]
dim(dfData.new)

p = predict(fit.lda, newdata = dfData.new)

dfSingleCellPred.lda = data.frame(p$posterior)

n = make.names(paste('Predictions for single cell classes using random forest and lda louisa rds'))
n2 = paste0('~/Data/MetaData/', n)

save(dfSingleCellPred.lda, file=n2)

## comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='Predictions for single cell classes using random forest and lda for louisa j project')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

## display bar plots for both data sets

plot.bar = function(mBar, title='', cols){
  # get the median to plot
  #p.old = par(mar=c(6,3,2,2)+0.1)
  l = barplot(mBar, beside=F, xaxt='n', main=title, col=cols)
  axis(side = 1, l, labels=F)
  text(l, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=colnames(mBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  #par(p.old)
}

## make an overlay barplot
n = rainbow(6)
m = t(as.matrix(dfSingleCellPred)); m = m+0.001
plot.bar(t(as.matrix(m[1,])), 'Classification of single cells', cols = n[1])
barplot(m[2,], col = n[2], yaxt='n', xaxt='n', add=T)
barplot(m[3,], col = n[3], yaxt='n', xaxt='n', add=T)
barplot(m[4,], col = n[4], yaxt='n', xaxt='n', add=T)
barplot(m[5,], col = n[5], yaxt='n', xaxt='n', add=T)
barplot(m[6,], col = n[6], yaxt='n', xaxt='n', add=T)
legend('topright', legend = rownames(m), fill = n, cex = 0.7)

par(mfrow=c(2,1))

m = t(as.matrix(dfSingleCellPred)); m = m+0.001
cs = colSums(m)
m = sweep(m, 2, cs, '/')

plot.bar(m, 'binomial', rainbow(6))

m = t(as.matrix(dfSingleCellPred.lda)); m = m+0.001
cs = colSums(m)
m = sweep(m, 2, cs, '/')

plot.bar(m, 'lda', rainbow(6))


