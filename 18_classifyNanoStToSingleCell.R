# File: 18_classifyNanoStToSingleCell.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: identify signature genes to classify cells
# Date: 26/06/2017


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

## load the gene list from previous script
load('Results/topSelected.rds')

# load the single cell normalised count matrix
mCounts.SC = exprs(oSce.F)

# load the nano string count matrix
names(lNanoString)
lNanoString$desc

mCounts.norm = lNanoString$normalizedData$normalized.data
mCounts.norm = mCounts.norm[mCounts.norm$Code.Class == 'Endogenous', ]
rownames(mCounts.norm) = mCounts.norm$Name

## subset this matrix based on the top genes selected in
## previous section
mCounts.norm = mCounts.norm[cvTopSelected,]


## match the names between the 2 data sets 
## i.e. symbols in nanostring and enterez id in single cell
dfSymbols = AnnotationDbi::select(org.Hs.eg.db, rownames(mCounts.norm), 
                                  columns = 'ENTREZID', keytype = 'SYMBOL')
dfSymbols = na.omit(dfSymbols)
# how many match b/w the 2 data sets
table(dfSymbols$ENTREZID %in% rownames(mCounts.SC))
# 
# FALSE  TRUE 
# 12     87 
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
if(length(i) > 0) mCounts = mCounts[-i,]
dim(mCounts)
mCounts = scale(t(mCounts))
head(apply(mCounts, 2, sd))

## create 2 groups
fG = as.character(lNanoString$metaData$group2)
fG[fG %in% c('GC', 'PB')] = 'GC-PB' 
fG[!(fG %in% 'GC-PB')] = 'Others' 
lNanoString$metaData$fCellID = factor(fG)

## prepare test data i.e. single cell data
mCounts.test = scale(t(mCounts.SC))
head(apply(mCounts.test, 2, var))
mCounts.test = mCounts.test[,colnames(mCounts)]
dim(mCounts.test)

##################### try with CCrossvalidation library
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## setup the data frames
dfData = data.frame(mCounts)
## add the cell type id
fCellID = lNanoString$metaData$fCellID

## create test vector
set.seed(123)
test = sample(1:nrow(dfData), size = 4, replace = F)
table(fCellID[test])

## perform nested random forest on test set
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData, fCellID, boot.num = 100)
# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:60]
cvTopGenes = gsub('X', '', cvTopGenes)
mCounts = mCounts[,cvTopGenes]
dim(mCounts)

## remove correlated variables first
## find correlated variables
mCor = cor(mCounts, use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.8, names=T)
i = which(colnames(mCounts) %in% n)
cvTopGenes = colnames(mCounts)[-i]

mCounts = mCounts[,cvTopGenes]
rm(mCor)
gc()

# use the top genes to find top combinations of genes
## setup the data frames
dfData = data.frame(mCounts)
## add the cell type id
fCellID = lNanoString$metaData$fCellID
dim(dfData)
oVar.sub = CVariableSelection.ReduceModel(dfData, fCellID, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

## k fold nested cross validation with various variable combinations
#par(mfrow=c(2,2))

cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 1)
oCV = CCrossValidation.LDA(test.dat = data.frame(gene=dfData[, cvTopGenes.sub]),
                           train.dat = data.frame(gene=dfData[, cvTopGenes.sub]),
                           test.groups = fCellID,
                           train.groups = fCellID, 
                           level.predict = 'GC-PB', boot.num = 5, k.fold = 3)
plot.cv.performance(oCV)

# try models of various sizes with CV
for (i in 2:4){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = dfData[, cvTopGenes.sub]
  
  dfData.test = dfData[test, cvTopGenes.sub]
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, 
                             test.groups = fCellID[test],
                             train.groups = fCellID, 
                             level.predict = 'GC-PB', boot.num = 5, k.fold = 3)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  x = getAUCVector(oCV)
  print(paste('Variable Count', i))
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

#### use a 2 variable model
dfData = data.frame(mCounts)
## add the cell type id

cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 2)
dfData = dfData[,cvTopGenes.sub]
dfData$fCellID = lNanoString$metaData$fCellID
dim(dfData)

fit.lda = lda(fCellID ~ ., data=dfData)
plot(fit.lda)

# prediction error on training data
p = predict(fit.lda)

dfPred.train.lda = data.frame(round(p$posterior, 3), actual=dfData$fCellID)

### predict on the new data from single cells
dfData.new = data.frame(mCounts.test)
dfData.new = dfData.new[,cvTopGenes.sub]
dim(dfData.new)

p = predict(fit.lda, newdata = dfData.new)

dfSingleCellPred.lda = data.frame(p$posterior)




############### follow the regression based approach
########## perform binomial regression on each category 
dfData = data.frame(mCounts)
dfData = dfData[,cvTopGenes.sub]
## add the cell type id
dfData$fCellID = lNanoString$metaData$fCellID
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'GC-PB'] = 1
fGroups[dfData$fCellID != 'GC-PB'] = 0

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

## calculate coefficients for the top variables and model size
f = paste('fCellID ~ ', paste(cvTopGenes.sub, collapse='+'), collapse = ' ')
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

## new data, the same training data
dfData.new = data.frame(mCounts)[,cvTopGenes.sub]
## perform prediction on this
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,names(post)[-1]]))
colnames(X) = names(post)
lData = list(mModMatrix=X)
p = mypred(post, lData)

dfPred.train = data.frame(round(p, 3), actual=dfData$fCellID)

## a 2 variable model appears to be generally good enough 
## use this 2 variable model to make predictions in single cell data
dfData.new = data.frame(mCounts.test)
dfData.new = dfData.new[,cvTopGenes.sub]
dim(dfData.new)
## perform prediction on this
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,cvTopGenes.sub]))
colnames(X) = names(post)
lData = list(mModMatrix=X)
p = mypred(post, lData)

dfSingleCellPred.bin = data.frame(GC_PB=round(p, 3), Others=round(1-p, 3))

## display bar plots for both data sets

plot.bar = function(mBar, title='', cols, ylab='Predicted Probability of Cell Type'){
  # get the median to plot
  #p.old = par(mar=c(6,3,2,2)+0.1)
  l = barplot(mBar, beside=F, xaxt='n', main=title, col=cols, ylab=ylab)
  axis(side = 1, l, labels=F)
  text(l, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=colnames(mBar), srt=45, adj=1, xpd=TRUE, cex=0.8)
  #par(p.old)
}

## make an overlay barplot
n = rainbow(2)
m = t(as.matrix(dfSingleCellPred.bin)); m = m+0.001
plot.bar(t(as.matrix(m[1,])), 'Classification of single cells - bin 2 var', cols = n[1])
barplot(m[2,], col = n[2], yaxt='n', xaxt='n', add=T)
legend('topright', legend = rownames(m), fill = n, cex = 0.7)

## classify the cells proportions
m = as.matrix(dfSingleCellPred.lda)
f = apply(m, 1, which.max)
#f2 = apply(m, 1, function(x) any(x > 0.80))
#f[!f2] = '7'
## convert to factor
fCellTypes = factor(f, labels = c(colnames(m)))
table(fCellTypes)


