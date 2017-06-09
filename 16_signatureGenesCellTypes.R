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
n = findCorrelation((mCor), cutoff = 0.7, names=T)
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



dfSample.names = lNanoString$metaData
identical(dfSample.names$title, colnames(mCounts.norm))
identical(dfSample.names$title, colnames(mCounts.raw))