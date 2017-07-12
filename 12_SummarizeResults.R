# File: 12_SummarizeResults.R
# Auth: umar.niazi@kcl.as.uk
# DESC: Create a summary report of the cell classification and assembled sequence
# Date: 24/04/2017


source('header.R')

library(Biostrings)
library(IRanges)
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')
## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12)')
dfSample.db = dbGetQuery(db, q)
dfSample.db
# close connection after getting data
dbDisconnect(db)
## load the count matrix, bam files object pre and post removing duplicates, 
## assembled sequences, and sample annotation with receptor type information
n = paste0(dfSample.db$location[dfSample.db$id==40], dfSample.db$name[dfSample.db$id==40])
load(n)
n = paste0(dfSample.db$location[dfSample.db$id==39], dfSample.db$name[dfSample.db$id==39])
load(n)
n = paste0(dfSample.db$location[dfSample.db$id==43], dfSample.db$name[dfSample.db$id==43])
load(n)
## comment out older analysis using UCSC database
# n = paste0(dfSample$location[dfSample$id==42], dfSample$name[dfSample$id==42])
# dfSample = read.csv(n, header=T, stringsAsFactors = F, row.names=1)
## new one using imtg database
n = paste0(dfSample.db$location[dfSample.db$id==44], dfSample.db$name[dfSample.db$id==44])
dfSample = read.csv(n, header=T, stringsAsFactors = F, row.names=1)

## make count matrix and get ercc vs normal genes read counts
## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = gsub('^(Plate.+)_S.+', '\\1', colnames(mCounts))
# reorder the count matrix columns according to the order in samples table
i = match(dfSample$title, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample$title, colnames(mCounts))
i2 = grepl(pattern = '^ERCC', rownames(mCounts))
table(i2)
m1 = mCounts[i2,]
m2 = mCounts[!i2,]
m1 = colSums(m1)
m2 = colSums(m2)
mCounts = rbind(m1, m2)
rownames(mCounts) = c('ERCC', 'Genes')

dfSample$ERCC = round(mCounts['ERCC',]/1e6, 3)
dfSample$Genes = round(mCounts['Genes',]/1e6, 3)
## proprtion of reads aligned to genes vs ercc
cs = colSums(mCounts)
mCounts = sweep(mCounts, 2, cs, '/')
dfSample$Genes_Proportion = round(mCounts['Genes',], 3)

## how many reads aligned after duplicate removal
lAllBams$desc = NULL
lAllBams$meta = NULL
# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)
## order the read counts as the sample annotation
i = match(dfSample$title, names(iReadCount))
identical(dfSample$title, names(iReadCount)[i])
dfSample$ReadsAligned_rmDup = iReadCount[i]

## load the bam before duplicate removal
n = paste0(dfSample.db$location[dfSample.db$id==45], dfSample.db$name[dfSample.db$id==45])
load(n)

lAllBams$desc = NULL
lAllBams$meta = NULL
# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)
## order the read counts as the sample annotation
i = match(dfSample$title, names(iReadCount))
identical(dfSample$title, names(iReadCount)[i])
dfSample$ReadsAligned = iReadCount[i]

## perform logistic regression to see if 
## any relationship between these parameters and ig class assignment
p = apply(as.matrix(dfSample[,7:10]), 2, is.na)
p = apply(p, 1, function(x) any(x == T))
table(p)
df = data.frame(present=as.numeric(!p), dfSample[,c(11, 12, 13, 14, 15)] )
fit = glm(present ~ Genes_Proportion, data=df, family='binomial')
fit = glm(present ~ ReadsAligned_rmDup, data=df, family='binomial')
summary(fit)
df2 = data.frame(ReadsAligned_rmDup = seq(0, 8, length.out = 100))
p = predict(fit, df2, type='response')
plot(df2$ReadsAligned_rmDup, p, type='l', ylim=c(0,1))
points(df$ReadsAligned_rmDup, df$present)

########## add sequences to the table
identical(names(lSeqs), dfSample$title)

writeSeqs = function(s){
  paste('>',names(s), '\n', s, '\n', sep='')
}

cvSeqs = sapply(lSeqs, function(x){
  if (is.null(x)) return(NA)
  x2 = sapply(x, as.character)
  paste(writeSeqs(x2), collapse = '')
})

dfSample$Sequences = cvSeqs

## new section added for the assigned cell classes based on 
## script 18_classifyNanoStToSingleCell.R
##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
q = paste0('select Sample.id, Sample.group3 from Sample
           where Sample.id=', dfSample$id)
dfSample.db = lapply(q, function(x) dbGetQuery(db, x))
dfSample.db = do.call(rbind, dfSample.db)
# close connection after getting data
dbDisconnect(db)
identical(dfSample$id, dfSample.db$id)
dfSample$group3 = dfSample.db$group3

n = paste0('Results/Receptor_summary_imgt', make.names(date()), '.csv')
write.csv(dfSample, file=n)

### make a pca plot with the scater object and only for cells with 
### classes defined
library(scater)
n = paste0(dfSample.db$location[dfSample.db$id==46], dfSample.db$name[dfSample.db$id==46])
load(n)

dfSamples = dfSamples[dfSamples$id %in% oSce.F$dbID_Sample,]
# sanity check
identical(dfSamples$id, oSce.F$dbID_Sample)

mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')
# set the factor for colours and drop NA samples
fSamples = dfSamples$IGH
i = which(is.na(fSamples))
mCounts.s = mCounts.s[,-i]
fSamples = factor(fSamples[-i])
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)

# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)


plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized, Heavy Constant')
text(pr.out$x[,1:2], labels = dfSamples$title[-i], pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

### light constant
mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')
# set the factor for colours and drop NA samples
fSamples = dfSamples$IGL
i = which(is.na(fSamples))
mCounts.s = mCounts.s[,-i]
fSamples = factor(fSamples[-i])
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)

# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)


plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized, Light Constant')
text(pr.out$x[,1:2], labels = dfSamples$title[-i], pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## heavy variable
mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')
# set the factor for colours and drop NA samples
fSamples = dfSamples$IGHVar
i = which(is.na(fSamples))
mCounts.s = mCounts.s[,-i]
fSamples = factor(fSamples[-i])
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)

# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)


plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized, Heavy Variable')
text(pr.out$x[,1:2], labels = dfSamples$title[-i], pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.6, ncol=2)


## light variable
mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')
# set the factor for colours and drop NA samples
fSamples = dfSamples$IGLightVar
i = which(is.na(fSamples))
mCounts.s = mCounts.s[,-i]
fSamples = factor(fSamples[-i])
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)

# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)


plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized, Light Variable')
text(pr.out$x[,1:2], labels = dfSamples$title[-i], pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.6, ncol=2)
