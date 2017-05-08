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
write.csv(dfSample, file='Results/Receptor_summary_imtg.csv')




