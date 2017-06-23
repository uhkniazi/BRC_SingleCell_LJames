# File: 17_nanoStringDE.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: perform DE on the genes doing an anova, cluster and select genes
# Date: 21/06/2017


## set variables and source libraries
source('header.R')

library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

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
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
           where (Sample.idData = 13) AND (File.idSample = Sample.id)')
dfSample.names = dbGetQuery(db, q)
dim(dfSample.names)
dfSample.names
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
mCounts.norm = as.matrix(mCounts.norm[,-c(1:3)])

mCounts = mCounts.norm[rowMeans(mCounts.norm) > 1,]

dfSample.names = lNanoString$metaData
identical(dfSample.names$title, colnames(mCounts.norm))



## create groups and perform anova after fitting random effects model
fGroups = factor(dfSample.names$group2)
fPat = paste0(dfSample.names$group1, ':', gsub('louisajames', '', dfSample.names$BatchCartridge))
fPat = factor(fPat)
table(fPat)

library(lmerTest)
library(parallel)
getAnovaPvalue = function(iIndex){
  x = mCounts[iIndex,]
  df = data.frame(x, fGroups, pat=fPat)
  fit = lmerTest::lmer(x ~ fGroups + (1 | pat), data=df)
  return(anova(fit)$`Pr(>F)`)
}

lPvals = mclapply(1:nrow(mCounts), function(iIndexSub) {
  tryCatch(getAnovaPvalue(iIndexSub), error=function(e) NULL)
})

names(lPvals) = rownames(mCounts)
ivPvals = unlist(lPvals)
ivPvals.adj = p.adjust(ivPvals, 'bonf')
summary(ivPvals.adj); hist(ivPvals.adj)
# use the top genes to cluster that samples
fTop = ivPvals.adj < 1e-2
table(fTop)

oDiag = CDiagnosticPlots(mCounts[fTop,], 'Top 102')
l = CDiagnosticPlotsGetParameters(oDiag)
# don't centre and scale subjects, but do it in variable space only
l$PCA.scaleSubjects=F; l$HC.scaleSubjects = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

boxplot.median.summary(oDiag, fGroups, legend.pos = 'topright')
plot.mean.summary(oDiag, fGroups)
plot.sigma.summary(oDiag, fGroups)
plot.missing.summary(oDiag, fGroups)

plot.PCA(oDiag, fGroups)
plot.PCA(oDiag.norm, fBatch.2)
#plot.PCA(oDiag.raw, fBatch.2)

plot.dendogram(oDiag.norm, fBatch, labels_cex = 0.8, cex.main=0.8)
plot.dendogram(oDiag.norm, fBatch.2, labels_cex = 0.8, cex.main=0.8)
#plot.dendogram(oDiag.raw, fBatch, labels_cex = 0.8, cex.main=0.8)

