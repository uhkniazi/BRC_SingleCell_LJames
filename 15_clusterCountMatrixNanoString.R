# File: 15_clusterCountMatrixNanoString.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: cluster the samples based on count table for single cell data
# Date: 05/06/2017


## set variables and source libraries
source('header.R')

library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

######################### functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}

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
mCounts.raw = lNanoString$rawData$x
mCounts.norm = lNanoString$normalizedData$normalized.data
mCounts.raw = mCounts.raw[mCounts.raw$CodeClass == 'Endogenous', ]
rownames(mCounts.raw) = mCounts.raw$Name
mCounts.raw = as.matrix(mCounts.raw[,-c(1:3)])

mCounts.norm = mCounts.norm[mCounts.norm$Code.Class == 'Endogenous', ]
rownames(mCounts.norm) = mCounts.norm$Name
mCounts.norm = as.matrix(mCounts.norm[,-c(1:3)])

dfSample.names = lNanoString$metaData
identical(dfSample.names$title, colnames(mCounts.norm))
identical(dfSample.names$title, colnames(mCounts.raw))

## create the normalized and raw matrices and perform diagnostics
oDiag.norm = CDiagnosticPlots(mCounts.norm, 'Normalised')
min(mCounts.raw)
#oDiag.raw = CDiagnosticPlots(log(mCounts.raw+1), 'Raw')

#par(mfrow=c(1,2))
## choose a few batches of choice
fBatch = dfSample.names$BatchCartridge
fBatch = gsub('louisajames', '', fBatch)
fBatch = factor(paste('C', fBatch, sep='-'))

fBatch.2 = factor(paste('P', dfSample.names$group1, sep='-'))

boxplot.median.summary(oDiag.norm, fBatch, legend.pos = 'topright')
boxplot.median.summary(oDiag.norm, fBatch.2, legend.pos = 'topright')
#boxplot.median.summary(oDiag.raw, fBatch, legend.pos = 'topright')

plot.mean.summary(oDiag.norm, fBatch)
plot.mean.summary(oDiag.norm, fBatch.2)
#plot.mean.summary(oDiag.raw, fBatch)

plot.sigma.summary(oDiag.norm, fBatch)
plot.sigma.summary(oDiag.norm, fBatch.2)
#plot.sigma.summary(oDiag.raw, fBatch)

plot.missing.summary(oDiag.norm, fBatch)
plot.missing.summary(oDiag.norm, fBatch.2)
#plot.missing.summary(oDiag.raw, fBatch)

plot.PCA(oDiag.norm, fBatch)
plot.PCA(oDiag.norm, fBatch.2)
#plot.PCA(oDiag.raw, fBatch.2)

plot.dendogram(oDiag.norm, fBatch, labels_cex = 0.8, cex.main=0.8)
plot.dendogram(oDiag.norm, fBatch.2, labels_cex = 0.8, cex.main=0.8)
#plot.dendogram(oDiag.raw, fBatch, labels_cex = 0.8, cex.main=0.8)

