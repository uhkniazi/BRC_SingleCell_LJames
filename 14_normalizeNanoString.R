# File: 14_normalizeNanoString
# Auth: umar.niazi@kcl.ac.uk
# DESC: normalize and some basic quality checks on the .rcc files from nanostring
# Date: 17/05/2017


## set variables and source libraries
source('header.R')
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did2
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 13) AND (File.idSample = Sample.id)')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### set working directory to appropriate location with bam files
setwd(gcRemoteDir)
setwd('nanoString/rcc/')
csFiles = list.files('.', pattern = '.RCC')
# check if these files match the file names in database
identical(dfSample$name, csFiles)

## open the nanostring norm library
library(NanoStringNorm)
oNanoNorm = read.markup.RCC()
oNanoNorm
## some basic stats of the data
dim(oNanoNorm$x)
dim(oNanoNorm$header)
rownames(oNanoNorm$header)
colnames(oNanoNorm$header)
identical(dfSample$title, colnames(oNanoNorm$header))
## add some sample covariates
df = oNanoNorm$header
dfSample$BatchCartridge = as.character(df['CartridgeID',])
dfSample$BatchDate = as.character(df['Date',])
dfSample$BatchLane = as.character(df['lane.id',])
dfSample$BindingDensity = as.numeric(df['BindingDensity',])

##### download the diagnostics library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')
## normalize the data matrix in 2 different ways and perform diagnostics
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean.2sd', 
                             SampleContent = 'housekeeping.geo.mean', round.values = T, 
                             take.log=T)

## use the CDiagnostics class to check some diagnostic plots
mCounts = oNanoNorm.n$normalized.data
table(mCounts$Code.Class)
mCounts = mCounts[mCounts$Code.Class == 'Endogenous',]
mCounts = mCounts[,-c(1,2,3)]
mCounts = as.matrix(mCounts)
min(mCounts)
max(mCounts)

## create the object and set a batch for plotting
oDiag = CDiagnosticPlots(mCounts, 'Normalised method 1')
fBatch = factor(dfSample$BatchCartridge)
setwd(gcswd)
pdf('Temp/norm_pos.geo.mean_bkgrnd.mean.2sd_house.pdf')
par(mfrow=c(2,2))
boxplot.median.summary(oDiag, fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
plot.missing.summary(oDiag, fBatch)
par(mfrow=c(1,1))
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 1)
l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l$PCA.scaleVariables = F
l$HC.scaleVaribles = F

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.2, fBatch, labels_cex = 1)
oDiag = calculateExtremeValues(oDiag)
m = mGetExtremeValues(oDiag)
## which sample has the most extreme values
apply(m, 2, function(x) sum(x > 0))
## which variable was extreme most times
v = apply(m, 1, function(x) sum(x > 0))
v[which(v > 0)]
mCounts[v>0,]
dev.off(dev.cur())


## diagnostic plots
setwd(gcswd)
pdf('Results/nanoStringDiag.pdf')
Plot.NanoStringNorm(oNanoNorm.n, plot.type='all')
dev.off(dev.cur())

lNanoString = list('normalizedData' = oNanoNorm.n, 'rawData'=oNanoNorm, 'metaData'=dfSample)

n = make.names(paste('NanoString Data for louisa single cell project rds'))
lNanoString$desc = paste('NanoString Data with sample annotations for louisa single cell project', date())
n2 = paste0('~/Data/MetaData/', n)
# save(lNanoString, file=n2)
# 
# # comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did2, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='NanoString Data with sample annotations for louisa single cell project')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

#### normalize using a different parameter
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean.2sd', 
                             SampleContent = 'top.geo.mean', round.values = T, 
                             take.log=T)

## use the CDiagnostics class to check some diagnostic plots
mCounts = oNanoNorm.n$normalized.data
table(mCounts$Code.Class)
mCounts = mCounts[mCounts$Code.Class == 'Endogenous',]
mCounts = mCounts[,-c(1,2,3)]
mCounts = as.matrix(mCounts)
min(mCounts)
max(mCounts)

## create the object and set a batch for plotting
oDiag = CDiagnosticPlots(mCounts, 'Normalised method 2')
fBatch = factor(dfSample$BatchCartridge)
setwd(gcswd)
pdf('Temp/norm_pos.geo.mean_bkgrnd.mean.2sd_top.geo.mean.pdf')
par(mfrow=c(2,2))
boxplot.median.summary(oDiag, fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
plot.missing.summary(oDiag, fBatch)
par(mfrow=c(1,1))
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 1)
l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l$PCA.scaleVariables = F
l$HC.scaleVaribles = F

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.2, fBatch, labels_cex = 1)
oDiag = calculateExtremeValues(oDiag)
m = mGetExtremeValues(oDiag)
## which sample has the most extreme values
apply(m, 2, function(x) sum(x > 0))
## which variable was extreme most times
v = apply(m, 1, function(x) sum(x > 0))
v[which(v > 0)]
mCounts[v>0,]
dev.off(dev.cur())

#### normalize using a different 3rd parameter
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean', 
                             SampleContent = 'top.geo.mean', round.values = T, 
                             take.log=T)

## use the CDiagnostics class to check some diagnostic plots
mCounts = oNanoNorm.n$normalized.data
table(mCounts$Code.Class)
mCounts = mCounts[mCounts$Code.Class == 'Endogenous',]
mCounts = mCounts[,-c(1,2,3)]
mCounts = as.matrix(mCounts)
min(mCounts)
max(mCounts)

## create the object and set a batch for plotting
oDiag = CDiagnosticPlots(mCounts, 'Normalised method 3')
fBatch = factor(dfSample$BatchCartridge)
setwd(gcswd)
pdf('Temp/norm_pos.geo.mean_bkgrnd.mean_top.geo.mean.pdf')
par(mfrow=c(2,2))
boxplot.median.summary(oDiag, fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
plot.missing.summary(oDiag, fBatch)
par(mfrow=c(1,1))
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 1)
l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l$PCA.scaleVariables = F
l$HC.scaleVaribles = F

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.2, fBatch, labels_cex = 1)
oDiag = calculateExtremeValues(oDiag)
m = mGetExtremeValues(oDiag)
## which sample has the most extreme values
apply(m, 2, function(x) sum(x > 0))
## which variable was extreme most times
v = apply(m, 1, function(x) sum(x > 0))
v[which(v > 0)]
mCounts[v>0,]
dev.off(dev.cur())


#### normalize using a different 4th parameter setting
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean', 
                             SampleContent = 'housekeeping.geo.mean', round.values = T, 
                             take.log=T)

## use the CDiagnostics class to check some diagnostic plots
mCounts = oNanoNorm.n$normalized.data
table(mCounts$Code.Class)
mCounts = mCounts[mCounts$Code.Class == 'Endogenous',]
mCounts = mCounts[,-c(1,2,3)]
mCounts = as.matrix(mCounts)
min(mCounts)
max(mCounts)

## create the object and set a batch for plotting
oDiag = CDiagnosticPlots(mCounts, 'Normalised method 4')
fBatch = factor(dfSample$BatchCartridge)
setwd(gcswd)
pdf('Temp/norm_pos.geo.mean_bkgrnd.mean_housekeeping.geo.mean.pdf')
par(mfrow=c(2,2))
boxplot.median.summary(oDiag, fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
plot.missing.summary(oDiag, fBatch)
par(mfrow=c(1,1))
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 1)
l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l$PCA.scaleVariables = F
l$HC.scaleVaribles = F

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.2, fBatch, labels_cex = 1)
oDiag = calculateExtremeValues(oDiag)
m = mGetExtremeValues(oDiag)
## which sample has the most extreme values
apply(m, 2, function(x) sum(x > 0))
## which variable was extreme most times
v = apply(m, 1, function(x) sum(x > 0))
v[which(v > 0)]
mCounts[v>0,]
dev.off(dev.cur())

################ looking at all the diagnostics, method 2 appears to be good
#### normalize using strong background correction but top 75 expressed genes for normalizing
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean.2sd', 
                             SampleContent = 'top.geo.mean', round.values = T, 
                             take.log=T)


## diagnostic plots
setwd(gcswd)
pdf('Results/nanoStringDiag_2.pdf')
Plot.NanoStringNorm(oNanoNorm.n, plot.type='all')
dev.off(dev.cur())

lNanoString = list('normalizedData' = oNanoNorm.n, 'rawData'=oNanoNorm, 'metaData'=dfSample)

n = make.names(paste('NanoString Data for louisa single cell project normalized top geo mean rds'))
lNanoString$desc = paste('NanoString Data with sample annotations for louisa single cell project top geo mean normalised', date())
n2 = paste0('~/Data/MetaData/', n)
save(lNanoString, file=n2)

# comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did2, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='NanoString Data with sample annotations for louisa single cell project normalized using top expressed genes')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)