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

## normalize the data matrix
oNanoNorm.n = NanoStringNorm(oNanoNorm, CodeCount = 'geo.mean', Background = 'mean.2sd', 
                             SampleContent = 'housekeeping.geo.mean', round.values = T, 
                             take.log=T)


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

