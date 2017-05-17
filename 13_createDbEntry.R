# File: 13_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries for the nanoString data set
# Date: 17/05/2017


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

setwd(gcRemoteDir)
setwd('nanoString/rcc/')

# list the files
cvFiles = list.files(pattern = '.RCC')

gr = gsub('.RCC', '', cvFiles)
## create groupings for the samples
fGroups = strsplit(gr, '_')
fGroups = do.call(rbind, fGroups)
## create the entry for samples
cSampleCol
dfSamples = data.frame(idProject=g_pid, idData=g_did2, title=gr, description='NanoString data set from Louisa James for Tonsil B Cells RNA',
                       group1=fGroups[,1], group2=fGroups[,2])
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did2
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 13;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
# get the names of the samples
temp = lapply(dfSamples$title, function(x){
  # get the file names
  df = data.frame(name=paste(x, 'RCC', sep='.'), type='.RCC', idSample=dfSamples[dfSamples$title == x, 'id'], group1='NanoString raw data files')
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
# dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)
dbDisconnect(db)

