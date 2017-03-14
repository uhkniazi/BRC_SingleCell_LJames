# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 10/03/2017


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

# setwd(gcRemoteDir)
setwd('Data_external/')
setwd('Fastq/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz')

# each sample has 2 files 
fSplit = gsub('^(Plate\\d-C\\d+)_.+', '\\1', cvFiles)

lFiles = split(cvFiles, fSplit)

## create the entry for samples
cSampleCol
dfSamples = data.frame(idProject=g_pid, idData=g_did, title=unique(fSplit), description='Single Cell Sequencing rna-seq dataset')
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 12;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
# get the names of the samples
temp = lapply(dfSamples$title, function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[dfSamples$title == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
# dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)
dbDisconnect(db)

##### adding some additional information for covariates and groupings
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 12;'))

# group 1
group1 = gsub('^Plate(\\d)-C\\d+', '\\1', dfSamples$title)
group1 = paste('P', group1, sep='')

# group 2
group2 = gsub('^Plate\\d-C(\\d+)', '\\1', dfSamples$title)
group2 = paste('Capture Site', group2)

dfSamples$group1 = group1
dfSamples$group2 = group2

## create an update statement for each row
queries = paste('Update Sample Set group1="', dfSamples$group1, '", group2="', dfSamples$group2, '" where Sample.id=', dfSamples$id, ';', sep='')

#sapply(queries, function(x) dbSendQuery(db, x))

dbDisconnect(db)
