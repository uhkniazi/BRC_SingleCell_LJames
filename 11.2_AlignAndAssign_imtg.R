# File: 11.2_AlignAndAssign_imtg.R
# Auth: umar.niazi@kcl.as.uk
# DESC: Align the assembled sequences for each sample with reference to assign a type to each cell
#       while using the imgt database
# Date: 03/05/2017


source('header.R')

library(Biostrings)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(downloader)
# NAME: f_oDNAStringSetConvertPWAMatrix
# ARGS: pwa.matrix = a pairwise alignment object in matrix form, generated via a call like
#       as.matrix(pairwiseAlignment(patten, subject))
# DESC: converts the sequences in the matrix to a DNAStringSet object which can be exported to a FASTA
#       file to be viewed in a viewer of choice e.g. jalview
# RETS: a DNAStringSet object with all sequences from the pwa matrix
f_oDNAStringSetConvertPWAMatrix= function(pwa.matrix){
  require(Biostrings)
  # sequence to return
  seq.export = DNAStringSet()
  # extract each sequence from the matrix
  for (i in 1:nrow(pwa.matrix)){
    s = DNAStringSet(pwa.matrix[i,])
    s = DNAStringSet(unlist(s))
    seq.export = append(seq.export, s)
  }
  # set names for sequences
  names(seq.export) = rownames(pwa.matrix)
  return(seq.export)
}

## database url for imtg database
url = 'http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP'
dir.create('Data_external', showWarnings = F)
csImgt = 'Data_external/imgt.fasta'
# download the reactome file if it doesnt exist
if (!file.exists(csImgt)) download(url, csImgt)
oSeqImgt = readBStringSet(csImgt)
## remove the non human sequences
f = grepl('Homo sapiens', names(oSeqImgt))
table(f)
oSeqImgt = oSeqImgt[f]
## types of sequences
fTypes = gsub(".+\\|(.+)\\|Homo .+", '\\1', names(oSeqImgt))
fTypes = gsub('^(.+)\\*\\w+', '\\1', fTypes)

############# heavy chain constant GRangesList coordinates are ready, extract sequence from the genome
fTypes.heavy = grep('IGHA|IGHE$|IGHG\\d|IGHD$|IGHM', fTypes)

oSeqIGHc = oSeqImgt[fTypes.heavy]

scoreHeavyConstant = function(oSeqPred){
  s = pairwiseAlignment(oSeqIGHc, subject=oSeqPred, type='local', scoreOnly=T)
  names(s) = names(oSeqIGHc)
  return(s)
}

############ repeat similar process for light chain constant regions
fTypes.light = grep('IGKC$|IGLC\\d+$', fTypes)
############# light chain constant GRangesList coordinates are ready, extract sequence from the genome
oSeqIGLc = oSeqImgt[fTypes.light] 

scoreLightConstant = function(oSeqPred){
  s = pairwiseAlignment(oSeqIGLc, subject=oSeqPred, type='local', scoreOnly=T)
  names(s) = names(oSeqIGLc)
  return(s)
}

##########################################################################################
## load the sequences for each BASIC result
setwd('Results/Results')
cmd = paste('find Plate* -iname "result.fasta"')
csFiles = system(cmd, intern = T)

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select id, title, group1, group2, group3 from Sample where idData=12;")
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)

## split the file names to map them back to the original annotations in the sample table from database
csFiles.Plate = gsub('(Plate1|2).+', '\\1', x = csFiles)
csFiles.Capture = gsub('Plate[1|2](C\\d*).+', '\\1', x = csFiles)
dfFile.loc = data.frame(loc=csFiles, title=paste(csFiles.Plate, csFiles.Capture, sep='-'))
i = match(dfSample$title, dfFile.loc$title)
dfSample$location = dfFile.loc$loc[i]
dfSample$location = as.character(dfSample$location)
dfSample$IGH = NA
dfSample$IGL = NA
### go through each sample and read the heavy chain
for (i in 1:nrow(dfSample)){
  if (is.na(dfSample$location[i])) next;
  ## import the assembled sequence
  seq.as = import(dfSample$location[i])
  ### search for the word heavy chain and light chain
  ## for heavy
  seq.as.h = seq.as[grepl(pattern = 'heavy', names(seq.as))]
  if (length(seq.as.h) > 1) seq.as.h = seq.as.h[grepl(pattern = 'constant', names(seq.as.h))]
  ## error check
  if (length(seq.as.h) != 1) {warning(paste(dfSample$location[i], 'some error while looking for heavy constant region'));
  } else {
    ## assign the heavy chain type
    sc = scoreHeavyConstant(seq.as.h)
    dfSample$IGH[i] = names(which.max(sc))
  }
  ## for light
  seq.as.l = seq.as[grepl(pattern = 'light', names(seq.as))]
  if (length(seq.as.l) > 1) seq.as.l = seq.as.l[grepl(pattern = 'constant', names(seq.as.l))]
  ## error check
  if (length(seq.as.l) != 1) {warning(paste(dfSample$location[i], 'some error while looking for light constant region'));
  } else {
    ## assign the light chain type
    sc = scoreLightConstant(seq.as.l)
    dfSample$IGL[i] = names(which.max(sc))
  }
}

setwd(gcswd)
temp = gsub(".+\\|(.+)\\|Homo .+", '\\1', dfSample$IGH)
dfSample$IGH = gsub('^(.+)\\*\\w+', '\\1', temp)
temp = gsub(".+\\|(.+)\\|Homo .+", '\\1', dfSample$IGL)
dfSample$IGL = gsub('^(.+)\\*\\w+', '\\1', temp)

######################################################################################################
####### construct the variable regions
fTypes.heavy = grep('IGHV', fTypes)
fTypes.light = grep('IGLV|IGKV', fTypes)
##### get the sequence for the heavy and light chain regions
oSeqIGHc = oSeqImgt[fTypes.heavy]
oSeqIGLc = oSeqImgt[fTypes.light]

alignHeavyVariable = function(oSeqPred, title){
  pwa = pairwiseAlignment(oSeqIGHc, subject=oSeqPred, type='local')
  m = as.matrix(pwa)
  c = unlist(strsplit(toString(oSeqPred), split = ''))
  m = rbind(m, c)
  se = f_oDNAStringSetConvertPWAMatrix(m)
  names(se) = c(names(oSeqIGHc), title)
  return(se)
}

alignLightVariable = function(oSeqPred, title){
  pwa = pairwiseAlignment(oSeqIGLc, subject=oSeqPred, type='local')
  m = as.matrix(pwa)
  c = unlist(strsplit(toString(oSeqPred), split = ''))
  m = rbind(m, c)
  se = f_oDNAStringSetConvertPWAMatrix(m)
  names(se) = c(names(oSeqIGLc), title)
  return(se)
}

scoreHeavyVariable = function(oSeqPred){
  s = pairwiseAlignment(oSeqIGHc, subject=oSeqPred, type='local', scoreOnly=T)
  names(s) = names(oSeqIGHc)
  return(s)
}

scoreLightVariable = function(oSeqPred){
  s = pairwiseAlignment(oSeqIGLc, subject=oSeqPred, type='local', scoreOnly=T)
  names(s) = names(oSeqIGLc)
  return(s)
}

setwd('Results/Results/')
dfSample$IGHVar = NA
dfSample$IGLightVar = NA
### go through each sample and read the heavy chain
for (i in 1:nrow(dfSample)){
  if (is.na(dfSample$location[i])) next;
  ## import the assembled sequence
  seq.as = import(dfSample$location[i])
  ### search for the word heavy chain and light chain
  ## for heavy
  seq.as.h = seq.as[grepl(pattern = 'heavy', names(seq.as))]
  if (length(seq.as.h) > 1) seq.as.h = seq.as.h[grepl(pattern = 'variable', names(seq.as.h))]
  ## error check
  if (length(seq.as.h) != 1) {warning(paste(dfSample$location[i], 'some error while looking for heavy variable region'));
  } else {
    ## assign the heavy chain type
    sc = scoreHeavyVariable(seq.as.h)
    dfSample$IGHVar[i] = names(which.max(sc))
  }
  ## for light
  seq.as.l = seq.as[grepl(pattern = 'light', names(seq.as))]
  if (length(seq.as.l) > 1) seq.as.l = seq.as.l[grepl(pattern = 'variable', names(seq.as.l))]
  ## error check
  if (length(seq.as.l) != 1) {warning(paste(dfSample$location[i], 'some error while looking for light variable region'));
  } else {
    ## assign the light chain type
    sc = scoreLightVariable(seq.as.l)
    dfSample$IGLightVar[i] = names(which.max(sc))
  }
}

setwd(gcswd)
temp = gsub(".+\\|(.+)\\|Homo .+", '\\1', dfSample$IGHVar)
dfSample$IGHVar = gsub('^(.+)\\*\\w+', '\\1', temp)
temp = gsub(".+\\|(.+)\\|Homo .+", '\\1', dfSample$IGLightVar)
dfSample$IGLightVar = gsub('^(.+)\\*\\w+', '\\1', temp)

# ## load and store each sequence in a biostrings object 
# setwd('Results/Results/')
# ### go through each sample and read the sequences and store in a biostrings object
# lSeqs = lapply(1:nrow(dfSample), function(i) {
#   if (is.na(dfSample$location[i])) return(NULL);
#   ## import the assembled sequence
#   seq.as = import(dfSample$location[i])
#   return(seq.as)
# })
# 
# names(lSeqs) = as.character(dfSample$title)
# setwd(gcswd)

#################### write results to data base
# n = make.names(paste('dfSingleCellAnnotation imtg joana single cell chain annotation csv'))
# n2 = paste0('~/Data/MetaData/', n)
# write.csv(dfSample, file=n2)
# 
# ## note: comment out as this entry has been made in db
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='csv', location='~/Data/MetaData/', comment='joana single cell sequencing project annotations for the heavy and light chains for each cell using the imtg database')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)

# ##
# n = make.names(paste('list of DNAStringSet objects joana b cell heavy light chains rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(lSeqs, file=n2)
# 
# ## note: comment out as this entry has been made in db
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='joana single cell sequencing project list of DNAStringSet objects for heavy and light chains')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)




########### if we want to generate fasta alignment files then follow this
### go through each sample and read the heavy chain and light chain
# for (i in 1:nrow(dfSample)){
#   if (is.na(dfSample$location[i])) next;
#   ## import the assembled sequence
#   setwd('Results/Results')
#   seq.as = import(dfSample$location[i])
#   ### search for the word heavy chain and light chain
#   ## for heavy
#   seq.as.h = seq.as[grepl(pattern = 'heavy', names(seq.as))]
#   if (length(seq.as.h) > 1) seq.as.h = seq.as.h[grepl(pattern = 'variable', names(seq.as.h))]
#   ## error check
#   if (length(seq.as.h) != 1) {warning(paste(dfSample$location[i], 'some error while looking for heavy variable region'));
#   } else {
#     ## align the heavy chain type
#     sc = alignHeavyVariable(seq.as.h, dfSample$title[i])
#     setwd(gcswd)
#     export(sc, paste0('Results/', dfSample$title[i], 'HeavyVariable.fasta'))
#   }
#   ## for light
#   seq.as.l = seq.as[grepl(pattern = 'light', names(seq.as))]
#   if (length(seq.as.l) > 1) seq.as.l = seq.as.l[grepl(pattern = 'variable', names(seq.as.l))]
#   ## error check
#   if (length(seq.as.l) != 1) {warning(paste(dfSample$location[i], 'some error while looking for light variable region'));
#   } else {
#     ## assign the light chain type
#     sc = alignLightVariable(seq.as.l, dfSample$title[i])
#     setwd(gcswd)
#     export(sc, paste0('Results/', dfSample$title[i], 'LightVariable.fasta'))
#   }
# }

