# File: 11_AlignAndAssign.R
# Auth: umar.niazi@kcl.as.uk
# DESC: Align the assembled sequences for each sample with reference to assign a type to each cell
# Date: 05/04/2017


source('header.R')

library(Biostrings)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

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

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'tx', use.names=T)
length(oGRLgenes)

# load the coordinates for the constant regions
df = read.csv('IGH_coordinates.csv', header=T, stringsAsFactors = F)
oGRighC = GRanges(df$Chromosome, IRanges(df$Start, df$End), strand=ifelse(df$Strand == 'minus', '-', '+'))

# find the transcript positions for these regions
f = overlapsAny(oGRLgenes, oGRighC)
table(f)
# the overlaps are with 19 transcripts, some may be alternative versions?
oGRLgenes.o = oGRLgenes[f]
countOverlaps(oGRighC, oGRLgenes.o)
oGRLighC = subsetByOverlaps(oGRLgenes.o, oGRighC[1])
oGRLighC = append(oGRLighC, subsetByOverlaps(oGRLgenes.o, oGRighC[2]))
# 4 exons according to ucsc
c = subsetByOverlaps(oGRLgenes.o, oGRighC[3], type = 'end')
oGRLighC = append(oGRLighC, subsetByOverlaps(oGRLgenes.o, oGRighC[3], type = 'end'))
# 6 exons according to ucsc
c = subsetByOverlaps(oGRLgenes.o, oGRighC[4])
oGRLighC = append(oGRLighC, c[1])
# 4 exons according to ucsc
c = subsetByOverlaps(oGRLgenes.o, oGRighC[5])
oGRLighC = append(oGRLighC, c[1])

c = subsetByOverlaps(oGRLgenes.o, oGRighC[6])
oGRLighC = append(oGRLighC, c)

c = subsetByOverlaps(oGRLgenes.o, oGRighC[7])
oGRLighC = append(oGRLighC, c[1])

c = subsetByOverlaps(oGRLgenes.o, oGRighC[8])
oGRLighC = append(oGRLighC, c[1])

c = subsetByOverlaps(oGRLgenes.o, oGRighC[9])
oGRLighC = append(oGRLighC, c[1])

names(oGRLighC) = df$Gene
############# heavy chain constant GRangesList coordinates are ready, extract sequence from the genome
oSeqIGHc = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRLighC)
## unlist each exon to make a contigous sequence
oSeqIGHc = sapply(oSeqIGHc, unlist)

scoreHeavyConstant = function(oSeqPred){
  s = pairwiseAlignment(oSeqIGHc, subject=oSeqPred, type='local', scoreOnly=T)
  names(s) = names(oSeqIGHc)
  return(s)
}

############ repeat similar process for light chain constant regions

df = read.csv('IGLambdaConstant_coordinates.csv', header=T, stringsAsFactors = F)
oGRiglC = GRanges(df$Chromosome, IRanges(df$Start, df$End), strand=ifelse(df$Strand == 'minus', '-', '+'))

# find the transcript positions for these regions
f = overlapsAny(oGRLgenes, oGRiglC)
table(f)
# the overlaps are with 18 transcripts, some may be alternative versions?
oGRLgenes.o = oGRLgenes[f]
countOverlaps(oGRiglC, oGRLgenes.o)
# 4 exons according to ucsc
c = subsetByOverlaps(oGRLgenes.o, oGRiglC[1])
# use the second one with the widest exons
oGRLigLightC = c[2]

oGRLigLightC = append(oGRLigLightC, subsetByOverlaps(oGRLgenes.o, oGRiglC[2]))
oGRLigLightC = append(oGRLigLightC, subsetByOverlaps(oGRLgenes.o, oGRiglC[3]))
oGRLigLightC = append(oGRLigLightC, subsetByOverlaps(oGRLgenes.o, oGRiglC[4]))
# 11 exons according to ucsc for last one
c = subsetByOverlaps(oGRLgenes.o, oGRiglC[5])
# use uc284ppo.1 with 5 exons and largest width
oGRLigLightC = append(oGRLigLightC, c[4])
# sanity check
countOverlaps(oGRiglC, oGRLigLightC)
countOverlaps(oGRighC, oGRLigLightC)
countOverlaps(oGRighC, oGRLighC)

names(oGRLigLightC) = df$Gene
############# light chain constant GRangesList coordinates are ready, extract sequence from the genome
oSeqIGLc = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRLigLightC)
## unlist each exon to make a contigous sequence
oSeqIGLc = sapply(oSeqIGLc, unlist)

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
######################################################################################################
####### construct the variable regions
chr14 = read.csv(gzfile('hg38_chr14_gencodecompV24.gz'), header=T, sep='\t', stringsAsFactors = F)
chr2 = read.csv(gzfile('hg38_chr2_gencodecompV24.gz'), header=T, sep='\t', stringsAsFactors = F)
chr22 = read.csv(gzfile('hg38_chr22_gencodecompV24.gz'), header=T, sep='\t', stringsAsFactors = F)

## coordinates for the heavy variables
l = grep('IGHV', chr14$name2); length(l)
chr = chr14[l,]
oGRighVar = GRanges(chr$chrom, IRanges(chr$txStart, chr$txEnd), strand=chr$strand) 
names(oGRighVar) = chr$name2

## coordinates for the light variables
l = grep('IGKV', chr2$name2); length(l)
chr = chr2[l,]
l = grep('IGLV', chr22$name2); length(l)
chr = rbind(chr, chr22[l,])

oGRigLightVar = GRanges(chr$chrom, IRanges(chr$txStart, chr$txEnd), strand=chr$strand) 
names(oGRigLightVar) = chr$name2

##### get the sequence for the heavy and light chain regions
oSeqIGHc = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRighVar)
oSeqIGLc = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRigLightVar)

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
  if (length(seq.as.l) != 1) {warning(paste(dfSample$location[i], 'some error while looking for light constant region'));
  } else {
    ## assign the light chain type
    sc = scoreLightVariable(seq.as.l)
    dfSample$IGLightVar[i] = names(which.max(sc))
  }
}

setwd(gcswd)

## load and store each sequence in a biostrings object 
setwd('Results/Results/')
### go through each sample and read the sequences and store in a biostrings object
lSeqs = lapply(1:nrow(dfSample), function(i) {
  if (is.na(dfSample$location[i])) return(NULL);
  ## import the assembled sequence
  seq.as = import(dfSample$location[i])
  return(seq.as)
})

names(lSeqs) = as.character(dfSample$title)
setwd(gcswd)

#################### write results to data base
# n = make.names(paste('dfSingleCellAnnotation joana single cell chain annotation csv'))
# n2 = paste0('~/Data/MetaData/', n)
# write.csv(dfSample, file=n2)
# 
# ## note: comment out as this entry has been made in db
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='csv', location='~/Data/MetaData/', comment='joana single cell sequencing project annotations for the heavy and light chains for each cell')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# 
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

