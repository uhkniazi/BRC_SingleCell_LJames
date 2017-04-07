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
df = read.csv('Data_external/IGH_coordinates.csv', header=T, stringsAsFactors = F)
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

df = read.csv('Data_external/IGLambdaConstant_coordinates.csv', header=T, stringsAsFactors = F)
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



## import the assembled sequence
g2 = import(file.choose())
pwa = pairwiseAlignment(oSeqIGLc, subject = seq.as.l, type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(seq.as.l), split = ''))
m = rbind(m, c)
dim(m)
se = f_oDNAStringSetConvertPWAMatrix(m)
names(se) = c(names(oSeqIGLc), 'predicted')
export(se, 'Temp/pwa.fasta')
writePairwiseAlignments(pwa, file='Temp/pwa_pwa.txt')

g1 = import('Temp/gene.fasta')

pwa = pairwiseAlignment(seq[[1]], g2[[2]], type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g2[[2]]), split = ''))
m = rbind(m, c)
dim(m)

se = f_oDNAStringSetConvertPWAMatrix(m)
export(se, 'Temp/pwa.fasta')
writePairwiseAlignments(pwa, file='Temp/pwa_pwa.txt')

louisa.igh = import('Temp/igh_louisa_seqs.fasta')
louisa.igh = reverseComplement(louisa.igh)
names(seq[[1]]) = c('exon1', 'exon2', 'exon3')
louisa.igh = append(louisa.igh, seq[[1]])
pwa = pairwiseAlignment(louisa.igh, g2[[2]], type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g2[[2]]), split = ''))
m = rbind(m, c)
dim(m)
se = f_oDNAStringSetConvertPWAMatrix(m)
names(se) = c(names(louisa.igh), 'predicted')
export(se, 'Temp/pwa_igh_louisa-norc.fasta')












pwa = pairwiseAlignment(g2, g1, type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g1), split = ''))
m = rbind(m, c)
dim(m)

se = f_oDNAStringSetConvertPWAMatrix(m)
export(se, 'Temp/pwa.fasta')




pairwiseAlignment(pattern = g1, subject = g2[[2]], type='local')

