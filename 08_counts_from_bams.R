# File: 08_counts_from_bams.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: generate count tables for transcripts from bam files
# Date: 15/03/2017


## set variables and source libraries
source('header.R')

## load the transcript db objects
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)
library(rtracklayer)

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

# load the ERCC gtf with spikein information
oGRercc = import('~/Data/MetaData/ERCC92.gtf')
strand(oGRercc) = '*'
# reformat metadata column
f = oGRercc$gene_id
df = DataFrame(exon_id=oGRercc$transcript_id, exon_name=NA)
mcols(oGRercc) = df
oGRLercc = split(oGRercc, f)

oGRLgenes = append(oGRLgenes, oGRLercc)


## create the bamfile list from database
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 12) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dfSample$group1 = gsub(" ", "", dfSample$group1, fixed = T)
# create a new file path variable 
dfSample$fp = paste0(dfSample$name)

#### set working directory to appropriate location with bam files
setwd('Data_external/Bams/')
csFiles = list.files('.', pattern = '*.bam$', recursive = T)
# check if these files match the file names in database
table(dfSample$fp %in% csFiles)

## order the data frame by sample
dfSample = dfSample[order(dfSample$title),]
# split the files by titles to reduce memory usage
lFiles = split(dfSample$fp, dfSample$title)

# for each of these bam file lists do the counting
lCounts = lapply(lFiles, function(bfl){
  ## create a bamfiles list object
  oBamFiles = BamFileList(bfl, index=paste0(bfl, '.bai'))
  return(assays(summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F))$counts)
})


## save the summarized experiment object
setwd(gcswd)
n = make.names(paste('lCounts rd louisa single cell rds'))
n2 = paste0('~/Data/MetaData/', n)
save(lCounts, file=n2)

## comment out after first time execution
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='list of Count matrix from louisa single cell sequencing run with quality 10 duplicates removed and strict settings on trimmomatic')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

names(lCounts)
mCounts = do.call(cbind, lCounts)
i = which(rowMeans(mCounts) > 100)

# select transcripts with coverage
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')
i = i[names(i) %in% names(oGRLgenes)]
oGRLgenes = oGRLgenes[names(i)]

# ### take a sample of exons/transcripts
# oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'tx')
# ## select transcripts with length close to  mean
# oGRLgenes = oGRLgenes[sample(1:length(oGRLgenes), 1000, replace = F)]
# length(oGRLgenes)

getTranscriptCoverage = function(bam){
  ###
  f_bin_vector = function(start, end, bins){
    s = floor(seq(start, end, length.out=bins+1))
    e = s-1
    e[length(e)] = s[length(s)]
    length(s) = length(s)-1
    e = e[2:length(e)]
    return(data.frame(start=s, end=e))
  }# f_bin_vector
  ###
  which = unlist(range(oGRLgenes))
  param = ScanBamParam(flag=scanBamFlag(), what = scanBamWhat(), which=which)
  # read the GAlignments object
  oGA = readGAlignments(bam, param=param)
  # get the coverage
  cov = coverageByTranscript(oGA, oGRLgenes, ignore.strand=FALSE)
  mCoverage = sapply(cov, function(temp) {
  # create a binned vector to create views on this coverage
  bins = f_bin_vector(1, sum(width(temp)), bins=100)
  # create views on these bins
  ivCoverage = viewMeans(Views(temp, bins$start, bins$end))
  ivCoverage = ivCoverage / max(ivCoverage)})
  mCoverage = t(mCoverage)
  mCoverage = colMeans(mCoverage, na.rm = T)
  mCoverage = mCoverage/max(mCoverage)
  return(mCoverage)
}


### coverage for each bam file
mCoverage = sapply(dfSample$fp, getTranscriptCoverage)

