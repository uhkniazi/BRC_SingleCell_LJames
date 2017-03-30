# File: 09_clustering_counts_matrix.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: cluster the samples based on count table for single cell data
# Date: 15/03/2017


## set variables and source libraries
source('header.R')

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
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12) AND (MetaFile.comment like "%count matrix%")')
dfSample = dbGetQuery(db, q)
dfSample
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
           where (Sample.idData = 12) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample.names = dbGetQuery(db, q)
dim(dfSample.names)
dfSample.names
# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location, dfSample$name)

load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))
colnames(mCounts) = dfSample.names$title
## samples with less than 1 million align reads, remove these first
cvSamplesRemove = c('Plate1-C50', 'Plate2-C03', 'Plate2-C26', 'Plate2-C63', 'Plate2-C66')
i = which(colnames(mCounts) %in% cvSamplesRemove)
mCounts = mCounts[,-i]
i = which(dfSample.names$title %in% cvSamplesRemove)
dfSample.names = dfSample.names[-i,]
# sanity check
identical(dfSample.names$title, colnames(mCounts))

##################################### single cell sequencing workflow
## see here for a working example https://www.bioconductor.org/help/workflows/simpleSingleCell/
library(scater)
oSce = newSCESet(countData = mCounts)
dim(oSce)

# Features  Samples 
# 25313       66

## extract spikein and mitochondrial rows
bSpike = grepl('^ERCC', rownames(oSce)); table(bSpike)

## calculate QC metrics 
oSce = calculateQCMetrics(oSce, feature_controls = list(ERCC=bSpike))
head(colnames(pData(oSce)))

## provide spike-in information
library(scran)
isSpike(oSce) = 'ERCC'

### cell QC

# Two common measures of cell quality are the library size and the number of expressed features in each library.
par(mfrow=c(1,2))
hist(oSce$total_counts/1e6, xlab="Library size in millions", main="QC Library Size", prob=T,
     ylab="")
hist(oSce$total_features, xlab="Number of expressed genes", main="QC Features with non-zero counts", 
      ylab="", prob=T)


# remove cells with log-library sizes that are more than 3 MADs 
# below the median log-library sizs
bDropLibsize = isOutlier(oSce$total_counts, nmads=3, type="lower", log=TRUE); table(bDropLibsize)
bDropFeature = isOutlier(oSce$total_features, nmads=3, type="lower", log=TRUE); table(bDropFeature)

# The quantity of spike-in RNA added to each cell should be constant, 
# which means that the proportion should increase upon loss of endogenous RNA in low-quality cells.
hist(oSce$pct_counts_feature_controls_ERCC, xlab="ERCC proportion (%)", 
     ylab="", main="QC Proportion of reads mapped to ERCC")

# cells with more endogenous RNA or that are assayed with protocols using less spike-in RNA will have lower spike-in proportions
bDropSpike = isOutlier(oSce$pct_counts_feature_controls_ERCC, nmads=3, type="higher"); table(bDropSpike)

# Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality.
oSce = oSce[,!(bDropFeature | bDropLibsize | bDropSpike)]
data.frame(ByLibSize=sum(bDropLibsize), ByFeature=sum(bDropFeature),
           BySpike=sum(bDropSpike), Remaining=ncol(oSce))
#             ByLibSize ByFeature BySpike Remaining
# Samples         1         0      9        57

par(p.old)

######### correlation b/w spike-in mix and ercc alignment
mCounts = counts(oSce)
dim(mCounts)
dfERCC = read.csv('~/Data/MetaData/ERCC_concentrations.csv', header=T, sep='\t', stringsAsFactors = F)
i2 = grepl(pattern = '^ERCC', rownames(mCounts))
table(i2)
mCounts = mCounts[i2,]
# order the matrix in the same order as the ERCC names in the dataframe
i = match(dfERCC$ERCC_ID, rownames(mCounts))
mCounts = mCounts[i,]
# sanity check
identical(rownames(mCounts), dfERCC$ERCC_ID)

x = dfERCC$Mix_1
pdf('Temp/mix1.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 1', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())

x = dfERCC$Mix_2
pdf('Temp/mix2.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 2', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())

### choose ERCCs subset that have a reasonable number of reads aligned
m = rowMeans(mCounts)
i = which(m > 10)
mCounts = mCounts[i,]
dfERCC = dfERCC[i,]

x = dfERCC$Mix_1
pdf('Temp/mix1.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 1', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())

x = dfERCC$Mix_2
pdf('Temp/mix2.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 2', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())

# perform a principal components analysis (PCA) based on the quality metrics for each cell
# drop the samples from the main annotation data.frame
i = match(colnames(oSce), dfSample.names$title)
dfSample.names = dfSample.names[i,]
f = dfSample.names$group1
f = as.factor(f)
c = c('black', 'red')
temp = oSce
temp$f = f
plotPCA(temp, pca_data_input="pdata", colour_by='f')

## PCA on non-normalized data 
mCounts = counts(oSce)
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)
plot(pr.out$sdev)

f_Plot3DPCA(pr.out$x, col, pch=20)

# classify cells into cell cycle phases based on the gene expression data
mm.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
library(org.Hs.eg.db)
anno = select(org.Hs.eg.db, keys=rownames(oSce), keytype="ENTREZID", column="ENSEMBL")
ensembl = anno$ENSEMBL[match(rownames(oSce), anno$ENTREZID)]
assignments = cyclone(oSce, mm.pairs, gene.names=ensembl)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16,
     main='Classification of cells in cell cycle phases', col=col)
text(assignments$scores[,c('G1', 'G2M')], labels = dfSample.names$title, pos = 1, cex=0.6)

table(assignments$phases)
# 
# G1 G2M   S 
# 18   4  35 

oSce$group1 = f
oSce$Phase = assignments$phases
oSce$Phase_score = assignments$scores
plotPCA(oSce, pca_data_input="pdata", colour_by='Phase')

## make another pca plot
dfSample.names$Phase = assignments$phases
fSamples = factor(dfSample.names$Phase)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)
f_Plot3DPCA(pr.out$x, col, pch=20)


# low-abundance genes are defined as those with an average count below a filter threshold of 1
ave.counts = rowMeans(counts(oSce))
bKeep = ave.counts >= 1
table(bKeep)
# FALSE  TRUE 
# 16983  8330 

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

# This should generally be dominated by constitutively expressed transcripts, such as those for 
# ribosomal or mitochondrial proteins. The presence of other classes of features may be cause for
# concern if they are not consistent with expected biology. For example, a top set containing many
# spike-in transcripts suggests that too much spike-in RNA was added during library preparation, 
# while the ABSENCE of ribosomal proteins and/or the presence of their pseudogenes are indicative of suboptimal alignment.

plotQC(oSce, type = "highest-expression", n=50)
anno = select(org.Hs.eg.db, keys=rownames(oSce), keytype="ENTREZID", column="SYMBOL")
symb = anno$SYMBOL[match(rownames(oSce), anno$ENTREZID)]
temp = oSce
i = grepl('^ERCC', rownames(temp))
table(i)
symb[i] = paste('ERCC', 1:92)
table(is.na(symb))
i = duplicated(symb)
symb[i] = paste('dup', 1:3)
symb[is.na(symb)] = 'missing'
rownames(temp) = symb
plotQC(temp, type = "highest-expression", n=50)

# An alternative approach to gene filtering is to select genes that have non-zero counts in at least n cells.
# This provides some more protection against genes with outlier expression patterns
numcells <- nexprs(oSce, byrow=TRUE)
bKeepAlt <- numcells >= 10
table(bKeepAlt)

# bKeepAlt
# FALSE  TRUE 
# 22578  2735 
plot(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), pch=20,
                    ylab="Number of expressing cells")

smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), 
              ylab="Number of expressing cells")
is.ercc <- isSpike(oSce, type="ERCC")
points(log10(ave.counts[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)

# apply the mean-based filter to the data 
oSce = oSce[bKeep,]

########### another round of PCA after removing low count genes
## PCA on non-normalized data 
mCounts = counts(oSce)
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)
f_Plot3DPCA(pr.out$x, col=col, pch=19)
########### normalization
# we pool counts from many cells to increase the count size for accurate size factor estimation 
# Pool-based size factors are then “deconvolved” into cell-based factors for cell-specific normalization.
ncol(oSce)
oSce = computeSumFactors(oSce, sizes=c(10, 20, 25))
summary(sizeFactors(oSce))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1221  0.5531  0.8318  1.0000  1.1320  4.4500

plot(sizeFactors(oSce), oSce$total_counts/1e6, log="xy", pch=20,
     ylab="Library size (millions)", xlab="Size factor", col=col, 
     main='Non linear trend between total count and size factor')
text(sizeFactors(oSce), oSce$total_counts/1e6, labels = dfSample.names$title, pos = 1, cex=0.6)

fSamples = factor(dfSample.names$Phase)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(sizeFactors(oSce), oSce$total_counts/1e6, log="xy", pch=20,
     ylab="Library size (millions)", xlab="Size factor", col=col)
text(sizeFactors(oSce), oSce$total_counts/1e6, labels = dfSample.names$title, pos = 1, cex=0.6)

# plot the coverage over the transcripts after loading the appropriate file
##### connect to mysql database to get file name
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12) AND (MetaFile.comment like "%coverage%")')
dfLoadedFile = dbGetQuery(db, q)
dfLoadedFile
# close connection after getting data
dbDisconnect(db)
n = paste0(dfLoadedFile$location, dfLoadedFile$name)
load(n)

dim(mCoverage)
## keep only the samples that are remaining
table(rownames(mCoverage) %in% dfSample.names$title)
mCoverage = mCoverage[rownames(mCoverage) %in% dfSample.names$title,]
dim(mCoverage)

pdf('Temp/coverage.pdf')
par(mfrow=c(3,3))
sapply(1:nrow(mCoverage), function(x) {
  plot(lowess(mCoverage[x,], f=1/10), xaxt='n', xlab='5 to 3 prime binned transcript', ylab='Coverage', type='l', lwd=2,
       main=paste(rownames(mCoverage)[x]), ylim=c(0, 1))
  axis(1, at = c(1, 2000), labels = c('5', '3'), tick = T)})
dev.off(dev.cur())




# we compute a separate set of size factors for the spike-in set. For each cell, 
# the spike-in-specific size factor is defined as the total count across all transcripts in the spike-in set.
# try both ways to see which method is better for normalization i.e. general.use=T and general.use=F
oSce.F = computeSpikeFactors(oSce, type='ERCC', general.use=F)
oSce.T = computeSpikeFactors(oSce, type='ERCC', general.use=T)

# The computed values are stored as an exprs matrix in addition to the other assay elements
# and are log-transformed
oSce.F = normalize(oSce.F)
oSce.T = normalize(oSce.T)

## how is the correlation b/w spike-in and gene based size factors
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pData(oSce.F)$size_factor, pData(oSce.F)$size_factor_ERCC, pch=20, xlab='Gene based size factor', col=col,
     ylab='ERCC size factor', main='Comparisons of two size factors')
text(pData(oSce.F)$size_factor, pData(oSce.F)$size_factor_ERCC, labels = dfSample.names$title, pos = 1, cex=0.6)
lines(lowess(pData(oSce.F)$size_factor, pData(oSce.F)$size_factor_ERCC), lwd=2)

## PCA on normalized data without spike-in 
mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized without spike-in')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)
f_Plot3DPCA(pr.out$x, col, pch=19)
plot(pr.out$sdev)

fSamples = factor(dfSample.names$Phase)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized without spike-in')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## PCA on normalized data WITH spike-in 
mCounts = exprs(oSce.T)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized with spike-in')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)
f_Plot3DPCA(pr.out$x, col, pch=19)

plot(pr.out$sdev)

fSamples = factor(dfSample.names$Phase)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized with spike-in')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

#################################################
### ERCC correlation check after normalization
mCounts = exprs(oSce.T)
dim(mCounts)
dfERCC = read.csv('~/Data/MetaData/ERCC_concentrations.csv', header=T, sep='\t', stringsAsFactors = F)
i2 = grepl(pattern = '^ERCC', rownames(mCounts))
table(i2)
mCounts = mCounts[i2,]
dfERCC = dfERCC[dfERCC$ERCC_ID %in% rownames(mCounts),]
# order the matrix in the same order as the ERCC names in the dataframe
i = match(dfERCC$ERCC_ID, rownames(mCounts))
mCounts = mCounts[i,]
# sanity check
identical(rownames(mCounts), dfERCC$ERCC_ID)

x = dfERCC$Mix_1
pdf('Temp/mix1_norm.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 1', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())

x = dfERCC$Mix_2
pdf('Temp/mix2_norm.pdf')
par(mfrow=c(3, 3))
sapply(1:ncol(mCounts), function(cn){
  cn = colnames(mCounts)[cn]
  plot(log(x), log(mCounts[,cn]+1), pch=20, main=paste(cn), xlab='Log ERCC concentration Mix 2', ylab='Log Aligned')
  lines(lowess(log(x), log(mCounts[,cn]+1)), lwd=2)
})
dev.off(dev.cur())


############## some additional QC checks with distribution of gene expressions
# get the mean vector and total vector for each sample
mCounts.s = exprs(oSce.T)
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$group1)
library(lattice)

densityplot(~ ivMean, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Phenotype',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ condition, data=dfData, auto.key=TRUE, main='Average Gene Expression in Each Sample, SPIKE-IN Normalised',
        xlab='Groups', ylab='Mean Expression', pch=20, cex.axis=0.7)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$phenotype)
dotplot(ivMean ~ condition, data=dfData, auto.key=TRUE, main='Average Gene Expression in Each Sample, SPIKE-IN Normalised',
        xlab='Groups', ylab='Mean Expression', pch=20, scales=list(x=list(cex=0.5, rot=45)))

#########################
i = which(ivMean > 1.4)
dfData[i,]

## repeat with gene normalized data
mCounts.s = exprs(oSce.F)
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$group1)
library(lattice)

densityplot(~ ivMean, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Phenotype',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ condition, data=dfData, auto.key=TRUE, main='Average Gene Expression in Each Sample, Normalised',
        xlab='Groups', ylab='Mean Expression', pch=20, cex.axis=0.7)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$phenotype)
c = c('black', 'red')
c = c[as.numeric(factor(dfSample.names$group1))]
dotplot(ivMean ~ condition, data=dfData, auto.key=TRUE, main='Average Gene Expression in Each Sample, Normalised',
        xlab='Groups', ylab='Mean Expression', pch=20, col=c, scales=list(x=list(cex=0.8, rot=45)))

#########################
i = which(ivMean > 1.4)
dfData[i,]

### samples that share the same capture site
m = as.matrix(xtabs( ~ group1 + phenotype, data=dfSample.names))
m = colSums(m)
i = which(m > 1)
i2 = which(dfSample.names$phenotype %in% names(i))

## PCA on this data but colour by capture sites that are shared
mCounts = exprs(oSce.F)
# n = dim(mCounts)[1] * dim(mCounts)[2]
# mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to FALSE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = F)

# set the factor for colours
## unique colours for matching plates
fLab = as.character(dfSample.names$phenotype)
i = which(as.numeric(table(fLab)) == 2)
n = names(table(fLab))[i]
## these particular names shoud be set to the same level as they are not duplicated
fLab[!fLab %in% n] = 'C00'

g1 = factor(fLab)
iCol = rainbow(nlevels(g1))
iCol[1] = 'lightgrey'
col = iCol[as.numeric(g1)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, Coloured on Shared Capture Sites')
text(pr.out$x[i2,1:2], labels = dfSample.names$title[i2], pos = 1, cex=0.8)




###### association between proportion of reads aligned to ERCC vs Genes
## get count matrix
mCounts = counts(oSce.T)

i2 = grepl(pattern = '^ERCC', rownames(mCounts))
table(i2)
m1 = mCounts[i2,]
m2 = mCounts[!i2,]
m1 = colSums(m1)
m2 = colSums(m2)
mCounts = rbind(m1, m2)
rownames(mCounts) = c('ERCC', 'Genes')
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(t(mCounts), pch=20, main='Reads aligned to genes vs Spike-in', col=col)
text(t(mCounts), labels = colnames(mCounts), pos = 1, cex=0.6)

cs = colSums(mCounts)
mCounts = sweep(mCounts, 2, cs, '/')
i = which(mCounts['ERCC',] < 0.041)

b = barplot((mCounts), xaxt='n', main='Proportion of reads aligned to ERCC (Black) and Genes (Grey)')
axis(1, at = b[-i], labels=colnames(mCounts)[-i], tick = F, las=2, cex.axis=0.6)
axis(1, at = b[i], labels=colnames(mCounts)[i], tick = F, las=2, cex.axis=0.6, col.axis='red')



############################################################################################################
######### Additional steps from tutorial - stop here for first QC
############################################################################################################
# 
# # For each gene, we calculate the percentage of the variance of the 
# # expression values that is explained by the spike-in totals
# plotExplanatoryVariables(oSce, variables=c("counts_feature_controls_ERCC", 
#                                           "log10_counts_feature_controls_ERCC"))
# 
# ## calculate technical variability using all genes
# var.fit = trendVar(oSce, trend="loess", use.spikes=FALSE, span=0.2)
# var.out = decomposeVar(oSce, var.fit)
# 
# ## trend fitted to the endogenous variances by examining whether it is consistent with the spike-in variances
# plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
#      ylab="Variance of log-expression")
# o = order(var.out$mean)
# lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
# cur.spike = isSpike(oSce)
# points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
# 
# # This strategy exploits the large number of endogenous genes to obtain a stable trend, 
# # with the spike-in transcripts used as diagnostic features rather than in the trend fitting 
# # itself. However, if our assumption did not hold, we would instead fit the trend directly to the spike-in
# # variances with the default use.spikes=TRUE. 
# var.fit = trendVar(oSce, trend="loess", use.spikes=TRUE, span=0.2)
# var.out = decomposeVar(oSce, var.fit)
# 
# ## trend fitted to the endogenous variances by examining whether it is consistent with the spike-in variances
# plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
#      ylab="Variance of log-expression")
# o = order(var.out$mean)
# lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
# cur.spike = isSpike(oSce)
# points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
# 
# # HVGs are defined as genes with biological components that are significantly 
# # greater than zero at a false discovery rate (FDR) of 5%.
# # we only consider a gene to be a HVG if it has a biological component greater than or equal to 0.5.
# hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
# hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
# nrow(hvg.out)
# 
# ## [1] 228
# 
# write.csv(hvg.out, file="Results/hvg.csv")
# head(hvg.out)

# plotExpression(oSce, rownames(hvg.out)[1:10])
# 
# # PCA plot from the normalized log-expression values of the correlated HVGs 
# plotPCA(oSce, exprs_values="exprs", colour_by="total_features", feature_set=rownames(oSce))
# plotTSNE(oSce, exprs_values="exprs", perplexity=20, colour_by="total_features", 
#          feature_set=rownames(oSce))






