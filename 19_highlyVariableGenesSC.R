# File: 19_highlyVariableGenesSC.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: identify highly variable genes in the single cell data and cluster using tsne
# Date: 31/07/2017


## set variables and source libraries
source('header.R')

## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12)')
dfSample = dbGetQuery(db, q)
head(dfSample)

## load the scater object with normalized data
n = paste0(dfSample$location[dfSample$id == 46], dfSample$name[dfSample$id == 46])
load(n)

## load the anti-body annotations assigned from previous analysis
n = paste0(dfSample$location[dfSample$id == 44], dfSample$name[dfSample$id == 44])
dfAntiBody = read.csv(n, header=T, row.names=1, stringsAsFactors = F) 

## load the sample annotations, particularly cell classification using classifier trained via nanoString data
q = paste0('SELECT * FROM Projects.Sample where Sample.idData = 12')
dfSample = dbGetQuery(db, q)
head(dfSample)
# close connection after getting data
dbDisconnect(db)

## drop samples that are not present in the scater object
i = match(pData(oSce.F)$dbID_Sample, dfSample$id)
identical(dfSample$id[i], pData(oSce.F)$dbID_Sample)
dfSample = dfSample[i,]

## drop samples from the antibody table
i = match(pData(oSce.F)$dbID_Sample, dfAntiBody$id)
identical(dfAntiBody$id[i], pData(oSce.F)$dbID_Sample)
dfAntiBody = dfAntiBody[i,]

##################################### single cell sequencing workflow
## see here for a working example https://www.bioconductor.org/help/workflows/simpleSingleCell/
library(scater)
library(scran)
dim(oSce.F)

# Features  Samples 
# 8330       57 

# For each gene, we calculate the percentage of the variance of the
# expression values that is explained by the spike-in totals
plotExplanatoryVariables(oSce.F, variables=c("counts_feature_controls_ERCC",
                                          "log10_counts_feature_controls_ERCC"))
#
## calculate technical variability using all genes
var.fit = trendVar(oSce.F, trend="loess", use.spikes=FALSE, span=0.2)
var.out = decomposeVar(oSce.F, var.fit)
#
## trend fitted to the endogenous variances by examining whether it is consistent with the spike-in variances
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o = order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
cur.spike = isSpike(oSce.F)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
#
# This strategy exploits the large number of endogenous genes to obtain a stable trend,
# with the spike-in transcripts used as diagnostic features rather than in the trend fitting
# itself. However, if our assumption did not hold, we would instead fit the trend directly to the spike-in
# variances with the default use.spikes=TRUE.
var.fit = trendVar(oSce.F, trend="loess", use.spikes=TRUE, span=0.2)
var.out = decomposeVar(oSce.F, var.fit)
#
## trend fitted to the endogenous variances by examining whether it is consistent with the spike-in variances
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o = order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
cur.spike = isSpike(oSce.F)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
#
# HVGs are defined as genes with biological components that are significantly
# greater than zero at a false discovery rate (FDR) of 5%.
# we only consider a gene to be a HVG if it has a biological component greater than or equal to 0.5.
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
#
## [1] 245

write.csv(hvg.out, file="Results/hvg.csv")
head(hvg.out)

plotExpression(oSce.F, rownames(hvg.out)[1:10])

# PCA plot from the normalized log-expression values of the correlated HVGs
plotPCA(oSce.F, exprs_values="exprs", colour_by="total_features", feature_set=rownames(hvg.out))
set.seed(123)
plotTSNE(oSce.F, exprs_values="exprs", colour_by="group1",
         feature_set=rownames(hvg.out))

oSce.F$cellType = dfSample$group3 
set.seed(123)
plotTSNE(oSce.F, exprs_values="exprs", colour_by="cellType",
         feature_set=rownames(hvg.out))

oSce.F$IGH = dfAntiBody$IGH
set.seed(123)
plotTSNE(oSce.F, exprs_values="exprs", colour_by="IGH",
         feature_set=rownames(hvg.out))

oSce.F$IGL = dfAntiBody$IGL
set.seed(123)
plotTSNE(oSce.F, exprs_values="exprs", colour_by="IGL",
         feature_set=rownames(hvg.out))


## use the cdiagnostics library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

mCounts = exprs(oSce.F)
mCounts = mCounts[rownames(hvg.out),]

mCounts = mCounts+min(mCounts)

oDiag = CDiagnosticPlots(mCounts, 'HVG')
fBatch = factor(oSce.F$group1)
fBatch = factor(oSce.F$cellType)
fBatch = factor(oSce.F$IGH)
table(fBatch); nlevels(fBatch)

boxplot.median.summary(oDiag, fBatch)
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
plot.missing.summary(oDiag, fBatch)
plot.PCA(oDiag, fBatch)
plot.dendogram(oDiag, fBatch, labels_cex = 0.7)
plot.heatmap(oDiag)
oDiag = calculateExtremeValues(oDiag)
m1 = mGetExtremeValues(oDiag)
## samples with most extreme values
apply(m1, 2, function(x) sum(x > 0))

## variables that are contributing to this
v1 = apply(m1, 1, function(x) sum(x > 0))
which(v1 > 0)

l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.scaleVariables = F;
l$HC.scaleVaribles = F;

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch, cex.main=1)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l = lapply(l, function(x) x = F)

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l)
plot.PCA(oDiag.2, fBatch, cex.main=1)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)













