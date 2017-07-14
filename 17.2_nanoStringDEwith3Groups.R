# File: 17.2_nanoStringDEwith3Groups.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: perform DE on the genes doing an anova, using only 3 groups decided based on previous analysis
# Date: 30/06/2017


## set variables and source libraries
source('header.R')

library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did2
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 13) AND (MetaFile.comment like "%NanoString%")')
dfSample = dbGetQuery(db, q)
dfSample
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
           where (Sample.idData = 13) AND (File.idSample = Sample.id)')
dfSample.names = dbGetQuery(db, q)
dim(dfSample.names)
dfSample.names
# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location[2], dfSample$name[2])

load(n)

## make count matrix
names(lNanoString)
lNanoString$desc
mCounts.norm = lNanoString$normalizedData$normalized.data

mCounts.norm = mCounts.norm[mCounts.norm$Code.Class == 'Endogenous', ]
rownames(mCounts.norm) = mCounts.norm$Name
mCounts.norm = as.matrix(mCounts.norm[,-c(1:3)])

mCounts = mCounts.norm[rowMeans(mCounts.norm) > 1,]

dfSample.names = lNanoString$metaData
identical(dfSample.names$title, colnames(mCounts.norm))



## create groups and perform anova after fitting random effects model
## create 3 groups
fG = as.character(lNanoString$metaData$group2)
fG[!(fG %in% c('GC', 'PB'))] = 'Others' 
lNanoString$metaData$fCellID = factor(fG, levels = c('Others', 'GC', 'PB'))
fGroups = lNanoString$metaData$fCellID
fPat = paste0(dfSample.names$group1, ':', gsub('louisajames', '', dfSample.names$BatchCartridge))
fPat = factor(fPat)
table(fPat)

library(lmerTest)
library(parallel)
getAnovaPvalue = function(iIndex){
  x = mCounts[iIndex,]
  df = data.frame(x, fGroups, pat=fPat)
  fit = lmerTest::lmer(x ~ fGroups + (1 | pat), data=df)
  return(anova(fit)$`Pr(>F)`)
}

getAdjustedData = function(iIndex){
  x = mCounts[iIndex,]
  df = data.frame(x, fGroups, pat=fPat)
  fit = lmerTest::lmer(x ~ fGroups + (1 | pat), data=df)
  return(predict(fit, re.form=NA))
}


lPvals = mclapply(1:nrow(mCounts), function(iIndexSub) {
  tryCatch(getAnovaPvalue(iIndexSub), error=function(e) NULL)
})

names(lPvals) = rownames(mCounts)
ivPvals = unlist(lPvals)
ivPvals.adj = p.adjust(ivPvals, 'bonf')
summary(ivPvals.adj); hist(ivPvals.adj)

## get the population level adjusted data
lAdjusted = mclapply(names(ivPvals.adj), getAdjustedData)
names(lAdjusted) = names(ivPvals.adj)
mCounts.adj = do.call(rbind, lAdjusted)

# use the top genes to cluster that samples
fTop = ivPvals.adj < 1e-2
table(fTop)
cvTopGenes = rownames(mCounts.adj[fTop,])
oDiag = CDiagnosticPlots(mCounts.adj[fTop,], 'Top 58 at 1% FDR')
l = CDiagnosticPlotsGetParameters(oDiag)
# don't centre and scale subjects, but do it in variable space only
l$PCA.scaleSubjects=F; l$HC.scaleSubjects = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

boxplot.median.summary(oDiag, fGroups, legend.pos = 'topright')
plot.mean.summary(oDiag, fGroups)
plot.sigma.summary(oDiag, fGroups)
plot.missing.summary(oDiag, fGroups)

plot.PCA(oDiag, fGroups)

plot.dendogram(oDiag, fGroups, labels_cex = 0.8, cex.main=0.8)

oDiag2 = CDiagnosticPlots(mCounts[fTop,], 'Top 58 at 1% FDR')
plot.heatmap(oDiag2, main='Top 58 at 1% FDR')

#### use these genes to cluster the single cell data
####### load single cell data
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 12) AND (MetaFile.comment like "%scater%")')
dfSample = dbGetQuery(db, q)
dfSample
# close connection after getting data
dbDisconnect(db)
library(scater)
library(org.Hs.eg.db)
n = paste0(dfSample$location, dfSample$name)
load(n)
# load the single cell normalised count matrix
mCounts.SC = exprs(oSce.F)

mCounts.norm = mCounts.norm[cvTopGenes,]
### match the gene names in the two data sets
## i.e. symbols in nanostring and enterez id in single cell
dfSymbols = AnnotationDbi::select(org.Hs.eg.db, rownames(mCounts.norm), 
                                  columns = 'ENTREZID', keytype = 'SYMBOL')
dfSymbols = na.omit(dfSymbols)
# how many match b/w the 2 data sets
table(dfSymbols$ENTREZID %in% rownames(mCounts.SC))
# 
# FALSE  TRUE 
# 2      52 
dfSymbols = dfSymbols[(dfSymbols$ENTREZID %in% rownames(mCounts.SC)),]
# match these symbols with nanostring data
i = match(dfSymbols$SYMBOL, rownames(mCounts.norm))
identical(dfSymbols$SYMBOL, rownames(mCounts.norm)[i])
mCounts.norm = mCounts.norm[i,]
# replace these names by enterez ids
rownames(mCounts.norm) = dfSymbols$ENTREZID

## prepare test data i.e. single cell data
mCounts.test = mCounts.SC[rownames(mCounts.norm),]
dim(mCounts.test)

mC = mCounts.test
fGroups = colnames(mC)
fGroups = factor(gsub('(Plate\\d)-.+', '\\1', fGroups))

pr.out = prcomp(t(mC), scale = F)
col.p = rainbow(nlevels(fGroups))
col = col.p[as.numeric(fGroups)]
plot(pr.out$x[,1:2], col=col, pch=20, xlab='Z1', ylab='Z2',
     main=paste('PCA comp 1 and 2', ' SC'))
csLabels = colnames(mC)
text(pr.out$x[,1:2], labels = csLabels, pos = 1, cex=0.6)

## make dendogram
hc = hclust(dist(t(mC)))
library(dendextend)
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = col.p[as.numeric(fGroups)][order.dendrogram(dend)]
labels_cex(dend) = 0.8
# Plotting the new dendrogram
plot(dend, main=paste('Hierarchical clustering of distance matrix for', ' SC'), xlab='', sub='Coloured on Batch')

oDiag.sc = CDiagnosticPlots(mC, 'Single Cell Data')
plot.heatmap(oDiag.sc, main='Single Cell Data')