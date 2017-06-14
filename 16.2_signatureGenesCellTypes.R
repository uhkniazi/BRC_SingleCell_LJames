
library(scater)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
anno = select(org.Hs.eg.db, keys=rownames(oSce), keytype="ENTREZID", column="SYMBOL")
symb = anno$SYMBOL[match(rownames(oSce), anno$ENTREZID)]

mCounts = exprs(oSce.F)

dfSymbols = AnnotationDbi::select(org.Hs.eg.db, cvTopGenes.oneVar, columns = 'ENTREZID', keytype = 'SYMBOL')
table(dfSymbols$ENTREZID %in% rownames(mCounts.SC))
dfSymbols[(dfSymbols$ENTREZID %in% rownames(mCounts.SC)),]
cvTopGenes.oneVar = dfSymbols[(dfSymbols$ENTREZID %in% rownames(mCounts.SC)),'SYMBOL']





## add the cell type id
dfData$fCellID = lNanoString$metaData$group2
table(dfData$fCellID)

## setup the appropriate grouping
fGroups = rep(NA, length.out=length(dfData$fCellID))
fGroups[dfData$fCellID == 'CD19'] = 1
fGroups[dfData$fCellID != 'CD19'] = 0

dfData$fCellID = factor(fGroups)

lData = list(resp=ifelse(dfData$fCellID == 0, 0, 1),
             mModMatrix=model.matrix(as.formula(f), data=dfData))

# set starting values for optimiser
start = c(rep(0, times=ncol(lData$mModMatrix)))
names(start) = colnames(lData$mModMatrix)
# fit model
fit.lap = laplace(mylogpost, start, lData)
### lets take a sample from this 
## parameters for the multivariate t density
tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
colnames(s) = colnames(lData$mModMatrix)
#fit.lap$sir = s

## averages of posterior from sir sample
post = apply(s, 2, mean)

# calculate AIC
iAIC = (lpd(post, lData) - length(start)) * -2

# calculate E(lpd(theta))
eLPD = mean(sapply(1:nrow(s), function(x) lpd(s[x,], lData)))

# calclate ilppd
ilppd = sum(log(sapply(seq_along(lData$resp), function(x) {
  d = list(resp=lData$resp[x], mModMatrix = lData$mModMatrix[x,])
  lppd(s, d)
})))

## effective numbers of parameters pWAIC1
pWAIC1 = 2 * (ilppd - eLPD)

iWAIC = -2 * (ilppd - pWAIC1)

fit.lap$modelCheck = list('AIC'=iAIC, 'pWAIC1'=pWAIC1, 'WAIC'=iWAIC)




















csFile = file.choose()

# load the bam file as GAlignment object
oGABam = readGAlignments(csFile)

# load the bam file using BamFile function or BamFileList
bf = BamFile(csFile)
# get names of sequences
sn = seqinfo(bf)
gr = as(sn, 'GRanges')
sort(seqlengths(bf))


CBamScaffold.getReadWidthSample = function(obj, size=1000) sample(obj@ivReadWidth, size = size, replace = T)

boots = rep(NA, times=1000)
for (i in 1:1000){
  s1 = sample(labs, 36)
  s2 = sample(labs, 36)
  boots[i] = sum(!(labs %in%(unique(c(s1,s2))))) 
}
hist(boots)

chr14:105,707,118-105,708,665






## import the assembled sequence
g2 = import(file.choose())
pwa = pairwiseAlignment(oSeqIGHc, subject = g2[2], type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g2[2]), split = ''))
m = rbind(m, c)
dim(m)
se = f_oDNAStringSetConvertPWAMatrix(m)
names(se) = c(names(oSeqIGHc), 'predicted')
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



# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))
colnames(mCounts) = dfSample.names$title
bSpike = grepl('^ERCC', rownames(oSce)); table(bSpike)

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




pwa = pairwiseAlignment(g2, g1, type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g1), split = ''))
m = rbind(m, c)
dim(m)

se = f_oDNAStringSetConvertPWAMatrix(m)
export(se, 'Temp/pwa.fasta')




pairwiseAlignment(pattern = g1, subject = g2[[2]], type='local')

