# File: 02_fastQC.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the fastq files before trimming
# Date: 14/03/2017


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/master/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select title, group1, group2 from Sample where idData=12;")
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)
#### get the names of the fastq files for first sequencing run
setwd('Data_external/Fastq/')
csFiles = list.files('.', pattern = '*.gz')

# match the files to samples by title
i = sapply(dfSample$title, function(x) grep(x, csFiles))

# index of list with no matching files
i2 = sapply(i, length)
# skip this step if i2 has no 0 values
# dfSample = dfSample[which(i2 != 0),]

# match the files to their respective index in the table
lFilesIndex = lapply(dfSample$title, function(x) grep(x, csFiles))
names(lFilesIndex) = dfSample$title


## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(fls, indir, title){
  wd = getwd()
  setwd(indir)
  ob = CFastqQuality(fls, title)
  setwd(wd)
  cat(paste('done', title, '\n'))
  return(ob)
}

ivFilesIndex = unlist(lFilesIndex)

lOb = lapply(ivFilesIndex, function(x){
  write.qa(csFiles[x], getwd(), csFiles[x])
})

setwd(gcswd)
n = make.names(paste('CFastqQuality joana single cell rds'))
lOb$meta = dfSample
lOb$desc = paste('CFastqQuality object from joana single cell sequencing project', date())
n2 = paste0('~/Data/MetaData/', n)
save(lOb, file=n2)

## note: comment out as this entry has been made in db
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='joana single cell sequencing project FASTQ file quality data')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### create the plots of interest
getwd()
lOb$desc = NULL
lOb$meta = NULL
pdf(file='Results/qa.fastq.pdf')

iReadCount = sapply(lOb, CFastqQuality.getReadCount)
iReadCount = iReadCount/1e+6

barplot(iReadCount, las=2, main='Pre-Trim Read Count', ylab = 'No. of Reads in Millions', cex.names =0.25, col=grey.colors(2))

mQuality = sapply(lOb, function(x){
  m = mGetReadQualityByCycle(x)
  m = colMeans(m, na.rm = T)
  return(m)
})

matplot(mQuality, type='l', main='Pre-trim base quality', ylab = 'Mean Score', xlab='Position in Read')

lReadWidth = lapply(lOb, iGetReadWidth)
boxplot(lReadWidth, las=2, main='Pre-trim Read Width', ylab = 'Read Width', col=grey.colors(2), outline=F, xaxt='n')
axis(1, at=1:length(lReadWidth), labels = names(lReadWidth), cex.axis=0.25, las=2)

## some samples may have quality drops which can be seen by clustering
hc = hclust(dist(t(mQuality)))
plot(hc, main='Clustering of Pre-trim data based on Per base quality', cex=0.25, xlab='', sub='')
abline(h = 15, col=2)

## colour by label groups
g3 = as.matrix(table(dfSample$group1, dfSample$group2))
g3 = colSums(g3)
i = match(dfSample$group2, names(g3))
dfSample$group3 = g3[i]
fLab = hc$labels
fLab = sapply(fLab, function(x) substr(x, 1, nchar(x)-1))
i = match(fLab, dfSample$title)
df = data.frame(fLab, title=dfSample$title[i], g1=dfSample$group1[i], g2=dfSample$group2[i], g3=dfSample$group3[i])
## choose a colour by using the g1 or g2 as a factor
g1 = factor(df$g1)
iCol = c('black', 'red'); #rainbow(nlevels(g1))
## plotting a hc object with colours is not straightforward so following example from 
## http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
library(dendextend)
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = iCol[as.numeric(g1)][order.dendrogram(dend)]
labels_cex(dend) = 0.25
# Plotting the new dendrogram
plot(dend, main='Clustering of Pre-trim data based on Per base quality', xlab='', sub='Coloured on Patients')

g1 = factor(df$g3)
iCol = c('black', 'red'); #rainbow(nlevels(g1))
## plotting a hc object with colours is not straightforward so following example from 
## http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
library(dendextend)
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = iCol[as.numeric(g1)][order.dendrogram(dend)]
labels_cex(dend) = 0.25
# Plotting the new dendrogram
plot(dend, main='Clustering of Pre-trim data based on Per base quality', xlab='', sub='Duplicated capture sites')

## unique colours for matching plates
fLab = as.character(df$fLab)
fLab = sapply(fLab, function(x) substr(x, 8, nchar(x)))
i = which(as.numeric(table(fLab)) == 2)
n = names(table(fLab))[i]
## these particular names shoud be set to the same level as they are not duplicated
fLab[fLab %in% n] = 'C00'

g1 = factor(fLab)
iCol = rainbow(nlevels(g1))
iCol[1] = 'black'
## plotting a hc object with colours is not straightforward so following example from 
## http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = iCol[as.numeric(g1)][order.dendrogram(dend)]
labels_cex(dend) = 0.25
# Plotting the new dendrogram
plot(dend, main='Clustering of Pre-trim data based on Per base quality', xlab='', sub='Duplicated capture sites')

####################
# ## A sub tree - so we can see better what we got:
# par(cex = 1)
# plot(dend[[2]], horiz = TRUE)

c = cutree(hc, h=15)
table(c)
# draw for each of the clusters
sapply(unique(c), function(x){
  m = mQuality[,c == x]
  matplot(m, type='l', main='Pre-trim Per base quality Clusters', ylab = 'Mean Score', xlab='Position in Read', col=1:ncol(m))
  legend('bottomleft', legend = c(colnames(m)), fill = 1:ncol(m), ncol=4, cex=0.7)
})
dev.off(dev.cur())




