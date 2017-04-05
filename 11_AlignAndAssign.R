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
oGRloc = GRanges('chr14', IRanges(105707168, 105708665), strand = '*')
oGRloc = GRanges('chr19', IRanges(58345178, 58347029), strand = '+')

f = overlapsAny(oGRLgenes, oGRloc)
table(f)
oGRLgenes[f]
seq = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRLgenes[f])


f = overlapsAny(oGRloc, oGRLgenes)

g1 = import('Temp/gene.fasta')
g2 = import(file.choose())

pwa = pairwiseAlignment(seq[[1]], g1, type='global')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g1), split = ''))
m = rbind(m, c)
dim(m)

se = f_oDNAStringSetConvertPWAMatrix(m)
export(se, 'Temp/pwa.fasta')
writePairwiseAlignments(pwa, file='Temp/pwa_pwa.txt')

pwa = pairwiseAlignment(g2, g1, type='local')
m = as.matrix(pwa)
dim(m)
c = unlist(strsplit(toString(g1), split = ''))
m = rbind(m, c)
dim(m)

se = f_oDNAStringSetConvertPWAMatrix(m)
export(se, 'Temp/pwa.fasta')




pairwiseAlignment(pattern = g1, subject = g2[[2]], type='local')

