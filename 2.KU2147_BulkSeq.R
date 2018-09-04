## read allele counts
library(data.table)
# alleleCounts_NLR <- fread('AC_RenSeq-Av1.txt', header = T, check.names = F, data.table = F)
alleleCounts_NLR <- fread('AC_RenSeq-Tv2.txt', header = T, check.names = F, data.table = F)
chrSizeCorrection <- read.table('../commonFiles/chrSizeCorrection.bed', header = T, as.is = T)

for (i in 1:nrow(chrSizeCorrection)) {
   indices <- grep(chrSizeCorrection$part_name[i], alleleCounts_NLR$CHR)
   alleleCounts_NLR$POS[indices] = alleleCounts_NLR$POS[indices] + chrSizeCorrection$addThis[i]
}; rm(indices)

## function to compute pValues
assnTest <- function(dat, countCols) {
   dat <- dat[,c(1:4,countCols)] # create a dataset with count and other info for a bulk
   dat <- dat[complete.cases(dat[,5:8]) & rowSums(dat[,5:8] >= 5) == 4, ] # remove SNPs if missing counts & keep SNPs if counts >= 5
   # pVals <- rep(NA, nrow(dat))
   pVals <- matrix(NA, nrow(dat), 3)
   for (i in 1:nrow(dat)) {
      counts = matrix(as.numeric(dat[i,5:8]),2,2)
      if (T
          & sum(counts[,1]) > 50 & sum(counts[,2]) > 50
          & which.max(counts[,1]) != which.max(counts[,2])
          # &  (max(counts[,1])/min(counts[,1])) > 1.25
          # &  (max(counts[,2])/min(counts[,2])) > 1.25
      ) # check if ref in one and alt in the other is higher
      {
         test <- fisher.test(counts)
         pVals[i,1] = test$p.value
         # pVals[i,2] = counts[2,1]/sum(counts[,1])
         # pVals[i,3] = counts[2,2]/sum(counts[,2])
         pVals[i,2] = min(counts[,1])/sum(counts[,1])
         pVals[i,3] = min(counts[,2])/sum(counts[,2])
         # pVals[i,2] = min(counts[,1])/sum(counts)
         # pVals[i,3] = min(counts[,2])/sum(counts)
      }
   }
   pVals=cbind(pVals, 'CHR', 'CHRHomoeo', dat)
   pVals=pVals[complete.cases(pVals),]
   return(pVals)
}

colnames(alleleCounts_NLR)

fam1 <- assnTest(dat = alleleCounts_NLR, countCols = c(6:9))
fam2 <- assnTest(dat = alleleCounts_NLR, countCols = c(10:13))

##
fam=fam1

allchr_bwa=fam
allchr_bwa=allchr_bwa[-grep('UN', allchr_bwa$CHR, ignore.case = T),]
colnames(allchr_bwa)[1:9]=c('P', 'RBlkAltFreq', 'SBlkAltFreq', 'CHR', 'CHRHomoeo', 'SNP', 'BP', 'RefAllele', 'AltAllele')
allchr_bwa$CHR=gsub('(...)(..)_(.+)', '\\2', allchr_bwa$SNP)
allchr_bwa$CHRHomoeo=allchr_bwa$CHR
allchr_bwa$SNP=paste(allchr_bwa$CHRHomoeo, allchr_bwa$BP, sep = '_')
allchr_bwa$P=as.numeric(as.character(allchr_bwa$P))

for (i in 1:7) {
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'A', sep='')]=i*3-2
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'B', sep='')]=i*3-1
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'D', sep='')]=i*3
}

allchr_bwa$CHR=as.integer(allchr_bwa$CHR)

chrNames=unique(allchr_bwa$CHRHomoeo)
cols=rep('black', length(chrNames)); cols[grep('B', chrNames)]='red'; cols[grep('D', chrNames)]='blue'

library(qqman)
manhattan(allchr_bwa, suggestiveline = F, genomewideline = -log10(0.05/nrow(allchr_bwa)), chrlabs = chrNames, col = cols)

genome.bulkseq.fam1 = allchr_bwa[, c(6, 1, 4, 7, 5)]

## fam2
fam=fam2

allchr_bwa=fam
allchr_bwa=allchr_bwa[-grep('UN', allchr_bwa$CHR, ignore.case = T),]
colnames(allchr_bwa)[1:9]=c('P', 'RBlkAltFreq', 'SBlkAltFreq', 'CHR', 'CHRHomoeo', 'SNP', 'BP', 'RefAllele', 'AltAllele')
allchr_bwa$CHR=gsub('(...)(..)_(.+)', '\\2', allchr_bwa$SNP)
allchr_bwa$CHRHomoeo=allchr_bwa$CHR
allchr_bwa$SNP=paste(allchr_bwa$CHRHomoeo, allchr_bwa$BP, sep = '_')
allchr_bwa$P=as.numeric(as.character(allchr_bwa$P))

for (i in 1:7) {
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'A', sep='')]=i*3-2
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'B', sep='')]=i*3-1
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'D', sep='')]=i*3
}

allchr_bwa$CHR=as.integer(allchr_bwa$CHR)

chrNames=unique(allchr_bwa$CHRHomoeo)
cols=rep('black', length(chrNames)); cols[grep('B', chrNames)]='red'; cols[grep('D', chrNames)]='blue'

library(qqman)
manhattan(allchr_bwa, suggestiveline = F, genomewideline = -log10(0.05/nrow(allchr_bwa)), chrlabs = chrNames, col = cols)

genome.bulkseq.fam2 = allchr_bwa[, c(6, 1, 4, 7, 5)]
