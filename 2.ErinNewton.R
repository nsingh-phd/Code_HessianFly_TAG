### gene H5

# setting up MAF and missing data threshold
alpha=0.001
maf=0.3
miss=0.2

hapgeno <- h5

res.parent <- 'Erin'
sus.parent <- 'Newton'

hap = cbind(hap.orig[, 1:11], hapgeno)

numlevels=apply(hapgeno, 1, function(x) nlevels(as.factor(x)))
hap=hap[numlevels > 1 & numlevels < 4,]
hapgeno=hapgeno[numlevels > 1 & numlevels < 4,]
cat('# biallelic SNPs =', nrow(hapgeno))

colnames(hap)[6:11] <- c('alleleA', 'alleleB', 'het', 'missing', 'maf', 'present')

a = substring(hap$alleles,1,1)
b = substring(hap$alleles,3,3)
hap$alleleA=rowSums(hapgeno==a, na.rm = T)
hap$alleleB=rowSums(hapgeno==b, na.rm = T)

hap$het=rowSums(hapgeno=='H', na.rm = T)

hap$missing=rowSums(is.na(hapgeno), na.rm = T)/(ncol(hapgeno))
cat('missing data range:', range(hap$missing))

hap$maf= apply(cbind((2*apply(cbind(hap$alleleA, hap$alleleB), 1, min)), hap$het), 1, sum) /
   (2*apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum))

hap$present=rowSums(!is.na(hap[,12:ncol(hap)]), na.rm = T)

hist(hap$missing)
hist(hap$maf)

hap=hap[hap$missing <= miss, ]
hap=hap[hap$maf > maf,]

#    # segregation distortion
#    seg.distort <- rep(NA, nrow(hap))
#       for (i in 1:nrow(hap)) {
#          chi.test <- chisq.test(x = c(hap[i,6], hap[i,8], hap[i,7]), p = c(0.25, 0.5, 0.25))
#          seg.distort[i] = chi.test$p.value
#       }
#
# hap=hap[seg.distort > alpha, ]

rm(hapgeno)

# hap=hap[hap$het <= 0.7, ]
# chiTest <- sapply(apply(hap[,6:8], 1, function(x) chisq.test(x, p=c(.25,.25,.5))$p.value), '[[', 1)
# threshold <- 0.05/length(chiTest)
# hap=hap[chiTest > threshold,]


# position update
for (i in 1:nrow(chrSizeCorrection)) {
   name <- chrSizeCorrection$part_name[i]
   indices <- grep(name, hap$`rs#`, ignore.case = T)
   hap$pos[indices] = hap$pos[indices] + chrSizeCorrection$addThis[i]
}

hap$`rs#`=gsub('(.+_.+)_(\\d+)', '\\1', hap$`rs#`)
hap$`rs#`=paste(hap$`rs#`, hap$pos, sep='_')

geno=as.data.frame(t(hap[,12:ncol(hap)]))
colnames(geno)=hap[,1]
geno=cbind(NA,geno)
geno[,1][grep('resis',rownames(geno), ignore.case = T)]=100
geno[,1][grep(res.parent,rownames(geno), ignore.case = T)]=100
geno[,1][grep('suscep',rownames(geno), ignore.case = T)]=0
geno[,1][grep(sus.parent,rownames(geno), ignore.case = T)]=0

geno=as.matrix(geno)
geno=geno[, geno[rownames(geno)==res.parent,] != geno[rownames(geno)==sus.parent,]]

geno=geno[, !is.na(geno[rownames(geno)==res.parent,])]
geno=geno[, !is.na(geno[rownames(geno)==sus.parent,])]
geno=geno[, geno[rownames(geno)==res.parent,]!='H']
geno=geno[, geno[rownames(geno)==sus.parent,]!='H']

sumNA=rowSums(is.na(geno[,2:ncol(geno)]))
geno=geno[sumNA < (ncol(geno)-1),]

numSNPs=apply(geno, 2, function(x) nlevels(as.factor(x)))
hist(numSNPs)
range(numSNPs)

rownames(geno)
dim(geno)

f.test=cbind(as.character(colnames(geno)[2:ncol(geno)]), NA, NA); colnames(f.test)=c('tag','p.value','p.value.round')
dim(f.test)

for (i in 2:ncol(geno)) {
   fit=aov(geno[,1]~geno[,i])
   temp=summary(fit)
   tryCatch( {
      f.test[i-1,2]=temp[[1]]$`Pr(>F)`[1]
      f.test[i-1,3]=round(temp[[1]]$`Pr(>F)`[1], 5)
   }, error=function(e){}   )
}

f.test <- f.test[complete.cases(f.test), ]

allchr_bwa=cbind(as.data.frame(f.test),NA,NA)
allchr_bwa=allchr_bwa[grep('UN', allchr_bwa[,1], invert = T),]
colnames(allchr_bwa)=c('SNP', 'P', 'CHR', 'BP', 'CHRHomoeo')
allchr_bwa$BP=as.integer(gsub('(.+)_(\\d+)', '\\2', allchr_bwa$SNP))
allchr_bwa$CHR=gsub('(.+)_(\\d+)', '\\1', allchr_bwa$SNP)
allchr_bwa$CHR=substr(allchr_bwa$CHR, 2,3)
allchr_bwa$CHRHomoeo=allchr_bwa$CHR
allchr_bwa$P=as.numeric(as.character(allchr_bwa$P))

for (i in 1:7) {
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'A', sep='')]=i*3-2
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'B', sep='')]=i*3-1
   allchr_bwa$CHR[allchr_bwa$CHR==paste(i,'D', sep='')]=i*3
}

allchr_bwa$CHR=as.integer(allchr_bwa$CHR)
allchr_bwa <- allchr_bwa[allchr_bwa$P > 0, ]; allchr_bwa <- droplevels(allchr_bwa)

chrNames=unique(allchr_bwa$CHRHomoeo)
cols=rep('black', length(chrNames)); cols[grep('B', chrNames)]='red'; cols[grep('D', chrNames)]='blue'

# Bonferroni threshold
cat('Bonferroni level:', alpha/nrow(allchr_bwa))

# range of pvalues
cat('P-value range:', range(allchr_bwa$P))

# significant SNPs
cat('Number of significant SNPs:', sum(allchr_bwa$P < alpha/nrow(allchr_bwa)))

## plot
p_load(qqman)
qq(allchr_bwa$P)

pdf('EN_genomewide.pdf', height = 6.5, width = 11)
manhattan(allchr_bwa, suggestiveline = F, genomewideline = -log10(alpha/ncol(geno)), chrlabs = chrNames, col = cols, ylab='', xlab='')
mtext('Chromosome', 1, line = 2.8, cex = 1.5)
mtext('-log10(p)', 2, line = 2.7, cex = 1.5)

dev.off()

# library(manhattanly)
# manhattanly(allchr_bwa, snp = "SNP", labelChr = chrNames, col = cols, xlab = '', ylab = '', genomewideline = -log10(0.05/nrow(allchr_bwa)), suggestiveline = F)
# mtext('Chromosome', 1, line = 2.8, cex = 2)
# mtext('-log10(p)', 2, line = 2.7, cex = 2)

# subset for one chromosome
chrX_bwa=allchr_bwa[allchr_bwa$CHR==17,] #6B = 17

pdf('EN_6BS.pdf', height = 6.5, width = 11)
manhattan(chrX_bwa, suggestiveline = F, genomewideline = -log10(alpha/nrow(allchr_bwa)), cex.lab=1.4)

abline(v=325245204, col='gray', lwd=7) #6B - Erin
text(x = 325245300, y = 8, labels = 'Centromere', srt=90, pos=4, offset = 0.7, cex = 1)
text(x = 420086000, y = 6.5, labels = 'Bonferroni threshold', pos=4, offset = 0.7, cex = 1)

# calculate possible introgression size
introgression=chrX_bwa[-log10(chrX_bwa$P) > -log10(alpha/nrow(allchr_bwa)), ]

# if(res.parent=='Erin') {
#    introgression=introgression[-nrow(introgression),] #run this only for Erin population
# }

(introgressionSize=range(introgression$BP)[2]-range(introgression$BP)[1])
introgressionRange=range(introgression$BP)/1000000

legend('topright', legend = paste('Introgression size =', round(introgressionSize/1000000, 2), 'Mb', sep = ' '), cex = 1)

dev.off()

# 6B of erin x newton
erin_6b = allchr_bwa

genome.erin.newton = allchr_bwa
#####################