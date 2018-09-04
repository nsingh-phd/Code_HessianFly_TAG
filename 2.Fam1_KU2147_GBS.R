#
alpha=0.001
maf=0.3
miss=0.5
library(data.table)

### start

hap <- fread("hFly_refNRGene_170925.hmp.txt", header = T, check.names = F, data.table = F)


hap <- hap[, -(grep('BLANK', colnames(hap)))]
colnames(hap)[5]='orig_rs'
hap$orig_rs=hap$`rs#`
hap$`rs#` <- sub('S', 'chr', hap$`rs#`, ignore.case=T)

colnames(hap)

hapgeno=hap[,12:ncol(hap)]
hapgeno[hapgeno=='N' | hapgeno=='0' | hapgeno=='+' | hapgeno=='-']='N'
hapgeno[hapgeno=='N']=NA
hapgeno[hapgeno=='R' | hapgeno=='Y' | hapgeno=='S' | hapgeno=='W' | hapgeno=='K' | hapgeno=='M']='H'

hap=cbind(hap[,1:11], hapgeno)
numlevels=apply(hapgeno, 1, function(x) nlevels(as.factor(x)))
hap=hap[numlevels > 1 & numlevels < 4,]
hapgeno=hapgeno[numlevels > 1 & numlevels < 4,]
cat('# biallelic SNPs =', nrow(hapgeno))

colnames(hap)[6:11]=c('alleleA', 'alleleB', 'het', 'missing', 'maf', 'present')
a = substring(hap$alleles,1,1)
b = substring(hap$alleles,3,3)
hap$alleleA=rowSums(hap==a, na.rm = T)
hap$alleleB=rowSums(hap==b, na.rm = T)
hap$het=rowSums(hap=='H', na.rm = T)
hap$missing=rowSums(is.na(hap[,12:ncol(hap)]), na.rm = T)/(ncol(hap)-11)
range(hap$missing)
hap$maf= apply(cbind((2*apply(cbind(hap$alleleA, hap$alleleB), 1, min)), hap$het), 1, sum) /
   (2*apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum))
hap$present=rowSums(!is.na(hap[,12:ncol(hap)]), na.rm = T)
hap$het=rowSums(hap=='H', na.rm = T)/hap$present
hist(hap$missing); hist(hap$maf)
hap=hap[hap$missing <= miss, ]
hap=hap[hap$maf > maf,]
rm(hapgeno)
# hap=hap[hap$het <= 0.7, ]
# chiTest <- sapply(apply(hap[,6:8], 1, function(x) chisq.test(x, p=c(.25,.25,.5))$p.value), '[[', 1)
# threshold <- 0.05/length(chiTest)
# hap=hap[chiTest > threshold,]

for (i in 1:nrow(chrSizeCorrection)) {
   name <- chrSizeCorrection$part_name[i]
   indices <- grep(name, hap$`rs#`, ignore.case = T)
   hap$pos[indices] = hap$pos[indices] + chrSizeCorrection$addThis[i]
}

hap$`rs#`=gsub('(.+_.+)_(\\d+)', '\\1', hap$`rs#`)
hap$`rs#`=paste(hap$`rs#`, hap$pos, sep='_')

geno=as.data.frame(t(hap[,12:ncol(hap)])); colnames(geno)=hap[,1]
geno=cbind(NA,geno)
geno[,1][grep('suscep',rownames(geno))]=0
geno[,1][grep('resis',rownames(geno))]=100
geno[,1][grep('KU21',rownames(geno))]=100
geno[,1][grep('Over',rownames(geno))]=0

geno=as.matrix(geno)
geno=geno[, geno[rownames(geno)=='Overley',]!=geno[rownames(geno)=='KU2147',]]

geno=geno[, !is.na(geno[rownames(geno)=='KU2147',])]
geno=geno[, !is.na(geno[rownames(geno)=='Overley',])]
geno=geno[, geno[rownames(geno)=='Overley',]!='H']
geno=geno[, geno[rownames(geno)=='KU2147',]!='H']

sumNA=rowSums(is.na(geno[,2:ncol(geno)]))
geno=geno[sumNA < (ncol(geno)-1),]

numSNPs=apply(geno, 2, function(x) nlevels(as.factor(x)))
hist(numSNPs)

rownames(geno)
dim(geno)

f.test=cbind(as.character(colnames(geno)[2:ncol(geno)]), NA, NA); colnames(f.test)=c('tag','p.value','p.value.round')
dim(f.test)

for (i in 2:ncol(geno)) {
   fit=aov(geno[,1]~geno[,i])
   temp=summary(fit)
   f.test[i-1,2]=temp[[1]]$`Pr(>F)`[1]
   f.test[i-1,3]=round(temp[[1]]$`Pr(>F)`[1], 5)
}

allchr_bwa=cbind(as.data.frame(f.test),NA,NA)
allchr_bwa=allchr_bwa[-grep('UN', allchr_bwa[,1]),]
colnames(allchr_bwa)=c('SNP', 'P', 'CHR', 'BP', 'CHRHomoeo')
allchr_bwa$BP=as.integer(gsub('(.+_.+)_(\\d+)', '\\2', allchr_bwa$SNP))
allchr_bwa$CHR=gsub('(.+_.+)_(\\d+)', '\\1', allchr_bwa$SNP)
allchr_bwa$CHR=substr(allchr_bwa$CHR, 4,5)
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

### end


# Bonferroni threshold
cat('Bonferroni level:', alpha/nrow(allchr_bwa))

# range of pvalues
cat('P-value range:', range(allchr_bwa$P))

# significant SNPs
cat('Number of significant SNPs:', sum(allchr_bwa$P < alpha/nrow(allchr_bwa)))


library(qqman)
qq(allchr_bwa$P)

pdf('genomewide_KU2147.pdf', height = 4, width = 11)
manhattan(allchr_bwa, suggestiveline = F, genomewideline = -log10(alpha/ncol(geno)),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex.axis = 0.75)
mtext('Chromosome', 1, line = 2.8, cex = 1.5)
mtext('-log10(p)', 2, line = 2.7, cex = 1.5)

dev.off()

# library(manhattanly)
# manhattanly(allchr_bwa, snp = "SNP", labelChr = chrNames, col = cols, xlab = '', ylab = '', genomewideline = -log10(0.05/nrow(allchr_bwa)), suggestiveline = F)
# mtext('Chromosome', 1, line = 2.8, cex = 2)
# mtext('-log10(p)', 2, line = 2.7, cex = 2)

# subset for one chromosome
chr3d_bwa=allchr_bwa[allchr_bwa$CHR==9,]

introgression=chr3d_bwa[-log10(chr3d_bwa$P) > -log10(alpha/nrow(allchr_bwa)), ]
(introgressionSize=range(introgression$BP)[2]-range(introgression$BP)[1])
introgressionRange=range(introgression$BP)/1000000

pdf("manhattanPlot_HF_3DL_bwa.pdf", width = 12, height = 10)
par(mfrow=c(2,1))
plot(chr3d_bwa$BP/1000000, -log10(chr3d_bwa$P), xlab = 'position in Mb', ylab='-log10(p)', pch=20, main = 'Manhattan Plot full chr 3D', xlim = c(0,620), cex.lab=1.4) #, col=col)
legend('topleft', legend = paste('Introgression size =', round(introgressionSize/1000000, 2), 'Mb', sep = ' '), cex = 1)
abline(v=242690774/1000000, col='gray', lwd=7)
abline(h=-log10(0.05/nrow(allchr_bwa)), lty=2, col = 'red')
text(x = 242690774/1000000, y = 15, labels = 'Centromere', srt=90, pos=4, offset = 0.7, cex = 2)
text(x = 50, y = 8, labels = 'Bonferroni threshold', pos=4, offset = 0.7, cex = 1.5)
axis(1, at=c(introgressionRange[1], introgressionRange[2]),labels=F,
     col.axis="red", col.ticks = 'red', las=2, cex.axis=1, tck=.05)
axis(1, at=c(573.87), labels=expression(italic('H26 / H32')),
     col.axis="blue", col.ticks = c('blue'), las=2, cex.axis=0.9, tck=-.05)

plot((chr3d_bwa$BP)/1000000, -log10(chr3d_bwa$P), xlab = 'position in Mb', ylab='-log10(p)', pch=20, main = '20 Mb segment of chr 3D long arm showing introgression', xlim = c(560, 615))
abline(h=-log10(0.001/nrow(allchr_bwa)), lty=2)
text(x = 560, y = 8, labels = 'Bonferroni threshold', pos=4, offset = 0.7, cex = 0.8)
axis(1, at=c(introgressionRange[1], introgressionRange[2]), labels=c('',''),
     col.axis="red", col.ticks = c('red'), cex.axis=0.7, tck=0.03)
text(x = c(introgressionRange[1], introgressionRange[2]), y = 1, labels = c(round(introgressionRange[1],2), round(introgressionRange[2],2)), cex = 0.8, col = 'red')
legend('topleft', legend = paste('Introgression size =', round(introgressionSize/1000000, 2), 'Mb', sep = ' '), cex = 1)
dev.off()

#qplot(chr3d_bwaNew$BP/1000000, -log10(chr3d_bwaNew$P), xlab = 'position in Mb', ylab='-log10(p)', main = 'Manhattan Plot full chr 3D', xlim = c(0,620), ylim = c(0, 70))

###
# bp_local=as.numeric(read.table('bp_local.txt', header = T)$bp); bp_local=bp_local/1000000
# bp_end2end=as.numeric(read.table('bp_end2end.txt', header = T)$bp); bp_end2end=bp_end2end/1000000
# par(mfrow=c(2,1))
# plot(bp_local, rep(1, length(bp_local)), xlim = c(1/1000000, 620), type = 'h')
# abline(v=242.690774, col='gray', lwd=1)
# plot(bp_end2end, rep(1, length(bp_end2end)), xlim = c(1/1000000, 620), type = 'h')
# abline(v=242.690774, col='gray', lwd=1)
# ###

genome.ku2147 = allchr_bwa
