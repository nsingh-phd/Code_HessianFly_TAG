# combine all datasets in one

colnames(genome.erin.newton)[2] = 'H5-EN'
colnames(genome.joy.newton)[2] = 'H10-JN'
colnames(genome.molly.newton)[2] = 'H13-MN'
colnames(genome.molly.overley)[2] = 'H13-MO'
colnames(genome.ku2147)[2] = 'ku2147'
colnames(genome.bulkseq.fam1)[2] = 'bulkseqfam1'
colnames(genome.bulkseq.fam2)[2] = 'bulkseqfam2'

all.manhattan.plots <- genome.erin.newton[, c(1,3,4,2,5)]

all.manhattan.plots <- merge(all.manhattan.plots, genome.joy.newton, by=1, all = T)
all.manhattan.plots <- merge(all.manhattan.plots, genome.molly.newton, by=1, all = T)
all.manhattan.plots <- merge(all.manhattan.plots, genome.molly.overley, by=1, all = T)
all.manhattan.plots <- merge(all.manhattan.plots, genome.ku2147, by=1, all = T)
all.manhattan.plots <- merge(all.manhattan.plots, genome.bulkseq.fam1, by=1, all = T)
all.manhattan.plots <- merge(all.manhattan.plots, genome.bulkseq.fam2, by=1, all = T)

unique.values <- function(x) {
   values <- unique(x)
   values <- values[!is.na(values)]
}

## consolidate homoeologous chromosome names in one column
homoeo.cols <- grep('chrhomoeo', colnames(all.manhattan.plots), ignore.case = T)
homoeo.names <- apply(all.manhattan.plots[, homoeo.cols], 1, unique.values)

all.manhattan.plots <- all.manhattan.plots[, -homoeo.cols]

all.manhattan.plots = cbind(all.manhattan.plots, homoeo.names)

## consolidate SNP positions in one column
snp.positions.cols <- grep('bp', colnames(all.manhattan.plots), ignore.case = T)
snp.positions <- apply(all.manhattan.plots[, snp.positions.cols], 1, unique.values)

all.manhattan.plots[, 3] = snp.positions
all.manhattan.plots = all.manhattan.plots[, -snp.positions.cols[-1]]

## consolidate the chromosome numbers in one column
chr.cols <- grep('chr', colnames(all.manhattan.plots), ignore.case = T)
chr.nums <- apply(all.manhattan.plots[, chr.cols], 1, unique.values)

all.manhattan.plots[, 2] = chr.nums
all.manhattan.plots = all.manhattan.plots[, -chr.cols[-1]]

## order based on chromosome and SNP positions
all.manhattan.plots <- all.manhattan.plots[order(all.manhattan.plots$homoeo.names, all.manhattan.plots$BP.x), ]

library(qqman)

alpha = 0.001

# assign chromosome colors
   chrNames=unique(as.character(all.manhattan.plots$homoeo.names))
   cols=rep('black', length(chrNames))
   cols[grep('B', chrNames)]='red'
   cols[grep('D', chrNames)]='blue'

###############
## GBS Plots ##
###############

pdf('Fig.3_iso_Mapping.pdf', height = 11, width = 8.5)

par(mfrow=c(4,1), oma = c(2, 4, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(3.5, 0, 0, 0))

# plot manhattan plot

# erin newton
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'H5-EN', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$`H5-EN`), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(A) H5-EN   ')), cex = 1.25, bg = 'gray90')
legend('top', legend = c('A-genome', 'B-genome', 'D-genome', 'Bonferroni threshold'),
       col = c('black', 'red', 'blue', 'red'), lty = c(NA, NA, NA, 3),
       pch = c(16, 16, 16, NA), cex = 1.5, bg = 'white')

# joy newton
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'H10-JN', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$`H10-JN`), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(B) H10-JN   ')), cex = 1.25, bg = 'gray90')

# molly newton
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'H13-MN', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$`H13-MN`), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(C) H13-MN   ')), cex = 1.25, bg = 'gray90')

# molly overley
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'H13-MO', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$`H13-MO`), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(D) H13-MO   ')), cex = 1.25, bg = 'gray90')
mtext('Chromosome', 1, line = 2.8, cex = 1)

dev.off()

### KU2147 mapping

pdf('Fig.5_ku2147_Mapping.pdf', height = 11, width = 8.5)

par(mfrow=c(4,1), oma = c(2, 4, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(3.5, 0, 0, 0))

# plot manhattan plot

# ku2147 fam1
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'ku2147', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$`ku2147`), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(A) KU2147 Fam1   ')), cex = 1.25, bg = 'gray90')
legend('topright', legend = c('A-genome', 'B-genome', 'D-genome', 'Bonferroni threshold'),
       col = c('black', 'red', 'blue', 'red'), lty = c(NA, NA, NA, 3),
       pch = c(16, 16, 16, NA), cex = 1.5, bg = 'white')

# bulkseq fam1
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'bulkseqfam1', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$bulkseqfam1), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
mtext('-log10(p)', 2, line = 2.7, cex = 1)
legend('topleft', legend = expression(bold('(B)  BulkSeq Fam1   ')), cex = 1.25, bg = 'gray90')


# bulkseq fam2
manhattan(all.manhattan.plots, chr = 'CHR.x', bp = 'BP.x', p = 'bulkseqfam2', suggestiveline = F, genomewideline = F,
          ylim = c(0, range(-log10(all.manhattan.plots$bulkseqfam2), na.rm = T)[2] + 2),
          chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
abline(h = -log10(alpha/nrow(all.manhattan.plots)), lty = 3, col = 'red')
legend('topleft', legend = expression(bold('(C)  BulkSeq Fam2   ')), cex = 1.25, bg = 'gray90')
mtext('Chromosome', 1, line = 2.8, cex = 1)

dev.off()

###################
