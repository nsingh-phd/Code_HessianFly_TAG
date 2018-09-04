
library(qqman)

alpha = 0.001

plot.manhattan <- function(allchr_bwa) {
   # assign chromosome colors
   chrNames=unique(allchr_bwa$CHRHomoeo)
   cols=rep('black', length(chrNames))
   cols[grep('B', chrNames)]='red'
   cols[grep('D', chrNames)]='blue'

   # plot manhattan plot
   manhattan(allchr_bwa, suggestiveline = F, genomewideline = F,
             ylim = c(0, -log10(range(allchr_bwa$P))[1] + 2),
             chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 0.7)
   abline(h = -log10(alpha/nrow(allchr_bwa)), lty = 3, col = 'red')
   mtext('-log10(p)', 2, line = 2.7, cex = 0.75)
}

###############
## GBS Plots ##
###############

pdf('All_manhattan_genomewide_gbs.pdf', height = 10, width = 8)

par(mfrow=c(5,1), oma = c(2, 4, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(3.5, 0, 0, 0))

# erin newton
plot.manhattan(genome.erin.newton)
legend('topleft', legend = expression(bold('(A) H5-EN   ')), cex = 0.75, bg = 'gray90')
legend('top', legend = c('A-genome', 'B-genome', 'D-genome', 'Bonferroni threshold'),
       col = c('black', 'red', 'blue', 'red'), lty = c(NA, NA, NA, 3), pch = c(16, 16, 16, NA), cex = 0.75, bg = 'white')

# joy newton
plot.manhattan(genome.joy.newton)
legend('topleft', legend = expression(bold('(B) H10-JN   ')), cex = 0.75, bg = 'gray90')

# molly newton
plot.manhattan(genome.molly.newton)
legend('topleft', legend = expression(bold('(C) H13-MN   ')), cex = 0.75, bg = 'gray90')

# molly overley
plot.manhattan(genome.molly.overley)
legend('topleft', legend = expression(bold('(D) H13-MO   ')), cex = 0.75, bg = 'gray90')

# ku2147 fam1
plot.manhattan(genome.ku2147)
legend('topleft', legend = expression(bold('(E) KU2147 Fam1   ')), cex = 0.75, bg = 'gray90')
mtext('Chromosome', 1, line = 2.8, cex = 0.75)

dev.off()

###################
## BulkSeq Plots ##
###################

pdf('All_manhattan_genomewide_bulkseq.pdf', height = 6, width = 8)

par(mfrow=c(2,1), oma = c(2, 4, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(3.5, 0, 0, 0))

# fam1
plot.manhattan(genome.bulkseq.fam1)
legend('topleft', legend = expression(bold('(A)  Fam1   ')), cex = 0.75, bg = 'gray90')
legend('topright', legend = c('A-genome', 'B-genome', 'D-genome', 'Bonferroni threshold'),
       col = c('black', 'red', 'blue', 'red'), lty = c(NA, NA, NA, 3), pch = c(16, 16, 16, NA), cex = 0.75, bg = 'white')

# fam2
plot.manhattan(genome.bulkseq.fam2)
legend('topleft', legend = expression(bold('(B)  Fam2   ')), cex = 0.75, bg = 'gray90')
mtext('Chromosome', 1, line = 2.8, cex = 0.75)

dev.off()

###################
