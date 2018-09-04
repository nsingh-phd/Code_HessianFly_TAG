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

alpha = 0.001

# assign chromosome colors
chrNames=unique(as.character(all.manhattan.plots$homoeo.names))
cols=rep('black', length(chrNames))
cols[grep('B', chrNames)]='red'
cols[grep('D', chrNames)]='blue'

## plot manhattan plots
library(CMplot)

for (i in 4:10) {
   CMplot(na.omit(all.manhattan.plots[,c(1:3,i)]), plot.type = 'd', col = c('black', 'red', 'blue'), band = 4,
          file = "pdf", cex.axis = 0.65, cex = 0.5, r = 2.5, cir.chr = F, cir.legend.cex = 0.5, amplify = F)
}

