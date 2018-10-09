
##################################### ##
## AC and GBS manhattan plots stacked ##
##################################### ##

# combine f.tests.AC.combined and f.tests.GBS.combined and order by chr and snp position
  f.tests.combined <- merge(f.tests.GBS.combined, f.tests.AC.combined, by = c('SNP', 'CHR', 'BP', 'CHRHomoeo'), all = T)
  f.tests.combined <- f.tests.combined[order(f.tests.combined$CHR, f.tests.combined$BP), ]
# replace x and y in colnames with gbs and AC
  colnames(f.tests.combined) <- sub(pattern = 'x$', replacement = 'gbs', colnames(f.tests.combined))
  colnames(f.tests.combined) <- sub(pattern = 'y$', replacement = 'AC', colnames(f.tests.combined))
  renseq.cols <- grep('renseq', colnames(f.tests.combined))
  colnames(f.tests.combined)[renseq.cols] <- paste0(colnames(f.tests.combined)[renseq.cols], '.AC')
# sort chromosomes
  f.tests.combined <- f.tests.combined[, c(1:4, order(colnames(f.tests.combined)[-c(1:4)]) + 4)]
# create a pdf file to hold plots
  pdf('output/manhattan_plots.AC.GBS.pdf', height = 11, width = 8.5)
# create some padding around the plots
  par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
# plot manhattan plot
  for (i in 5:ncol(f.tests.combined)) {
    plotManhattan(dat = f.tests.combined)
  }
  dev.off()  
  
  
## ############################ ##
## Single chrom manhattan plots ##
## ############################ ##
  
  pdf('output/H13.H26.combined.pdf', width = 8.5, height = 7)
  # create some padding around the plots
  par(mfrow=c(2,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  
  plotManhattan_chrom(dat = f.tests.combined, chrom = "6D", gene = "h13", legend.pos = 'topright')
  
  plotManhattan_chrom(dat = f.tests.combined, chrom = "3D", gene = "h26", legend.pos = 'topleft')
  
  dev.off()
  
## ### ##
## END ##
## ### ##
  