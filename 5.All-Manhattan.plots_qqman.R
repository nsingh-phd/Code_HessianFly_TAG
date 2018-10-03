## ##############
## GBS SNP calls
## ##############

  # get list of data frames to merge
    f.tests.GBS <- ls(pattern = 'f.test.GBS')
  # change P values colname to gene to aviod conflict
    for (i in 1:length(f.tests.GBS)) {
      dat <- get(f.tests.GBS[i])
      colnames(dat)[2] <- substring(text = f.tests.GBS[i], first = 12)
      assign(paste0(f.tests.GBS[i]), dat, envir = .GlobalEnv)
    }
  # merge all data frames and order based on chromosome and SNP positions
    f.tests.GBS.combined <- Reduce(f = function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP', 'CHR', 'BP', 'CHRHomoeo'), all = T),
                                   x = mget(f.tests.GBS))
    f.tests.GBS.combined <- f.tests.GBS.combined[order(f.tests.GBS.combined$CHR, f.tests.GBS.combined$BP), ]

## ##############
## Allele counts
## ##############
    
  # get list of data frames to merge
    f.tests.AC <- ls(pattern = 'f.test.AC')
  # change P values colname to gene to aviod conflict
    for (i in 1:length(f.tests.AC)) {
      dat <- get(f.tests.AC[i])
      colnames(dat)[2] <- substring(text = f.tests.AC[i], first = 11)
      assign(paste0(f.tests.AC[i]), dat, envir = .GlobalEnv)
    }
  # merge all data frames and order based on chromosome and SNP positions
    f.tests.AC.combined <- Reduce(f = function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP', 'CHR', 'BP', 'CHRHomoeo'), all = T),
                                  x = mget(f.tests.AC))
    f.tests.AC.combined <- f.tests.AC.combined[order(f.tests.AC.combined$CHR, f.tests.AC.combined$BP), ]

###############
## GBS Plots ##
###############

  pdf('manhattan_plots.GBS.pdf', height = 11, width = 8.5)
  # create some padding around the plots
    par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  
  # plot manhattan plot
    for (i in 5:ncol(f.tests.GBS.combined)) {
      plotManhattan(dat = f.tests.GBS.combined)
    }
  dev.off()

##############
## AC Plots ##
##############
  
  pdf('manhattan_plots.AC.pdf', height = 11, width = 8.5)
  # create some padding around the plots
  par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  
  # plot manhattan plot
  for (i in 5:ncol(f.tests.AC.combined)) {
    plotManhattan(dat = f.tests.AC.combined)
  }
  dev.off()
  