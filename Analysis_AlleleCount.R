## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Mapping with Allele counts
## ###########################

# load required functions
  source("functions.R")

# read allele counts from GBS data
  alleleCounts <- fread('data/AlleleCounts.txt', header = T, check.names = F, data.table = F)
  
# read allele counts from RenSeq Data
  alleleCounts_RenSeq_Tv2 <- fread('data/AC_RenSeq-Tv2.txt', header = T, check.names = F, data.table = F)
  
## association tests
  
  # # Carol x Newton - H3
  # associationTest_AC(data = alleleCounts, pop.id = 'CN', gene = 'h3')
  
  # Erin x Newton - H5
  associationTest_AC(data = alleleCounts, pop.id = 'EN', gene = 'h5')
  
  # Flynn x Newton - H6
  associationTest_AC(data = alleleCounts, pop.id = 'FN', gene = 'h6')
  
  # Joy x Newton - H10
  # associationTest_AC(data = alleleCounts, pop.id = 'JN', gene = 'h10')
  
  # Lola x Newton - H12
  associationTest_AC(data = alleleCounts, pop.id = 'LN', gene = 'h12')
  
  # Molly x Newton - H13
  associationTest_AC(data = alleleCounts, pop.id = 'MN', gene = 'h13.newton')
  
  # Molly x Overley - H13
  associationTest_AC(data = alleleCounts, pop.id = 'MO', gene = 'h13.overley')
  
  ## Gene H26 families
      # KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts, pop.id = 'KO', gene = 'h26')
      # RenSeq Tv2 Family 1 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_1', gene = 'h26.renseq.tv2.fam1')
      # RenSeq Tv2 Family 2 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_2', gene = 'h26.renseq.tv2.fam2')
      
  # SynOP DH - H32
  associationTest_AC(data = alleleCounts, pop.id = 'SynOp', gene = 'h32')
  
## ########################### ##
## Manhattan plots in one file ##
## ########################### ##
  
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
  # remove individual dataframes
  rm(list = f.tests.AC)
  # create a pdf file to hold plots
  pdf('output/manhattan_plots.AC.pdf', height = 11, width = 8.5)
  # create some padding around the plots
  par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  
  # plot manhattan plot
  for (i in 5:ncol(f.tests.AC.combined)) {
    plotManhattan(dat = f.tests.AC.combined)
  }
  dev.off()
  
## ### ##
## END ##
## ### ##
  