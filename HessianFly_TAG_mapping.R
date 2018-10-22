## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## ###########################

## ###########################
## Mapping with GBS SNP calls
## ###########################

# load required functions
  source("HessianFly_TAG_functions.R")
  
# read hapfile
  hap.orig <- fread("data/HessianFly_TAG.hmp.txt", header = T, check.names = F, stringsAsFactors = F, data.table = F)
  # remove 'S' from SNP name and rename the col
    colnames(hap.orig)[1] <- 'rs'
    hap.orig$rs <- sub(pattern = '^S', replacement = '', hap.orig$rs)
# check column names
  colnames(hap.orig)
# utilize extra columns
  colnames(hap.orig)[5:11] <- c('alleleA', 'alleleB', 'het', 'het_prop', 'missing', 'propmiss', 'maf')

# check blank wells to see if they are blank or not
# excessive counts of snp call data meaning that they are not blank
  colSums(hap.orig[, grep('blank', colnames(hap.orig), ignore.case = T)] != 'N')

# remove unanchored snps and blank sample wells
  hap.orig <- hap.orig[(grep('UN', hap.orig$rs, ignore.case = T, invert = T)),
                       (grep('BLANK', colnames(hap.orig), ignore.case = T, invert = T))]
# check column names
  colnames(hap.orig)

# subset of hapfile to convert IUPAC hets to H (convert to matrix for faster computation)
  hapgeno.orig <- as.matrix(hap.orig[,12:ncol(hap.orig)])
  # check unique allelic calls
    sort(unique(c(hapgeno.orig)))
  # replace missing and other ambiguous calls with NA
  # replace all het IUPAC codes to H
    hapgeno.orig[hapgeno.orig %in% c('N', '0', '+', '-')] = NA
    hapgeno.orig[hapgeno.orig %in% c('K', 'M', 'R', 'S', 'W', 'Y')] = 'H'
  # check again if conversion done properly
    sort(unique(c(hapgeno.orig)))
  # convert back to data frame
    hapgeno.orig <- data.frame(hapgeno.orig, stringsAsFactors = F, check.names = F)

# combine hapgeno.orig with original hap.orig 11 cols
    hap.orig <- data.frame(hap.orig[, 1:11], hapgeno.orig, stringsAsFactors = F)

# print out hapgeno.orig colnames to see populations
  sort(colnames(hapgeno.orig))

# separate out the populations
  assign('H5-EN', parse_pop('EN_', 'Erin', 'Newton'))
  assign('H6-FN', parse_pop('FN_', 'Flynn', 'Newton'))
  assign('H12-LN', parse_pop('LN_', 'Lola', 'Newton'))
  assign('H13-MN', parse_pop('MN_', 'Molly', 'Newton'))
  assign('H13-MO', parse_pop('MO_', 'Molly', 'Overley'))
  assign('H26-Fam1', parse_pop('KO_', 'KU2147', 'Overley'))
  assign('H32-SynOpDH', parse_pop('SynOp', 'SyntheticW7984', 'OpataM85'))
  
  # SynOP Phenotypes
  syn.op.pheno <- read.table(file = 'required_files/SynOpDH_hfly_score.txt', header = T, as.is = T)
  
  # remove samples with no phenotype from h32 population
  `H32-SynOpDH` <- `H32-SynOpDH`[, colnames(`H32-SynOpDH`) %in% syn.op.pheno$line]

  # merge and replace h32 sample names with phenotypes
  syn.op.pheno$phenotype[syn.op.pheno$line == 'OpataM85'] = 'OpataM85'
  syn.op.pheno$phenotype[syn.op.pheno$line == 'SyntheticW7984'] = 'SyntheticW7984'
  
  syn.op.pheno <- merge(colnames(`H32-SynOpDH`), syn.op.pheno, by = 1, sort = F)
  all(syn.op.pheno$x == colnames(`H32-SynOpDH`))
  
  colnames(`H32-SynOpDH`) <- syn.op.pheno$phenotype
  colnames(`H32-SynOpDH`)
  
##
## Assiciation tests
##
  
  # Erin x Newton - H5
    associationTest_GBS(dat = `H5-EN`, res.parent = 'Erin', sus.parent = 'Newton', miss = 0.2)
  # Flynn x Newton - H6
    associationTest_GBS(dat = `H6-FN`, res.parent = 'Flynn', sus.parent = 'Newton', miss = 0.2)
  # Lola x Newton - H12
    associationTest_GBS(dat = `H12-LN`, res.parent = 'Lola', sus.parent = 'Newton', miss = 0.2)
  # Molly x Newton - H13
    associationTest_GBS(dat = `H13-MN`, res.parent = 'Molly', sus.parent = 'Newton', miss = 0.2)
  # Molly x Overley - H13
    associationTest_GBS(dat = `H13-MO`, res.parent = 'Molly', sus.parent = 'Overley', miss = 0.2)
  # KU2147 x Overley - H26
    associationTest_GBS(dat = `H26-Fam1`, res.parent = 'KU2147', sus.parent = 'Overley', miss = 0.4)
  # SynOP DH - H32
    associationTest_GBS(dat = `H32-SynOpDH`, res.parent = 'SyntheticW7984', sus.parent = 'OpataM85', miss = 0.2)
  
  # estimate introgression size  
    est.introgression(pop.code = `H5-EN`, chrom = '6B', strategy = 'GBS')
    est.introgression(pop.code = `H6-FN`, chrom = '1A', strategy = 'GBS')
    est.introgression(pop.code = `H12-LN`, chrom = '3B', strategy = 'GBS')
    est.introgression(pop.code = `H13-MN`, chrom = '6D', strategy = 'GBS')
    est.introgression(pop.code = `H13-MO`, chrom = '6D', strategy = 'GBS')
    est.introgression(pop.code = `H26-Fam1`, chrom = '3D', strategy = 'GBS')
    est.introgression(pop.code = `H32-SynOpDH`, chrom = '3D', strategy = 'GBS')

##
## Manhattan plots in one file
##
  
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
  colnames(f.tests.GBS.combined)[-c(1:4)] <- sub(pattern = '$', replacement = '_GBS', colnames(f.tests.GBS.combined)[-c(1:4)])
  # create a pdf file for plotting
  pdf('output/Fig.S3_GWAS_GBS.pdf', height = 11, width = 7)
  # create some padding around the plots
  par(mfrow=c(7,1), oma = c(2, 4, 1, 1), mar = c(2.75, 0, 0, 0)) 
  # plot manhattan plot
  for (i in 5:ncol(f.tests.GBS.combined)) {
    plotManhattan(dat = f.tests.GBS.combined)
  }
  dev.off()
  
##
## Check isolines with Newton
##
  
  # create isolines matrix
    isolines <- c('carol', 'erin', 'flynn', 'iris', 'joy', 'karen', 'lola', 'molly', 'newton')
    isolines.mat <- hap.orig[, grep(paste(isolines, collapse = '|'), colnames(hap.orig), ignore.case = T)]
  # missing data data.frame
    isolines.missing <- colSums(is.na(isolines.mat), na.rm = T)
    isolines.missing <- data.frame('before' = isolines.missing, 'after' = NA, stringsAsFactors = F)
  # replace SNP calls to NA for the chromosomes with introgression
    isolines.mat$Erin[grep('6B', hap.orig$chrom)] = NA
    isolines.mat$Flynn[grep('1A', hap.orig$chrom)] = NA
    isolines.mat$Lola[grep('3B', hap.orig$chrom)] = NA
    isolines.mat$Molly[grep('6D', hap.orig$chrom)] = NA
  # check before and after missing data
    isolines.missing$after <- colSums(is.na(isolines.mat), na.rm = T)
  # allele matching
    id <- alleleMatching(allele.match = isolines.mat)


## ########################## ##
## Mapping with Allele counts ##
## ########################## ##

# read allele counts from GBS data
  alleleCounts <- fread('data/AlleleCounts.txt', header = T, check.names = F, data.table = F)

# read allele counts from RenSeq Data
  alleleCounts_RenSeq_Tv2 <- fread('data/AC_RenSeq-Tv2.txt', header = T, check.names = F, data.table = F)

##
## Assiciation tests
##
  
  # Erin x Newton - H5
    associationTest_BSA(data = alleleCounts, pop.id = 'EN', gene = 'H5-EN')
  # Flynn x Newton - H6
    associationTest_BSA(data = alleleCounts, pop.id = 'FN', gene = 'H6-FN')
  # Lola x Newton - H12
    associationTest_BSA(data = alleleCounts, pop.id = 'LN', gene = 'H12-LN')
  # Molly x Newton - H13
    associationTest_BSA(data = alleleCounts, pop.id = 'MN', gene = 'H13-MN')
  # Molly x Overley - H13
    associationTest_BSA(data = alleleCounts, pop.id = 'MO', gene = 'H13-MO')
  ## Gene H26 families
    # KU2147 x Overley - H26
      associationTest_BSA(data = alleleCounts, pop.id = 'KO', gene = 'H26-Fam1')
      # RenSeq Tv2 Family 1 - KU2147 x Overley - H26
      associationTest_BSA(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_1', gene = 'H26-Fam1_RenSeq')
      # RenSeq Tv2 Family 2 - KU2147 x Overley - H26
      associationTest_BSA(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_2', gene = 'H26-Fam2_RenSeq')
  # SynOP DH - H32
    associationTest_BSA(data = alleleCounts, pop.id = 'SynOp', gene = 'H32-SynOpDH')
    
  # estimate introgression size  
    est.introgression(pop.code = `H5-EN`, chrom = '6B', strategy = 'BSA')
    est.introgression(pop.code = `H6-FN`, chrom = '1A', strategy = 'BSA')
    est.introgression(pop.code = `H12-LN`, chrom = '3B', strategy = 'BSA')
    est.introgression(pop.code = `H13-MN`, chrom = '6D', strategy = 'BSA')
    est.introgression(pop.code = `H13-MO`, chrom = '6D', strategy = 'BSA')
    est.introgression(pop.code = `H26-Fam1`, chrom = '3D', strategy = 'BSA')
    est.introgression(pop.code = `H26-Fam1_RenSeq`, chrom = '3D', strategy = 'BSA')
    est.introgression(pop.code = `H26-Fam2_RenSeq`, chrom = '3D', strategy = 'BSA')
    est.introgression(pop.code = `H32-SynOpDH`, chrom = '3D', strategy = 'BSA')
    
##
## Manhattan plots in one file
##

  # get list of data frames to merge
    f.tests.BSA <- ls(pattern = 'f.test.BSA')
  # change P values colname to gene to aviod conflict
    for (i in 1:length(f.tests.BSA)) {
      dat <- get(f.tests.BSA[i])
      colnames(dat)[2] <- substring(text = f.tests.BSA[i], first = 12)
      assign(paste0(f.tests.BSA[i]), dat, envir = .GlobalEnv)
    }
  # merge all data frames and order based on chromosome and SNP positions
    f.tests.BSA.combined <- Reduce(f = function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP', 'CHR', 'BP', 'CHRHomoeo'), all = T),
                                  x = mget(f.tests.BSA))
    f.tests.BSA.combined <- f.tests.BSA.combined[order(f.tests.BSA.combined$CHR, f.tests.BSA.combined$BP), ]
    colnames(f.tests.BSA.combined)[-c(1:4)] <- sub(pattern = '$', replacement = '_BSA-GBS', 
                                                       colnames(f.tests.BSA.combined)[-c(1:4)])
    colnames(f.tests.BSA.combined) <- sub(pattern = '_RenSeq_BSA-GBS', replacement = '_BSA-RenSeq', colnames(f.tests.BSA.combined))
    f.tests.BSA.combined.no.renseq <- f.tests.BSA.combined[, grep('RenSeq', colnames(f.tests.BSA.combined), invert = T)]
  # create a pdf file to hold plots
    pdf('output/Fig.S4_GWAS_BSA.pdf', height = 11, width = 7)
    # create some padding around the plots
    par(mfrow=c(7,1), oma = c(2, 4, 1, 1), mar = c(2.75, 0, 0, 0)) 
  # plot manhattan plot
    for (i in 5:ncol(f.tests.BSA.combined.no.renseq)) {
      plotManhattan(dat = f.tests.BSA.combined.no.renseq)
    }
    dev.off()


##
## RenSeq manhattan plots
##
    # take renseq samples out in a separate data.frame
    f.tests.BSA.renseq.only <- f.tests.BSA.combined[, c(1:4, grep('RenSeq', colnames(f.tests.BSA.combined)))]
    # update colnames
    colnames(f.tests.BSA.renseq.only) <- sub(pattern = '_BSA-RenSeq', replacement = '', colnames(f.tests.BSA.renseq.only))
    # create a pdf file to hold plots
    pdf('output/Fig.5_GWAS_RenSeq_H26.pdf', height = 6, width = 7)
    # create some padding around the plots
    par(mfrow=c(2,1), oma = c(2, 4, 1, 1), mar = c(2.75, 0, 0, 0)) 
    # plot manhattan plot
    for (i in 5:ncol(f.tests.BSA.renseq.only)) {
      plotManhattan(dat = f.tests.BSA.renseq.only)
    }
    dev.off()
    
##################################### ##
## AC and GBS manhattan plots stacked ##
##################################### ##

# combine f.tests.BSA.combined and f.tests.GBS.combined and order by chr and snp position
  f.tests.combined <- merge(f.tests.GBS.combined, f.tests.BSA.combined, by = c('SNP', 'CHR', 'BP', 'CHRHomoeo'), all = T)
  f.tests.combined <- f.tests.combined[order(f.tests.combined$CHR, f.tests.combined$BP), ]
# sort chromosomes
  f.tests.combined <- f.tests.combined[, c(1:4, order(colnames(f.tests.combined)[-c(1:4)]) + 4)]
# create a pdf file to hold plots
  pdf('output/GWAS_GBS_BSA_Combined.pdf', height = 11, width = 7)
  # create some padding around the plots
  par(mfrow=c(8,1), oma = c(2, 4, 1, 1), mar = c(2.75, 0, 0, 0)) 
  # plot manhattan plot
  for (i in 5:ncol(f.tests.combined)) {
    plotManhattan(dat = f.tests.combined)
  }
  dev.off()  


## ###################################################### ##
## Combine both strategies and plot single manhattan plot ##
## ###################################################### ##
  
  bsa.gbs.only <- f.tests.combined[, grep('renseq', colnames(f.tests.combined), ignore.case = T, invert = T)]
  pop.gene.chrom <- data.frame('pop.code' = c('H5-EN', 'H6-FN', 'H12-LN', 'H13-MN', 'H13-MO', 'H26-Fam1', 'H32-SynOpDH'),
                               'gene' = c('H5', 'H6', 'H12', 'H13', 'H13', 'H26', 'H32'),
                               'chrom' = c('6B', '1A', '3B', '6D', '6D', '3D', '3D'), stringsAsFactors = F)
  # data frame to combine p values from both strategies
  bsa.gbs.combined.p <- data.frame(bsa.gbs.only[, 1:4], 'H5-EN'=NA, 'H6-FN'=NA, 'H12-LN'=NA, 'H13-MN'=NA, 
                                   'H13-MO'=NA, 'H26-Fam1'=NA, 'H32-SynOpDH'=NA, check.names = F, stringsAsFactors = F)
  # pdf file
    pdf('output/Fig.3_GWAS_SigSNPs.pdf', height = 11, width = 7)
  # create some padding around the plots
    par(mfrow=c(7,1), oma = c(2, 4, 1, 1), mar = c(2.75, 0, 0, 0)) 
  # plot manhattan
    for (i in 1:nrow(pop.gene.chrom)) {
    # data frame for specific population columns
      dat <- bsa.gbs.only[, c(1:4, grep(pop.gene.chrom$pop.code[i], colnames(bsa.gbs.only), ignore.case = T))]
    # bonf threshold for each strategy 
      bon.thresh.gbs = alpha / nrow(get(paste0('f.test.GBS.', pop.gene.chrom$pop.code[i])))
      bon.thresh.bsa = alpha / nrow(get(paste0('f.test.BSA.', pop.gene.chrom$pop.code[i])))
    # keep only those snps that are significant in both
      sig.in.both <- (dat[, grep('_BSA-GBS', colnames(dat))] <= bon.thresh.bsa) & 
                     (dat[, grep('_GBS', colnames(dat))] <= bon.thresh.gbs)
    # get the mean of pvals from both strategies
      dat$P = rowMeans(dat[, -c(1:4)])
    # replace non-significant snps to NA
      dat$P[!sig.in.both] = 1
    # append dat to bsa.gbs.combined.p
      bsa.gbs.combined.p[, i+4] <- dat$P
    # set up ylim.max limit
      ylim.max = max(-log10(na.omit(dat$P))) + 3
    # plot manhattan
      manhattan(x = dat, ylim = c(0, ylim.max), suggestiveline = F, genomewideline = F,
                chrlabs = chrNames, col = cols, xlab = '', ylab = '')
      abline(h = -log10(mean(bon.thresh.bsa, bon.thresh.gbs)), lty = 3, col = 'red')
      mtext(text = '-log10(p)', 2, line = 2.7, cex = 1)
      legend('topright', legend = bquote(bold(.(paste0(pop.gene.chrom$pop.code[i], '    ')))), cex = 1, bg = 'gray90')
      if (i == 1) (legend('topleft', legend = paste0(c('A', 'B', 'D'), ' sub-genome'), 
                          pch = c(16), col = c('black', 'red', 'blue'), cex = 1.25))
    }
    dev.off()

##
## Estimate introgressions sizes
##
    
  
## ############################ ##
## Single chrom manhattan plots ##
## ############################ ##

  pdf('output/Fig.4_Combined_H13_H26.pdf', height = 7, width = 8.5)
  # create some padding around the plots
  par(mfrow=c(2,1), oma = c(1, 4, 2, 1), mar = c(4, 0, 0, 0))
  # plot H13
  plotManhattan_chrom(dat = f.tests.combined, chrom = "6D", gene = "H13", legend.pos = 'topright', col.density = 0.4)
  # plot H26 and H32
  plotManhattan_chrom(dat = f.tests.combined, chrom = "3D", gene = c("H26", "H32"), legend.pos = 'topleft', col.density = 0.4)
  dev.off()

## ### ##
## END ##
## ### ##