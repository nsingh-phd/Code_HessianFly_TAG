## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## ###########################

## ###########################
## Mapping with GBS SNP calls
## ###########################

# load required functions
  source("functions.R")
  
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
  associationTest_GBS(dat = `H5-EN`, res.parent = 'Erin', sus.parent = 'Newton')
  est.introgression(dat = `H5-EN`, chrom = '6B')
  
  # Flynn x Newton - H6
  associationTest_GBS(dat = `H6-FN`, res.parent = 'Flynn', sus.parent = 'Newton')
  est.introgression(dat = `H6-FN`, chrom = '1A')
  
  # Lola x Newton - H12
  associationTest_GBS(dat = `H12-LN`, res.parent = 'Lola', sus.parent = 'Newton')
  est.introgression(dat = `H12-LN`, chrom = '3B')
  
  # Molly x Newton - H13
  associationTest_GBS(dat = `H13-MN`, res.parent = 'Molly', sus.parent = 'Newton')
  est.introgression(dat = `H13-MN`, chrom = '6D')
  
  # Molly x Overley - H13
  associationTest_GBS(dat = `H13-MO`, res.parent = 'Molly', sus.parent = 'Overley')
  est.introgression(dat = `H13-MO`, chrom = '6D')
  
  # KU2147 x Overley - H26
  associationTest_GBS(dat = `H26-Fam1`, res.parent = 'KU2147', sus.parent = 'Overley', miss = 0.5)
  est.introgression(dat = `H26-Fam1`, chrom = '3D')
  
  # SynOP DH - H32
  associationTest_GBS(dat = `H32-SynOpDH`, res.parent = 'SyntheticW7984', sus.parent = 'OpataM85')
  est.introgression(dat = `H32-SynOpDH`, chrom = '3D')
  
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
  pdf('output/GWAS_GBS.pdf', height = 11, width = 8.5)
  # create some padding around the plots
  par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
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
    colnames(f.tests.BSA.combined)[-c(1:4)] <- sub(pattern = '$', replacement = '_BSA', 
                                                       colnames(f.tests.BSA.combined)[-c(1:4)])
  # create a pdf file to hold plots
    pdf('output/GWAS_BSA.pdf', height = 11, width = 8.5)
    # create some padding around the plots
    par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  
  # plot manhattan plot
    for (i in 5:ncol(f.tests.BSA.combined)) {
      plotManhattan(dat = f.tests.BSA.combined)
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
  pdf('output/GWAS_GBS_BSA_Combined.pdf', height = 11, width = 8.5)
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

  pdf('output/Combined_H13_H26.pdf', width = 8.5, height = 7)
  # create some padding around the plots
  par(mfrow=c(2,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0))
  plotManhattan_chrom(dat = f.tests.combined, chrom = "6D", gene = "H13", legend.pos = 'topright')
  plotManhattan_chrom(dat = f.tests.combined, chrom = "3D", gene = "H26", legend.pos = 'topleft')
  dev.off()

## ### ##
## END ##
## ### ##