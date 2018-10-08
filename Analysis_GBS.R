## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Mapping with GBS SNP calls
## ###########################

# load required functions
  source("functions.R")
  
# read hapfile
  hap.orig <- fread("data/HessianFly_TAG.hmp.txt", header = T, check.names = F, stringsAsFactors = F, data.table = F)
# check column names
  colnames(hap.orig)
# utilize extra columns
  colnames(hap.orig)[5:11] <- c('alleleA', 'alleleB', 'het', 'het_prop', 'missing', 'propmiss', 'maf')

# check blank wells to see if they are blank or not
# excessive counts of snp call data meaning that they are not blank
  colSums(hap.orig[, grep('blank', colnames(hap.orig), ignore.case = T)] != 'N')

# remove unanchored snps and blank sample wells
  hap.orig <- hap.orig[(grep('UN', hap.orig$`rs#`, ignore.case = T, invert = T)),
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

# function to parse out populations
  parse_pop <- function(id=NULL, r.parent=NULL, s.parent=NULL) {
    datf <- hapgeno.orig[, grep(paste(c(id, r.parent, s.parent), collapse = '|'), colnames(hapgeno.orig))]
    rownames(datf) <- hap.orig$`rs#`
    return(datf)
  }

# separate out the populations
  # h3 <- parse_pop('CN_', 'Carol', 'Newton')
  h5 <- parse_pop('EN_', 'Erin', 'Newton')
  h6 <- parse_pop('FN_', 'Flynn', 'Newton')
  # h10 <- parse_pop('JN_', 'Joy', 'Newton')
  h12 <- parse_pop('LN_', 'Lola', 'Newton')
  h13.newton <- parse_pop('MN_', 'Molly', 'Newton')
  h13.overley <- parse_pop('MO_', 'Molly', 'Overley')
  h26 <- parse_pop('KO_', 'KU2147', 'Overley')
  h32 <- parse_pop('SynOp', 'SyntheticW7984', 'OpataM85')
  
  # SynOP Phenotypes
  syn.op.pheno <- read.table(file = 'required_files/SynOpDH_hfly_score.txt', header = T, as.is = T)
  
  # remove samples with no phenotype from h32 population
  h32 <- h32[, colnames(h32) %in% syn.op.pheno$line]

  # merge and replace h32 sample names with phenotypes
  syn.op.pheno$phenotype[syn.op.pheno$line == 'OpataM85'] = 'OpataM85'
  syn.op.pheno$phenotype[syn.op.pheno$line == 'SyntheticW7984'] = 'SyntheticW7984'
  
  syn.op.pheno <- merge(colnames(h32), syn.op.pheno, by = 1, sort = F)
  all(syn.op.pheno$x == colnames(h32))
  
  colnames(h32) <- syn.op.pheno$phenotype
  colnames(h32)
  
## ##################
## Assiciation tests
## ##################

  # Carol x Newton - H3
  associationTest_GBS(dat = h3, res.parent = 'Carol', sus.parent = 'Newton')
  # est.introgression(dat = , chrom = '')
  
  # Erin x Newton - H5
  associationTest_GBS(dat = h5, res.parent = 'Erin', sus.parent = 'Newton')
  est.introgression(dat = h5, chrom = '6B')
  
  # Flynn x Newton - H6
  associationTest_GBS(dat = h6, res.parent = 'Flynn', sus.parent = 'Newton')
  est.introgression(dat = h6, chrom = '1A')
  
  # Joy x Newton - H10
  associationTest_GBS(dat = h10, res.parent = 'Joy', sus.parent = 'Newton')
  est.introgression(dat = h10, chrom = '6D')
  
  # Lola x Newton - H12
  associationTest_GBS(dat = h12, res.parent = 'Lola', sus.parent = 'Newton')
  est.introgression(dat = h12, chrom = '3B')
  
  # Molly x Newton - H13
  associationTest_GBS(dat = h13.newton, res.parent = 'Molly', sus.parent = 'Newton')
  est.introgression(dat = h13.newton, chrom = '6D')
  
  # Molly x Overley - H13
  associationTest_GBS(dat = h13.overley, res.parent = 'Molly', sus.parent = 'Overley')
  est.introgression(dat = h13.overley, chrom = '6D')
  
  # KU2147 x Overley - H26
  associationTest_GBS(dat = h26, res.parent = 'KU2147', sus.parent = 'Overley', miss = 0.5)
  est.introgression(dat = h26, chrom = '3D')
  
  # SynOP DH - H32
  associationTest_GBS(dat = h32, res.parent = 'SyntheticW7984', sus.parent = 'OpataM85')
  est.introgression(dat = h32, chrom = '3D')
  
## ########################### ##
## Manhattan plots in one file ##
## ########################### ##
  
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
  # remove individual dataframes
  rm(list = f.tests.GBS)
  # create a pdf file for plotting
  pdf('output/manhattan_plots.GBS.pdf', height = 11, width = 8.5)
  # create some padding around the plots
  par(mfrow=c(4,1), oma = c(2, 4, 1, 1), mar = c(3.5, 0, 0, 0)) 
  # plot manhattan plot
  for (i in 5:ncol(f.tests.GBS.combined)) {
    plotManhattan(dat = f.tests.GBS.combined)
  }
  dev.off()
  
## ########################## ##
## Check isolines with Newton ##
## ########################## ##
  
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

## ### ##
## END ##
## ### ##
