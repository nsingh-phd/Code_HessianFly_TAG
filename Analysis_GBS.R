## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## ###########################

# load required functions, files, and packages
  chrom.info <- read.table('data/chrom_info.txt', header = T, as.is = T)
  source("functions.R")
  if(!require(pacman)) install.packages(pacman); require(pacman)
  p_load(data.table, qqman)

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

# print out hapgeno.orig colnames to see populations
  sort(colnames(hapgeno.orig))

# function to parse out populations
  parse_pop <- function(id=NULL, r.parent=NULL, s.parent=NULL) {
    return(hapgeno.orig[, grep(paste(c(id, r.parent, s.parent), collapse = '|'), colnames(hapgeno.orig))])
  }

# separate out the populations
  h3 <- parse_pop('CN_', 'Carol', 'Newton')
  h5 <- parse_pop('EN_', 'Erin', 'Newton')
  h6 <- parse_pop('FN_', 'Flynn', 'Newton')
  h10 <- parse_pop('JN_', 'Joy', 'Newton')
  h12 <- parse_pop('LN_', 'Lola', 'Newton')
  h13.newton <- parse_pop('MN_', 'Molly', 'Newton')
  h13.overley <- parse_pop('MO_', 'Molly', 'Overley')
  h26 <- parse_pop('KO_', 'KU2147', 'Overley')
  h32 <- parse_pop('SynOp', 'SyntheticW7984', 'OpataM85')
  
  # SynOP Phenotypes
  syn.op.pheno <- read.table(file = 'data/SynOpDH_hfly_score.txt', header = T, as.is = T)
  
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
  associationTest(dat = h3, res.parent = 'Carol', sus.parent = 'Newton')
  # est.introgression(dat = , chrom = '')
  
  # Erin x Newton - H5
  associationTest(dat = h5, res.parent = 'Erin', sus.parent = 'Newton')
  est.introgression(dat = h5, chrom = '6B')
  
  # Flynn x Newton - H6
  associationTest(dat = h6, res.parent = 'Flynn', sus.parent = 'Newton')
  est.introgression(dat = h6, chrom = '1A')
  
  # Joy x Newton - H10
  associationTest(dat = h10, res.parent = 'Joy', sus.parent = 'Newton')
  est.introgression(dat = h10, chrom = '6D')
  
  # Lola x Newton - H12
  associationTest(dat = h12, res.parent = 'Lola', sus.parent = 'Newton')
  est.introgression(dat = h12, chrom = '3B')
  
  # Molly x Newton - H13
  associationTest(dat = h13.newton, res.parent = 'Molly', sus.parent = 'Newton')
  est.introgression(dat = h13.newton, chrom = '6D')
  
  # Molly x Overley - H13
  associationTest(dat = h13.overley, res.parent = 'Molly', sus.parent = 'Overley')
  est.introgression(dat = h13.overley, chrom = '6D')
  
  # KU2147 x Overley - H26
  associationTest(dat = h26, res.parent = 'KU2147', sus.parent = 'Overley', miss = 0.5)
  est.introgression(dat = h26, chrom = '3D')
  
  # SynOP DH - H32
  associationTest(dat = h32, res.parent = 'SyntheticW7984', sus.parent = 'OpataM85')
  est.introgression(dat = h32, chrom = '3D')
  
  