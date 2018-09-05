## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## ###########################

# load required packages
library(pacman)
p_load(data.table)

# read hapfile
hap.orig <- fread("data/HessianFly_TAG.hmp.txt", header = T, check.names = F, data.table = F)
colnames(hap.orig) # check column names

# check blank wells to see if they are blank or not
# excessive counts of snp call data meaning that they are not blank
colSums(hap.orig[, grep('blank', colnames(hap.orig), ignore.case = T)] != 'N')

# remove unanchored snps and blank sample wells
hap.orig <- hap.orig[(grep('UN', hap.orig$`rs#`, ignore.case = T, invert = T)),
                     (grep('BLANK', colnames(hap.orig), ignore.case = T, invert = T))]
colnames(hap.orig) # check column names

# subset of hapfile to convert IUPAC hets to H
hapgeno.orig <- as.matrix(hap.orig[,12:ncol(hap.orig)])
  # check unique allelic calls
  sort(unique(c(hapgeno.orig)))
  
  hapgeno.orig[hapgeno.orig %in% c('N', '0', '+', '-')] = NA
  hapgeno.orig[hapgeno.orig %in% c('K', 'M', 'R', 'S', 'W', 'Y')] = 'H'
  
  sort(unique(c(hapgeno.orig)))
  hapgeno.orig <- data.frame(hapgeno.orig, stringsAsFactors = F, check.names = F)

# print out hapgeno.orig colnames to see populations
sort(colnames(hapgeno.orig))

# separate out the populations
h3 <- hapgeno.orig[, grep(paste(c('CN_', 'Carol', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h5 <- hapgeno.orig[, grep(paste(c('EN_', 'Erin', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h6 <- hapgeno.orig[, grep(paste(c('FN_', 'Flynn', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h10 <- hapgeno.orig[, grep(paste(c('JN_', 'Joy', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h12 <- hapgeno.orig[, grep(paste(c('LN_', 'Lola', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h13.newton <- hapgeno.orig[, grep(paste(c('MN_', 'Molly', 'Newton'), collapse = '|'), colnames(hapgeno.orig))]
h13.overley <- hapgeno.orig[, grep(paste(c('MO_', 'Molly', 'Overley'), collapse = '|'), colnames(hapgeno.orig))]
h26 <- hapgeno.orig[, grep(paste(c('KO_', 'KU2147', 'Overley'), collapse = '|'), colnames(hapgeno.orig))]
h32 <- hapgeno.orig[, grep(paste(c('SynOp', 'SyntheticW7984', 'OpataM85'), collapse = '|'), colnames(hapgeno.orig))]

#####################