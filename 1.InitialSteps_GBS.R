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
