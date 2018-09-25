## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## ###########################

# load required packages
library(pacman)
p_load(data.table)

# read hapfile
allele.count <- fread("data/AlleleCount.txt", header = T, check.names = F, data.table = F)
colnames(allele.count) # check column names

# check blank wells to see if they are blank or not
# excessive counts of allele counts meaning that they are not blank
colSums(allele.count[, grep('blank', colnames(allele.count), ignore.case = T)])

# remove unanchored snps and blank sample wells
allele.count <- allele.count[(grep('UN', allele.count$CHR, ignore.case = T, invert = T)),
                     (grep('BLANK', colnames(allele.count), ignore.case = T, invert = T))]
colnames(allele.count) # check column names

# print out allele.count colnames to see populations
sort(colnames(allele.count)[-c(1:5)])

# function to parse out populations
parse_pop <- function(id=NULL, r.parent=NULL, s.parent=NULL) {
  return(allele.count[, grep(paste(c(id, r.parent, s.parent), collapse = '|'), colnames(allele.count))])
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
