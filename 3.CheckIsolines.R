# load data.table for faster data loading
library(data.table)

### start
# chrSizeCorrection <- read.table('../commonFiles/chrSizeCorrection.bed', header = T, as.is = T)

hap.orig <- fread("HF_isolines_170925.hmp.txt", header = T, check.names = F, data.table = F)

hap.orig <- hap.orig[(grep('SUN', hap.orig[, 1], ignore.case = T, invert = T)), (grep('BLANK', colnames(hap.orig), ignore.case = T, invert = T))]

hap.orig[hap.orig == 'N' | hap.orig == '0' | hap.orig == '+' | hap.orig == '-'] = 'N'
hap.orig[hap.orig == 'R' | hap.orig == 'Y' | hap.orig == 'S' | hap.orig == 'W' | hap.orig == 'K' | hap.orig == 'M'] = 'H'

# check if isolines or not
isolines <- c('carol', 'erin', 'flynn', 'iris', 'joy', 'karen', 'lola', 'molly', 'newton')
check.isolines = hap.orig[, c(1:11, grep(paste(isolines, collapse = '|'), colnames(hap.orig), ignore.case = T))]

# missing data data.frame
isolines.missing <- colSums(check.isolines[,12:ncol(check.isolines)] == 'N', na.rm = T)
isolines.missing <- data.frame('before'=isolines.missing, 'after'=NA)

# replace SNP calls to N for the chromosomes with introgression
check.isolines$Erin[grep('6B', check.isolines$chrom)] = 'N' # setting chr 6B SNPs to N
check.isolines$Joy[grep('6D', check.isolines$chrom)] = 'N' # setting chr 6D SNPs to N
check.isolines$Molly[grep('6D', check.isolines$chrom)] = 'N' # setting chr 6D SNPs to N

# check before and after missing data
isolines.missing$after <- colSums(check.isolines[,12:ncol(check.isolines)] == 'N', na.rm = T)



# read the data in correct format
source('../../git/gbsPackage/gbs/R/hap.read.R')
source('../../git/gbsPackage/gbs/R/allele.match.R')

hap = hap.read(hap.obj = check.isolines, data.col = 12)
id = allele.match(hap = hap)
id

write.csv(id, file = 'isolines.id.csv', quote = F)

#############