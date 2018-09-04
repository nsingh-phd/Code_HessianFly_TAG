
# load data.table for faster data loading
library(data.table)

### start
chrSizeCorrection <- read.table('../commonFiles/chrSizeCorrection.bed', header = T, as.is = T)

# hap.orig <- fread("HF_isolines.hmp.txt", header = T, check.names = F, data.table = F)
hap.orig <- fread("HF_isolines_170925.hmp.txt", header = T, check.names = F, data.table = F)
hap.orig <- hap.orig[(grep('SUN', hap.orig[, 1], ignore.case = T, invert = T)),
                     (grep('BLANK', colnames(hap.orig), ignore.case = T, invert = T))]
colnames(hap.orig)
colnames(hap.orig)[5]='orig_rs'
hap.orig$orig_rs=hap.orig$`rs#`
hap.orig$`rs#` <- sub('S', 'chr', hap.orig$`rs#`, ignore.case=T)

colnames(hap.orig)

hapgeno.orig=hap.orig[,12:ncol(hap.orig)]
hapgeno.orig[hapgeno.orig == 'N' | hapgeno.orig == '0' | hapgeno.orig == '+' | hapgeno.orig == '-'] = NA
hapgeno.orig[hapgeno.orig == 'R' | hapgeno.orig == 'Y' | hapgeno.orig == 'S' | hapgeno.orig == 'W' |
           hapgeno.orig == 'K' | hapgeno.orig == 'M'] = 'H'
colnames(hapgeno.orig)

# separate out the populations
erin <- hapgeno.orig[, grep(paste(c('en', 'erin', 'newton'), collapse = '|'), colnames(hapgeno.orig), ignore.case = T)]
joy <- hapgeno.orig[, grep(paste(c('jn', 'joy', 'newton'), collapse = '|'), colnames(hapgeno.orig), ignore.case = T)]
molly.newton <- hapgeno.orig[, grep(paste(c('mn', 'molly', 'newton'), collapse = '|'), colnames(hapgeno.orig), ignore.case = T)]
molly.overley <- hapgeno.orig[, grep(paste(c('mo', 'molly', 'overley'), collapse = '|'), colnames(hapgeno.orig), ignore.case = T)]

#####################