## ################### ##
## General R functions ##
## Narinder Singh      ##
## ################### ##

## load packages
  library(data.table)
  library(qqman)
  library(reshape)
  library(tidyverse)

## default alpha
  alpha = 0.001

## chrom info file
  chrom.info <- read.table('required_files/chrom_info.txt', header = T, as.is = T)
  
## assign chromosome colors
  chrNames <- sort(unique(chrom.info$chrom_simple))
  cols <- rep('black', length(chrNames))
  cols[grep('B', chrNames)] = 'red'
  cols[grep('D', chrNames)] = 'blue'

## ################################# ##
## function to parse out populations ##
## ################################# ##
parse_pop <- function(id=NULL, r.parent=NULL, s.parent=NULL) {
  datf <- hapgeno.orig[, grep(paste(c(id, r.parent, s.parent), collapse = '|'), colnames(hapgeno.orig))]
  rownames(datf) <- hap.orig$rs
  return(datf)
}
  
## ################### ##
## Compute basic stats ##
## ################### ##

basic.stats <- function(hap = x) {
  # subsetting alleles
    a = substring(hap$alleles,1,1)   
    b = substring(hap$alleles,3,3)
    sum(a == b)
      # compute number of homozygous individuals
        hap$alleleA = alleleA = rowSums(hap[, 12:ncol(hap)] == a, na.rm = T)
        hap$alleleB = alleleB = rowSums(hap[, 12:ncol(hap)] == b, na.rm = T)
      # updating major minor alleles
        a[alleleA<alleleB] = substring(hap$alleles,3,3)[alleleA<alleleB]
        hap$alleleA[alleleA<alleleB] = alleleB[alleleA<alleleB]
        b[alleleA<alleleB] = substring(hap$alleles,1,1)[alleleA<alleleB]
        hap$alleleB[alleleA<alleleB] = alleleA[alleleA<alleleB]
        if (sum(a == b) == 0)
        if (all(hap$alleleA >= hap$alleleB))  {
            hap$alleles = paste(a, b, sep = '/')
        } else cat ('Allele counts are irrespective of major/minor state. \n')
  # compute per SNP heterozygosity
    hap$het = rowSums(hap[, 12:ncol(hap)] == 'H', na.rm = T) 
    hap$het_prop = (hap$het / (ncol(hap) - 11))
  # compute MAF
    hap$maf = (2*hap$alleleB + hap$het) / (2*(hap$alleleA + hap$alleleB + hap$het))
    cat('MAF computation done. \n')
  # compute missing data
    cat('Computing missing data...')
    hap$missing = rowSums(is.na(hap[, 12:ncol(hap)]), na.rm = T) # compute per SNP missing data
    hap$propmiss = hap$missing/(ncol(hap)-11)
    cat('Missing data range:', range(hap$propmiss))
    cat('Done.')
  return(hap)
}

## #################### ##
## Association Test GBS ##
## #################### ##

associationTest_GBS <- function(dat = NULL, res.parent = NULL, sus.parent = NULL, miss = 0.2, maf = 0.3, alpha = 0.001) {
  # change names to strings
    pop.id <- deparse(substitute(dat))

  # print out unique calls
    cat('unique SNP calls =', sort(unique(c(as.matrix(dat)))))
    cat('\n')
  # get alleles for each SNP
    alleles <- apply(X = dat, MARGIN = 1, function(snp) {
                      unique.calls <- unique(c(as.matrix(snp)))
                      paste0(unique.calls[!unique.calls %in% c("H", NA)], collapse = '/') 
                      })
  # replace sites with no allele with NA
    alleles <- gsub(" ", "", alleles)
    alleles[alleles == ""] <- NA
  # combine genotype data with first 11 columns
    hap <- cbind(hap.orig[, 1:11], dat)
  # replace the alleles col with above-computed alleles
    hap$alleles <- alleles
  # check and retain only biallelic snps
    hap <- hap[nchar(hap$alleles) == 3 & !is.na(hap$alleles), ]
    cat('# biallelic SNPs =', nrow(hap))
    cat('\n')
  # compute basic stats for hap file
    hap <- basic.stats(hap = hap)
    cat('\n')
  # filter SNPs
    hap <- hap[hap$propmiss <= miss, ]
    hap <- hap[hap$maf > maf,]
    cat('\n')
    cat('Biallelic SNPs after miss/maf filtering =', nrow(hap))
    cat('\n')
  # geno data  
    geno <- t(hap[,12:ncol(hap)])
    # set colnames of geno
    colnames(geno) <- hap[,1]
    # add an extra column to hold genotype names
    geno <- cbind(NA,geno)
    # create a phenotypic value for resis and susceptible samples
    res.parent.index = grep(res.parent, rownames(geno), ignore.case = T)
    sus.parent.index = grep(sus.parent, rownames(geno), ignore.case = T)
    geno[,1][grep('resis', rownames(geno), ignore.case = T)] = 100
    geno[,1][grep('suscep', rownames(geno), ignore.case = T)] = 0
    geno[,1][res.parent.index] = 100
    geno[,1][sus.parent.index] = 0
    # remove snps where res and susceptible parents have same allele or missing data
    geno <- geno[, geno[res.parent.index,] != geno[sus.parent.index,]]
    geno <- geno[, !(is.na(geno[res.parent.index,]) | is.na(geno[sus.parent.index,]))]
    geno <- geno[, !(geno[res.parent.index,] == 'H' | geno[sus.parent.index,] == 'H')]
    cat('Biallelic SNPs after parental filtering =', ncol(geno) - 1)
    cat('\n')
    # compute and remove if full row is missing
    sumNA <- rowSums(is.na(geno[, 2:ncol(geno)]))
    geno <- geno[sumNA < (ncol(geno) - 1),]
  
  # run association test and get p-value
    f.test <- data.frame('tag' = colnames(geno)[-1], 'p.value' = NA, stringsAsFactors = F)

    cat('\n')
    cat('Performing association...')
    for (i in 2:ncol(geno)) {
      fit <- aov(geno[, 1] ~ geno[, i])
      temp <- summary(fit)
      tryCatch( {
        f.test[i-1, 2] <- temp[[1]]$`Pr(>F)`[1]
      }, error = function(e){}   )
    }
    
    cat('\n')
  
    # remove incomplete rows
    f.test <- f.test[complete.cases(f.test), ]
    # add extra cols to f.test data frame
    f.test <- cbind(as.data.frame(f.test), NA, NA, NA)
    f.test <- f.test[grep('UN', f.test[, 1], invert = T), ]
    colnames(f.test) <- c('SNP', 'P', 'CHR', 'BP', 'CHRHomoeo')
    f.test$BP <- as.integer(gsub('(.+)_(\\d+)', '\\2', f.test$SNP))
    f.test$CHR <- gsub('(.+)_(\\d+)', '\\1', f.test$SNP)
    f.test$CHRHomoeo <- f.test$CHR
    f.test$P <- as.numeric(as.character(f.test$P))
    
  # change the chrom names to numeric
    for (i in 1:7) {
      f.test$CHR[f.test$CHR == paste(i,'A', sep = '')] = i*3-2
      f.test$CHR[f.test$CHR == paste(i,'B', sep = '')] = i*3-1
      f.test$CHR[f.test$CHR == paste(i,'D', sep = '')] = i*3
    }
  
    f.test$CHR <- as.integer(f.test$CHR)
    f.test <- f.test[f.test$P > 0, ]
    f.test <- droplevels(f.test)
    
    chrNames <- unique(f.test$CHRHomoeo)
    cols <- rep('black', length(chrNames))
    cols[grep('B', chrNames)] = 'red'
    cols[grep('D', chrNames)] = 'blue'
  
  # Bonferroni threshold
  cat('Bonferroni level:', alpha/nrow(f.test))
  cat('\n')
  # range of pvalues
  cat('P-value range:', range(f.test$P))
  cat('\n')
  # significant SNPs
  cat('Number of significant SNPs:', sum(f.test$P < alpha/nrow(f.test)))
  cat('\n')
  
  ## plot
  pdf(file = paste0('output/', pop.id, '.GBS.QQ.Plot.pdf'), height = 6.5, width = 11)
  qq(f.test$P)
  dev.off()
  
  pdf(file = paste0('output/', pop.id, '.GBS.genomewide.pdf'), height = 6.5, width = 11)
  manhattan(f.test, suggestiveline = F, genomewideline = -log10(alpha/nrow(f.test)), 
            chrlabs = chrNames, col = cols, ylab='', xlab='')
  mtext('Chromosome', 1, line = 2.8, cex = 1.5)
  mtext('-log10(p)', 2, line = 2.7, cex = 1.5)
  
  dev.off()
  
  # assign f.test to a new object for output
  assign(paste0('f.test.GBS.', pop.id), f.test, envir = .GlobalEnv)
  
}


## ############################# ##
## Association Test Allele count ##
## ############################# ##

associationTest_BSA <- function(data = NULL, pop.id = NULL, gene = NULL, alpha = 0.001) {
  # create a dataset with population specific columns
    dat <- data[, c(1:4, grep(paste0(pop.id, '_'), colnames(data)))]
    dat <- dat[, c(1:4, order(colnames(dat)[-c(1:4)]) + 4)]
    dat <- dat[grep('UN', dat$CHR, ignore.case = T, invert = T), ]
  # remove SNPs if total missing data and each allele count is less than 1
    dat <- dat[complete.cases(dat) & rowSums(dat[, 5:6] == 0) < 2 & rowSums(dat[, 7:8] == 0) < 2, ]
      cat('# Total SNPs:', nrow(dat)); cat('\n')
    dat <- dat[rowSums(dat[, 5:6]) >= 10 & rowSums(dat[, 7:8]) >= 10, ]
    
  # create a data frame to store p.values
    pVals <- as.matrix(data.frame('SNP' = paste0(dat$CHR, '_', dat$POS), 'P' = NA,
                                  'CHR' = dat$CHR, 'BP' = dat$POS, stringsAsFactors = F))
  # compute p.value
    for (i in 1:nrow(dat)) {
      counts = matrix(data = as.numeric(dat[i, 5:8]), nrow = 2, ncol = 2)
      if (T
          # remove if both alleles counts within bulk are equal
            & counts[1, 1] != counts[2, 1] & counts[1, 2] != counts[2, 2]
          # keep snps with higher count of alt in res and ref in sus
            & which.max(counts[, 1]) == 1 & which.max(counts[, 2]) == 2
          # compute ratio of each allele state (ref and alt) across R and S bulk
          # example alt count in R bulk / alt count in S bulk and same for the ref allele
          # kepp snps with ratio of max / min is greater than or equal to 2
            & (max(counts[1, ]) / min(counts[1, ])) >= 2
            & (max(counts[2, ]) / min(counts[2, ])) >= 2
          ) {
              # perform fisher test and assign pvalue to dataframe
                test <- fisher.test(counts)
                pVals[i, 2] <- test$p.value
          }
      }
      # convert to data frame and remove snps with no pvalue
        pVals <- as.data.frame(pVals, stringsAsFactors = F)
        pVals <- pVals[complete.cases(pVals), ]; cat('# Filtered SNPs:', nrow(pVals)); cat('\n')
      # change chrom from text to numeric
        pVals$CHRHomoeo <- pVals$CHR
          for (i in 1:7) {
            pVals$CHR[pVals$CHR == paste(i, 'A', sep = '')] = i * 3 - 2
            pVals$CHR[pVals$CHR == paste(i, 'B', sep = '')] = i * 3 - 1
            pVals$CHR[pVals$CHR == paste(i, 'D', sep = '')] = i * 3
          }
        pVals$CHR <- as.numeric(pVals$CHR)
      # change other cols to numeric
        pVals$BP <- as.numeric(pVals$BP)
        pVals$P <- as.numeric(pVals$P)
      
      ## plot
        pdf(file = paste0('output/', gene, '.BSA.QQ.Plot.pdf'), height = 6.5, width = 11)
        qq(pVals$P)
        dev.off()
        
        chrNames <- unique(pVals$CHRHomoeo)
        cols <- rep('black', length(chrNames))
        cols[grep('B', chrNames)] = 'red'
        cols[grep('D', chrNames)] = 'blue'
        
        pdf(file = paste0('output/', gene, '.BSA.genomewide.pdf'), height = 6.5, width = 11)
        manhattan(pVals, suggestiveline = F, genomewideline = -log10(alpha / nrow(pVals)), 
                  chrlabs = chrNames, col = cols, ylab='', xlab='')
        mtext('Chromosome', 1, line = 2.8, cex = 1.5)
        mtext('-log10(p)', 2, line = 2.7, cex = 1.5)
        
        dev.off()
        
      # assign pvalue dataframe to global variable
        assign(paste0('f.test.BSA.', gene), pVals, envir = .GlobalEnv)
}

## ###################### ##
## Estimate introgression ##
## ###################### ##

est.introgression <- function(pop.code = NULL, chrom = NULL, strategy = NULL, alpha = 0.001) {
  # assign a particular f.test data to pop.code variable
  pop.id <- deparse(substitute(pop.code))
  assign('pop.code', get(paste0('f.test.', strategy, '.', pop.id)))
  
  # chrom to numeric
  chrom.name <- chrom
  chrom <- chrom.info$chrom_num[chrom.info$chrom_simple == chrom]
  # centromere position
  cent <- chrom.info$centromere[chrom.info$chrom_num == chrom] * 10^6

  # get chrom specific snps
  chr_snps <- pop.code[pop.code$CHR == chrom, ]
  
  pdf(file = paste0('output/', pop.id, '.', strategy, '.introgression.pdf'), height = 6.5, width = 11)
  manhattan(chr_snps, suggestiveline = F, genomewideline = -log10(alpha/nrow(pop.code)), 
            cex.lab = 1, xlab = '', ylab = '', xaxt = 'n')
  axis(side = 1, at = seq(0, 900, 100) * 10^6, labels = seq(0, 900, 100))
  mtext(text = '-log10(p)', side = 2, line = 2.5, cex = 1.25)
  mtext(text = paste0('Chromosome ', chrom.name, ' (Mb)'), side = 1, line = 3, cex = 1.25)
  
  points(x = cent, y = 0, pch = '|', cex = 3, col = 'blue') # plot centromere

  # calculate possible introgression size
  introgression <- chr_snps[-log10(chr_snps$P) > -log10(alpha/nrow(pop.code)), ]
  introgressionSize <- range(introgression$BP)[2] - range(introgression$BP)[1]
  introgressionRange <- range(introgression$BP) / 10^6
  
  legend('topright', legend = paste('Introgression size =', round(introgressionSize/1000000, 2), 'Mb', sep = ' '), cex = 1)
  
  dev.off()
}

## ############## ##
## Plot Manhattan ##
## ############## ##

plotManhattan <- function(dat = NULL, strategy = NULL) {
  manhattan(dat, p = colnames(dat)[i], suggestiveline = F, genomewideline = F,
            ylim = c(0, range(-log10(dat[, i]), na.rm = T)[2] + 2),
            chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1)
  abline(h = -log10(alpha / sum(!is.na(dat[, i]))), lty = 3, col = 'red')
  mtext(text = '-log10(p)', 2, line = 2.7, cex = 1)
  legend('topleft', legend = bquote(bold(.(paste0(colnames(dat)[i], strategy, '    ')))), cex = 1, bg = 'gray90')
}

## ################################ ##
## Plot Manhattan single chromosome ##
## ################################ ##

plotManhattan_chrom <- function(dat = NULL, chrom = NULL, gene = NULL, legend.pos = NULL, col.density = 1, pch = 16) {
  # subset the dataframe with specific population and specific chromosome
    data <- dat[dat$CHRHomoeo == chrom, c(1:4, grep(paste(gene, collapse = '|'), colnames(dat), ignore.case = T))]
  # set up ylim.max limit
    ylim.max = min(apply(data[, -c(1:4)], 2, function(p) max(-log10(na.omit(p)))))
  # plot manhattan
    manhattan(data, p = colnames(data)[5], suggestiveline = F, genomewideline = F, ylim = c(0, ylim.max + 30),
              ylab='', xlab='', cex = 1, cex.axis = 1, xaxt = 'n')
    abline(h = -log10(alpha / nrow(data)), lty = 3, col = 'red')
    axis(side = 1, at = seq(0, 900, 100) * 10^6, labels = seq(0, 900, 100))
    mtext(text = '-log10(p)', side = 2, line = 2.25, cex = 1.25)
    mtext(text = paste0('Chromosome ', chrom, ' (Mb)'), side = 1, line = 2.25, cex = 1.25)

    for (i in 6:ncol(data)) {
      points(data$BP, -log10(data[, i]), col = alpha(i-4, col.density), cex = 0.75, pch = pch)
    }
    
    legend(legend.pos, legend = paste0(colnames(data)[-c(1:4)]), pch = c(16), col = c(1:(ncol(data)-4)), cex = 1.25)
}

## ############### ##
## Identity matrix ##
## ############### ##

alleleMatching <- function(allele.match = NULL) {
  nS <- ncol(allele.match)
  id <- matrix(NA, nrow = nS, ncol = nS)
  for (i in 1:nrow(id)){
    id_ct <- rep(NA, length(i:nrow(id)))
    id_pc <- rep(NA, length(i:nrow(id)))
    for (j in i:ncol(id)){
      line1 <- as.character(allele.match[, i])
      line2 <- as.character(allele.match[, j])
      shared <- !is.na(line1) & !is.na(line2) & line1 != "H" & line2 != "H"
      common <- line1[shared] == line2[shared]
      id_pc[j-i+1] <- sum(common) / sum(shared)
      id_ct[j-i+1] <- sum(shared)
    }
    id[i,i:ncol(id)] <- id_ct
    id[i:ncol(id),i] <- id_pc
  }
  rownames(id) <- colnames(allele.match)
  colnames(id) <- colnames(allele.match)
  write.table(id, file = "output/idFull.mat", quote = F, sep = "\t")
  cat('Done computing. Identity matrix has been saved as a tab-delimited file named idFull.mat in the output folder')
  return(id)
}

## ########### ##
## Pheno plots ##
## ########### ##

# barplots for HF and plant color data
phenoPlots <- function(values = NULL, name = NULL, alpha = 0.05, seg = NULL, trait = NULL) {
  if (seg == '1:2:1') (chiTest <- chisq.test(values, p = c(0.25, 0.5, 0.25)))
  if (seg == '1:1') (chiTest <- chisq.test(values[-2], p = c(0.5, 0.5)))
  midpoints <- barplot(values, plot = F)
  if(trait == 'HF') {cols <- c("forestgreen","orange","tomato")
                     xlabs = c("RR", "Rr", "rr")}
  if(trait == 'plant_color') {cols <- c('forestgreen', 'gray', 'green')
                              xlabs = c("'Overley'-like", 'Het', 'Ae. tauschii like')}
  barplot(values, col = cols, ylab = "", 
          cex.names = 1.5, cex.axis = 1.5, cex.lab = 1.5,
          main = name, ylim = c(0, max(values) + 10),
          cex.main = 2, names.arg = xlabs)
  mtext(text = "Number of lines", side = 2, line = 2.5)
  if (seg == '1:2:1') text(midpoints, values - values/2, labels = values, cex = 1.5)
  if (seg == '1:1') text(midpoints[-2], values[-2] - values[-2]/2, labels = values[-2], cex = 1.5)
  
  p = round(chiTest$p.value, digits = 3)
  if (p > alpha) p2 = paste0('p = ', p)
  if (p < alpha) p2 = "p < 0.05"
  if (p < 0.01) p2 = "p < 0.01"
  if (p < 0.001) p2 = "p < 0.001"
  if (p < 0.0001) p2 = "p < 0.0001"
  
  legend('topleft', legend = p2, cex = 1.25, bty = 'n')
  # legend('topleft', legend = p2, cex = 1.25, bg = 'lightgray', bty = 'o')
}

# function to plot bar plots for agronomic data
trait.plot <- function(data = NULL, trait = NULL, trait.name = NULL, plot = NULL) {
  # get mean and std deviation for the trait
    trait.val = data[, grep(trait, colnames(data), ignore.case = T)]
    mean.trait = tapply(as.numeric(trait.val), data$res_sus, mean, na.rm = T)
    sd.trait = tapply(as.numeric(trait.val), data$res_sus, sd, na.rm = T)
  # barplot
    if (plot == 'barplot') {
      bar.trait <- barplot(mean.trait, ylim = c(0, max(as.numeric(trait.val))), xlim = c(0, 30),
                           names.arg = c('RR', 'rr', "'Overley'"), width = 8, space = 0.1,
                           ylab = trait.name, cex.axis = 1.25, cex.names = 1.25, cex.lab = 1.25)
      arrows(bar.trait, mean.trait + sd.trait, bar.trait, mean.trait - sd.trait,
             col = 2, code = 3, length = 0.1, angle = 90)
    }
  # histogram
    if (plot == 'hist') {
      hist(as.numeric(trait.val), xlab = trait.name, breaks = 15,
           main = '', ylab = '# lines', cex.lab = 1.25)
    }
}
