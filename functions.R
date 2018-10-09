## ################### ##
## General R functions ##
## Narinder Singh      ##
## ################### ##

## load packages
  if(!require(pacman)) install.packages(pacman); require(pacman)
  p_load(data.table, qqman, CMplot)

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

associationTest_AC <- function(data = NULL, pop.id = NULL, gene = NULL, alpha = 0.001) {
  # create a dataset with population specific columns
    dat <- data[, c(1:4, grep(paste0(pop.id, '_'), colnames(data)))]
    dat <- dat[, c(1:4, order(colnames(dat)[-c(1:4)]) + 4)]
    dat <- dat[grep('UN', dat$CHR, ignore.case = T, invert = T), ]
  # remove SNPs if total missing data and each allele count is less than 3
    dat <- dat[complete.cases(dat) & rowSums(dat[, 5:8] >= 1) == 4, ]
  # create a data frame to store p.values
    pVals <- as.matrix(data.frame('SNP' = paste0(dat$CHR, '_', dat$POS), 'P' = NA,
                                  'CHR' = dat$CHR, 'BP' = dat$POS, stringsAsFactors = F))
  # compute p.value
    for (i in 1:nrow(dat)) {
      counts = matrix(data = as.numeric(dat[i, 5:8]), nrow = 2, ncol = 2)
      if (T
          # keep snps with higher count of alt in res and ref in sus
            & which.max(counts[, 1]) == 1 & which.max(counts[, 2]) == 2
          # kepp snps with ratio of max / min is greater than 1.25
            & (max(counts[, 1]) / min(counts[, 1])) >= 1.25
            & (max(counts[, 2]) / min(counts[, 2])) >= 1.25) {
              # perform fisher test and assign pvalue to dataframe
                test <- fisher.test(counts)
                pVals[i, 2] <- test$p.value
          }
      }
      # convert to data frame and remove snps with no pvalue
        pVals <- as.data.frame(pVals, stringsAsFactors = F)
        pVals <- pVals[complete.cases(pVals), ]
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
        pdf(file = paste0('output/', gene, '.AC.QQ.Plot.pdf'), height = 6.5, width = 11)
        qq(pVals$P)
        dev.off()
        
        chrNames <- unique(pVals$CHRHomoeo)
        cols <- rep('black', length(chrNames))
        cols[grep('B', chrNames)] = 'red'
        cols[grep('D', chrNames)] = 'blue'
        
        pdf(file = paste0('output/', gene, '.AC.genomewide.pdf'), height = 6.5, width = 11)
        manhattan(pVals, suggestiveline = F, genomewideline = -log10(alpha / nrow(pVals)), 
                  chrlabs = chrNames, col = cols, ylab='', xlab='')
        mtext('Chromosome', 1, line = 2.8, cex = 1.5)
        mtext('-log10(p)', 2, line = 2.7, cex = 1.5)
        
        dev.off()
        
      # assign pvalue dataframe to global variable
        assign(paste0('f.test.AC.', gene), pVals, envir = .GlobalEnv)
}

## ###################### ##
## Estimate introgression ##
## ###################### ##

est.introgression <- function(dat = NULL, chrom = NULL, alpha = 0.001) {
  # assign a particular f.test data to dat variable
  pop.id <- deparse(substitute(dat))
  assign('dat', get(paste0('f.test.GBS.', pop.id)))
  
  # chrom to numeric
  chrom.name <- chrom
  chrom <- chrom.info$chrom_num[chrom.info$chrom_simple == chrom]
  # centromere position
  cent <- chrom.info$centromere[chrom.info$chrom_num == chrom] * 10^6

  # get chrom specific snps
  chr_snps <- dat[dat$CHR == chrom, ]
  
  pdf(file = paste0('output/', pop.id, '.GBS.introgression.pdf'), height = 6.5, width = 11)
  manhattan(chr_snps, suggestiveline = F, genomewideline = -log10(alpha/nrow(dat)), 
            cex.lab = 1, xlab = '', ylab = '', xaxt = 'n')
  axis(side = 1, at = seq(0, 900, 100) * 10^6, labels = seq(0, 900, 100))
  mtext(text = '-log10(p)', side = 2, line = 2.5, cex = 1.25)
  mtext(text = paste0('Chromosome ', chrom.name, ' (Mb)'), side = 1, line = 3, cex = 1.25)
  
  points(x = cent, y = 0, pch = '|', cex = 3, col = 'blue') # plot centromere

  # calculate possible introgression size
  introgression <- chr_snps[-log10(chr_snps$P) > -log10(alpha/nrow(dat)), ]
  introgressionSize <- range(introgression$BP)[2] - range(introgression$BP)[1]
  introgressionRange <- range(introgression$BP) / 10^6
  
  legend('topright', legend = paste('Introgression size =', round(introgressionSize/1000000, 2), 'Mb', sep = ' '), cex = 1)
  
  dev.off()
}

## ############## ##
## Plot Manhattan ##
## ############## ##

plotManhattan <- function(dat = NULL) {
  manhattan(dat, p = colnames(dat)[i], suggestiveline = F, genomewideline = F,
            ylim = c(0, range(-log10(dat[, i]), na.rm = T)[2] + 2),
            chrlabs = chrNames, col = cols, ylab='', xlab='', cex = 1, cex.axis = 1.25)
  abline(h = -log10(alpha / sum(!is.na(dat[, i]))), lty = 3, col = 'red')
  mtext(text = '-log10(p)', 2, line = 2.7, cex = 1)
  legend('topleft', legend = bquote(bold(.(toupper(paste0(colnames(dat)[i],'    '))))), cex = 1.25, bg = 'gray90')
}

## ################################ ##
## Plot Manhattan single chromosome ##
## ################################ ##

plotManhattan_chrom <- function(dat = NULL, chrom = NULL, gene = NULL, legend.pos = NULL) {
  # subset the dataframe with specific population and specific chromosome
    data <- dat[dat$CHRHomoeo == chrom, c(1:4, grep(gene, colnames(dat)))]
  # set up ylim.max limit
    ylim.max = min(apply(h13.combined[, -c(1:4)], 2, function(p) max(-log10(na.omit(p)))))
  # plot manhattan
    manhattan(data, p = colnames(data)[5], suggestiveline = F, genomewideline = F, ylim = c(0, ylim.max + 15),
              ylab='', xlab='', cex = 1, cex.axis = 1, xaxt = 'n')
    abline(h = -log10(alpha / nrow(data)), lty = 3, col = 'red')
    axis(side = 1, at = seq(0, 900, 100) * 10^6, labels = seq(0, 900, 100))
    mtext(text = '-log10(p)', 2, line = 2.7, cex = 1)
    mtext(text = paste0('Chromosome ', chrom, ' (Mb)'), side = 1, line = 3, cex = 1)

    for (i in 6:ncol(data)) {
      points(data$BP, -log10(data[, i]), col = i-4)
    }
    
    legend(legend.pos, legend = toupper(paste0(colnames(data)[-c(1:4)],'    ')), 
           pch = c(16,1,1,1), col = c(1,2,3,4), cex = 0.75)
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
