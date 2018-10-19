## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Phenotype plots
## ###########################

## load functions
  source('HessianFly_TAG_functions.R')


##
## Hessian fly screening in greenhouse
##
  pdf(file="output/Fig.1_phenotypic_distribution.pdf", width = 8.5, height = 11)
  par(mfrow=c(4,2))

## H5
  phenoPlots(values = c(18,63,65), name = "H5-EN", seg = "1:2:1", trait = 'HF')
## H6
  phenoPlots(values = c(23,73,26), name = "H6-FN", seg = "1:2:1", trait = 'HF')
## H12
  phenoPlots(values = c(7,46,20), name = "H12-LN", seg = "1:2:1", trait = 'HF')
## H13.Newton
  phenoPlots(values = c(31,56,26), name = "H13-MN", seg = "1:2:1", trait = 'HF')
## H13.Overley
  phenoPlots(values = c(19,39,29), name = "H13-MO", seg = "1:2:1", trait = 'HF')
## H26.Fam1
  phenoPlots(values = c(68,177,71), name = "H26-Fam1", seg = "1:2:1", trait = 'HF')
## H26.Fam2
  phenoPlots(values = c(66,137,70), name = "H26-Fam2", seg = "1:2:1", trait = 'HF')
## H32
  phenoPlots(values = c(37,0,43), name = "H32-SynOpDH", seg = "1:1", trait = 'HF')

  dev.off()

##
## Plant color distribution in field
##
  pdf(file="output/S1_HF_PlantColor.pdf", width = 11, height = 5.5)
  par(mfrow=c(1,2))
  phenoPlots(values = c(72,151,81), name = "H26-Fam1", seg = "1:2:1", trait = 'plant_color')

  # test of independence of color and hessian fly
  hfColData <- read.csv('data/HF_color_data.csv', header = T)
  categData <- table(hfColData$hf, hfColData$color)
  categData
  
  chisq.test(categData)
  spineplot(categData, xaxlabels = c('Rr', 'RR', 'rr'), 
            yaxlabels = c("'Overley'-like", 'Het', 'Ae. tauschii like'),
            col = c('forestgreen', 'gray', 'green'), 
            main = 'Color distribution within different genotypes')
  mtext(paste('Chi-sq p-val =', round(chisq.test(categData)$p.value, 3)), side = 3, cex = 1.25)
  
  dev.off()

##
## ANOVA for agronomic traits
##
  
  #read the field trial data
  pheno <- read.csv('data/2016-17_field_data.csv', header = T, as.is = T)
  
  #use extra cols
  colnames(pheno)[4] = 'family'
  colnames(pheno)[5] = 'res_sus'
  pheno$family = pheno$res_sus = NA
  
  #update extra cols
  pheno$family[grep('fam1', pheno$entry, ignore.case = T)] = '1'
  pheno$family[grep('fam2', pheno$entry, ignore.case = T)] = '2'
  
  pheno$res_sus[grep('res', pheno$entry, ignore.case = T)] = 'HF Res'
  pheno$res_sus[grep('sus', pheno$entry, ignore.case = T)] = 'HF Sus'
  
  pheno$res_sus[grep('Everest', pheno$entry, ignore.case = T)] = 'Everest'
  pheno$res_sus[grep('Overley', pheno$entry, ignore.case = T)] = 'Overley'
  pheno$res_sus[grep('SY', pheno$entry, ignore.case = T)] = 'SY Flint'
  
  # family info
  fam.info <- pheno[!duplicated(pheno$entity_id), c(1, 4:5)]
  
  # rehsape data
  pheno <- reshape(data = pheno[, 1:3], idvar = "entity_id", timevar = "trait_id", direction = "wide")
  
  pheno <- merge(fam.info, pheno, by = 1, sort = F)
  
  ## trait analysis
  unique(pheno$trait_id)
  
  # adjusting grain weight based on missing rows
  pheno <- cbind(pheno, adj.grainwt=NA)
  pheno$adj.grainwt = (as.numeric(pheno$phenotype_value.GRWT) / 
                         (6 - as.numeric(pheno$phenotype_value.NAROWS))) * 6
  pheno$adj.grainwt = pheno$adj.grainwt / 1000 # g to kg
  
  
  # adjusting days to heading (relative to Jan 1, 2017)
  pheno$phenotype_value.DTHD = as.numeric(pheno$phenotype_value.DTHD) - 73
  
  # remove Everest and SY Flint
  pheno <- pheno[!pheno$res_sus %in% c('Everest', 'SY Flint'), ]
  
  # ANOVA 
  summary.trait <- function(trait=trait, predictor=predictor, dataset=data) {
    cat(trait)
    trait.mean.sd <- group_by(.data = dataset, get(predictor)) %>%
      summarise(
        count = n(),
        mean = mean(as.numeric(get(trait)), na.rm = TRUE),
        sd = sd(as.numeric(get(trait)), na.rm = TRUE))
    fit.aov <- aov(formula = as.numeric(get(trait)) ~ get(predictor), data = dataset)
    return.list <- list("anova.summary"=summary(fit.aov), "mean.sd"=trait.mean.sd)
    return(return.list)
  }
  
  summary.trait(trait = 'adj.grainwt', predictor = 'res_sus', dataset = pheno)
  summary.trait(trait = 'phenotype_value.TESTWT', predictor = 'res_sus', dataset = pheno)
  summary.trait(trait = 'phenotype_value.MOIST', predictor = 'res_sus', dataset = pheno)
  summary.trait(trait = 'phenotype_value.PTHT_avg', predictor = 'res_sus', dataset = pheno)
  summary.trait(trait = 'phenotype_value.DTHD', predictor = 'res_sus', dataset = pheno)
  
  pdf('output/Fig.2_FieldTest.pdf', height = 4.5, width = 8.5)
  par(mfrow=c(1,2), oma=c(0,2,0,2))
  trait.plot(data = pheno, trait = "adj.grainwt", trait.name = 'Adjusted grain weight (kg)', plot = 'barplot')
  trait.plot(data = pheno, trait = "adj.grainwt", trait.name = 'Adjusted grain weight (kg)', plot = 'hist')
  dev.off()
  
  pdf('output/S2_FieldTests.pdf', height = 11, width = 8.5)
  par(mfrow=c(4,2), oma=c(0,6,0,6))
  trait.plot(data = pheno, trait = "moist", trait.name = 'Grain moisture (%)', plot = 'barplot')
  trait.plot(data = pheno, trait = "moist", trait.name = 'Grain moisture (%)', plot = 'hist')
  trait.plot(data = pheno, trait = "ptht_avg", trait.name = 'Plant height (cm)', plot = 'barplot')
  trait.plot(data = pheno, trait = "ptht_avg", trait.name = 'Plant height (cm)', plot = 'hist')
  trait.plot(data = pheno, trait = "testwt", trait.name = 'Test weight (lbs./bu.)', plot = 'barplot')
  trait.plot(data = pheno, trait = "testwt", trait.name = 'Test weight (lbs./bu.)', plot = 'hist')
  trait.plot(data = pheno, trait = "dthd", trait.name = 'Days to heading (days)', plot = 'barplot')
  trait.plot(data = pheno, trait = "dthd", trait.name = 'Days to heading (days)', plot = 'hist')
  dev.off()
  
  