# load packages
library(pacman)
p_load(reshape, dplyr)

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

