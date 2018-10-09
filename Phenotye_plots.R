## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Phenotype plots
## ###########################

## load functions
  source('functions.R')

## Hessian fly screening in greenhouse
  pdf(file="output/phenotypic_distribution.pdf", width = 8.5, height = 11)
  par(mfrow=c(4,2))

## H5
  phenoPlots(c(18,63,65), "H5")
## H6
  phenoPlots(c(23,73,26), "H6")
## H12
  phenoPlots(c(7,46,20), "H12")
## H13.Newton
  phenoPlots(c(31,56,26), "H13 Newton")
## H13.Overley
  phenoPlots(c(19,39,29), "H13 Overley")
## H26.Fam1
  phenoPlots(c(68,177,71), "H26 Fam1")
## H26.Fam2
  phenoPlots(c(66,137,70), "H26 Fam2")
## H32
  phenoPlots(c(37,0,43), "H32")

  dev.off()

## Plant color distribution in field
  pdf(file="output/plant_color_field.pdf", width = 11, height = 5.5)
  par(mfrow=c(1,2))
  barplot(c(72,151,81), col = c('forestgreen', 'gray', 'green'), 
          names.arg = c('Overley like', 'Segregating', 'Ae. tauschii like'), 
          cex.axis = 1.5, cex.names = 1.2,
          main = "Fam1 lines distribution for spike color")

# test of independence of color and hessian fly
  hfColData <- read.csv('data/HF_color_data.csv', header = T)
  categData <- table(hfColData$hf, hfColData$color)
  categData
  
  chisq.test(categData)
  spineplot(categData, xaxlabels = c('Rr', 'RR', 'rr'), yaxlabels = c('Overley like', 'Segregating', 'Ae. tauschii like'),
            col = c('forestgreen', 'gray', 'green'), main = 'Color distribution within different genotypes')
  mtext(paste('Chi-sq p-val =', round(chisq.test(categData)$p.value, 3)), side = 3, cex = 1.25)
  
  dev.off()

  