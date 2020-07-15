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
  pdf(file="output/Fig.S1_HF_PlantColor.pdf", width = 11, height = 5.5)
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

