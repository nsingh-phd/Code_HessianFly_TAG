## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Mapping with Allele counts
## ###########################

# load required functions, files, and packages
  chrom.info <- read.table('data/chrom_info.txt', header = T, as.is = T)
  source("functions.R")

# read allele counts
  alleleCounts <- fread('data/AlleleCounts.txt', header = T, check.names = F, data.table = F)

## association tests
  
  # # Carol x Newton - H3
  associationTest_AC(data = alleleCounts, pop.id = 'CN', gene = 'h3')
  
  # Erin x Newton - H5
  associationTest_AC(data = alleleCounts, pop.id = 'EN', gene = 'h5')
  
  # Flynn x Newton - H6
  associationTest_AC(data = alleleCounts, pop.id = 'FN', gene = 'h6')
  
  # Joy x Newton - H10
  associationTest_AC(data = alleleCounts, pop.id = 'JN', gene = 'h10')
  
  # Lola x Newton - H12
  associationTest_AC(data = alleleCounts, pop.id = 'LN', gene = 'h12')
  
  # Molly x Newton - H13
  associationTest_AC(data = alleleCounts, pop.id = 'MN', gene = 'h13.newton')
  
  # Molly x Overley - H13
  associationTest_AC(data = alleleCounts, pop.id = 'MO', gene = 'h13.overley')
  
  # KU2147 x Overley - H26
  associationTest_AC(data = alleleCounts, pop.id = 'KO', gene = 'h26')
  
  # SynOP DH - H32
  associationTest_AC(data = alleleCounts, pop.id = 'SynOp', gene = 'h32')
  