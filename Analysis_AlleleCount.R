## ###########################
## Hessian Fly TAG Manuscript
## Author - Narinder Singh
## Mapping with Allele counts
## ###########################

# load required functions
  source("functions.R")

# read allele counts from GBS data
  alleleCounts <- fread('data/AlleleCounts.txt', header = T, check.names = F, data.table = F)
  
# read allele counts from RenSeq Data
  alleleCounts_RenSeq_Av1 <- fread('data/AC_RenSeq-Av1.txt', header = T, check.names = F, data.table = F)
  alleleCounts_RenSeq_Tv2 <- fread('data/AC_RenSeq-Tv2.txt', header = T, check.names = F, data.table = F)
  
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
  
  ## Gene H26 families
      # KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts, pop.id = 'KO', gene = 'h26')
      # RenSeq Av1 Family 1 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Av1, pop.id = 'fam_1', gene = 'h26.renseq.av1.fam1')
      # RenSeq Av1 Family 2 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Av1, pop.id = 'fam_2', gene = 'h26.renseq.av1.fam2')
      # RenSeq Tv2 Family 1 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_1', gene = 'h26.renseq.tv2.fam1')
      # RenSeq Tv2 Family 2 - KU2147 x Overley - H26
        associationTest_AC(data = alleleCounts_RenSeq_Tv2, pop.id = 'fam_2', gene = 'h26.renseq.tv2.fam2')
      
  # SynOP DH - H32
  associationTest_AC(data = alleleCounts, pop.id = 'SynOp', gene = 'h32')
  