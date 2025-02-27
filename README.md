## This repo contains code for mapping and analyzing Hessian fly resistance genes

### software dependencies ===
- GNU parallel
- samtools
- bcftools
- bowtie2
- Tassel 5 GBSv2 Pipeline
- R and perl

### workflow ===
1. The `*.txt` files in `computation_scripts` folder are the keyfiles containing sample information for their corresponding `*.sh` shell scripted pipelines. These pipelines are used to get variant calls and allele counts for GBS data that are further used in mapping. File `03_script-HFly_RenSeq_AC_Tv2.sh` is used to get allele counts from RenSeq data. The `R` script `03_script-ToFixChromPos_RenSeq_AC_Tv2.R` is called from within the previous script to fix chromosome number and positions from using the pseudomolecule assembly. Finally, the perl script `04_vcf2AC_v2.pl` is used to extract allele count information from vcf file.

2. Once the variants and allele count data is generated, the `*.R` scripts in `src` folder can be used for further analyses.
