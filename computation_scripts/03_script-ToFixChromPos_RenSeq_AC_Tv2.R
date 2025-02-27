library(data.table)

chrom.info <- read.table('required_files/chrom_info.txt', header = T, as.is = T)

RenSeq_AC_parts <- fread("variants/HFly_RenSeq_AC_Tv2_parts.txt", 
                         header = T, check.names = F, stringsAsFactors = F, data.table = F)

RenSeq_AC_parts <- RenSeq_AC_parts[RenSeq_AC_parts$CHR != "chrUn", ]
RenSeq_AC_parts$CHR <- sub("_part1", "", RenSeq_AC_parts$CHR)

RenSeq_AC_full <- merge(RenSeq_AC_parts[, c("CHR", "POS")], 
                        chrom.info[, c("part_name", "addThis")], 
                        by.x = "CHR", by.y = "part_name", all.x = T, 
                        incomparables = chrom.info$chrom)

RenSeq_AC_full$addThis[is.na(RenSeq_AC_full$addThis)] <- 0L

RenSeq_AC_full$newPOS <- RenSeq_AC_full$POS + RenSeq_AC_full$addThis

RenSeq_AC_full <- RenSeq_AC_full[, c("CHR", "POS", "newPOS")]

RenSeq_AC_full <- merge(RenSeq_AC_full, RenSeq_AC_parts)

RenSeq_AC_full <- RenSeq_AC_full[, !colnames(RenSeq_AC_full) %in% "POS"]

colnames(RenSeq_AC_full)[colnames(RenSeq_AC_full) %in% "newPOS"] <- "POS"

RenSeq_AC_full <- RenSeq_AC_full[order(RenSeq_AC_full[, "CHR"], RenSeq_AC_full[, "POS"]), ]

RenSeq_AC_full$CHR <- sub("_part2", "", RenSeq_AC_full$CHR)

RenSeq_AC_full$CHR <- sub("chr", "", RenSeq_AC_full$CHR)

fwrite(RenSeq_AC_full, file = "variants/HFly_RenSeq_AC_Tv2.vcf2AC.txt", sep = "\t")
