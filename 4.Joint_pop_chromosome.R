####################
# full chromosome 6D
####################

pdf('Fig.4_H10_H13.pdf', width = 11, height = 6.5)

molly_6d.newton2 = molly_6d.newton
molly_6d.newton2$BP = molly_6d.newton2$BP/1000000

# molly newton
manhattan(molly_6d.newton2, suggestiveline = F, genomewideline = F,
          xlab = 'Chromosome 6D position in Mb', xlim = range(0, 510),
          cex.lab=1.4)
abline(h = -log10(alpha/nrow(allchr_bwa)), lty = 3, col = 'red')

# molly overley

points(molly_6d.overley$BP/1000000, -log10(molly_6d.overley$P), col=2)

# joy newton
points(joy_6d$BP/1000000, -log10(joy_6d$P), col=3)

# centromere
abline(v=214085311/1000000, col='gray', lwd=7) #6D - Molly and Joy
text(x = 214085400/1000000, y = 9, labels = 'Centromere', srt=90, pos=4, offset = 0.7, cex = 1)
text(x = 220086000/1000000, y = 7.5, labels = 'Bonferroni threshold', pos=4, offset = 0.7, cex = 1)

# legend
legend('topright', legend = c('H13-MN', 'H13-MO', 'H10-JN', 'Bonferroni threshold'), col = c(1, 2, 3, 2),
       pch = c(20,1,1,NA), lty = c(NA,NA,NA,3), cex = 1)

dev.off()


###########################
# portion of the chromosome
###########################
pdf('molly_joy_6DS_part.pdf', width = 11, height = 6.5)

# molly newton
manhattan(molly_6d.newton, suggestiveline = F, genomewideline = -log10(alpha/nrow(allchr_bwa)), cex.lab=1.4, xlim=c(15000000, 33039120))

# molly overley
points(molly_6d.overley$BP, -log10(molly_6d.overley$P), col=2)

# joy newton
points(joy_6d$BP, -log10(joy_6d$P), col=3)

# legend
legend('topright', legend = c('H13-MN', 'H13-MO', 'H10-JN'), col = c(1, 2, 3), pch = c(20,1,1), cex = 1)

dev.off()

#####################