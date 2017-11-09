require(limma)

cbred <- 2 # '#D55E00'
cbblue <- '#0072B2'
cborange <- '#E69F00'
cbgreen <- '#009E73'
cbpink <- '#CC79A7'
cblblue <- '#56B4E9'
cdyellow <- '#F0E442'

args <- commandArgs(trailingOnly = TRUE)
sub_analyse = paste(args[1])
FDR = as.numeric(paste(args[2]))
logFC = as.numeric(paste(args[3]))

# sub_analyse <- 'B_FB'
# FDR <- 0.1
# logFC <- 0

outpath <- paste("~/git/feminisation_direction/output/courtship_bias_analysis/", sep="")

mcgraw_data <- read.table("~/git/feminisation_direction/input/McGraw_genes.txt", header=T)
str(mcgraw_data)

dpse_dmel_data <- read.table('/Users/zabameos/git/feminisation_direction/input/dpseDmel.txt', header=T)
str(dpse_dmel_data)

c_de_data <- read.table(paste("/Users/zabameos/git/feminisation_direction/output/courtship_bias_analysis/", sub_analyse, "/LogCPM_", FDR, '_', sub_analyse, ".txt", sep=''), header=T) 
str(c_de_data)

mcgraw_annotation <- read.table('/Users/zabameos/git/feminisation_direction/input/McGraw_annotation.txt', header=T)

merged <- merge(mcgraw_data, dpse_dmel_data, by.x='FBID', by.y='dmel', all=F)
merged2 <- merge(merged, c_de_data, by.x='dpse', by.y='gene', all=T)
deTable <- subset(merged2, merged2$logCPM!='NA') 

deTable$category <- as.character(deTable$category) 
deTable$category[is.na(deTable$category)] <- 'NA' 
deTable$category <- as.factor(deTable$category)

merged2$category <- as.character(merged2$category) 
merged2$category[is.na(merged2$category)] <- 'NA' 
merged2$category <- as.factor(merged2$category)

source('~/git/feminisation_direction/scripts/vennFunctions.r')

set1 <- as.character(factor(merged2$dpse))
set2 <- as.character(factor(merged2$dpse[merged2[14] < FDR & abs(merged2[10]) > logFC ]))
set3 <- as.character(factor(merged2$dpse[merged2$category!='NA']))

set2_up <- as.character(factor(merged2$dpse[merged2[14] < FDR & abs(merged2[10]) > logFC & merged2[10] > 0 ]))
set2_down <- as.character(factor(merged2$dpse[merged2[14] < FDR & abs(merged2[10]) > logFC & merged2[10] < 0 ]))

set3_up <- as.character(factor(merged2$dpse[merged2$category!='NA' & merged2$updown=='up']))
set3_up <- subset(set3_up, set3_up!='NA')
set3_down <- as.character(factor(merged2$dpse[merged2$category!='NA' & merged2$updown=='down']))
set3_down <- subset(set3_down, set3_down!='NA')



pdf(file.path(outpath, paste('Venn_3_McGraw_FDR_',FDR, '_logFC_', logFC, '_', sub_analyse,'.pdf', sep="")), width=12, height=12)
par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3), oma=c(0,0,4,0))
Venn3(set1, set2, set3, c('All', 'DE', 'McGraw'), 'All genes')
Venn4(set2_up, set2_down, set3_up, set3_down, c('Virgin up', 'Courted up', 'McGraw up', 'McGraw down'), 'DE only')
Venn5(set1, set2_up, set2_down, set3_up, set3_down, c('All', 'Virgin up', 'Courted up', 'McGraw up', 'McGraw down'), 'All by up and down')

mtext(paste(sub_analyse, ' FDR = ', FDR, ' - logFC = ', logFC, sep=""), outer = TRUE, cex = 1)
dev.off()


set1 <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='NA' ]))
set2 <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='sperm' ]))
set3 <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='acp' ]))
set4 <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='other' ]))

set1_up <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='NA' & deTable[10] > 0]))
set2_up <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='sperm' & deTable[10] > 0]))
set3_up <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='acp' & deTable[10] > 0]))
set4_up <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='other' & deTable[10] > 0]))

set1_down <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='NA' & deTable[10] < 0]))
set2_down <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='sperm' & deTable[10] < 0]))
set3_down <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='acp' & deTable[10] < 0]))
set4_down <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category=='other' & deTable[10] < 0]))

set2_upup <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category!='NA' & deTable[10] > 0 & deTable$updown=='up' ]))
set2_updown <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category!='NA' & deTable[10] > 0 & deTable$updown=='down' ]))
set2_downup <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category!='NA' & deTable[10] < 0 & deTable$updown=='up' ]))
set2_downdown <- as.character(factor(deTable$dpse[deTable[14] < FDR & abs(deTable[10]) > logFC & deTable$category!='NA' & deTable[10] < 0 & deTable$updown=='down' ]))

pdf(file.path(outpath, paste('Venn_4_McGraw_FDR_',FDR, '_logFC_', logFC, '_', sub_analyse,'.pdf', sep="")), width=16, height=6)
par(mfrow=c(1,3)) 
par(mar=c(5,5,4,3), oma=c(0,0,4,0))
Venn4(set1, set2, set3, set4, c('notMcGraw', 'sperm', 'acp', 'other'), 'DE - all')
Venn4(set1_up, set2_up, set3_up, set4_up, c('notMcGraw', 'sperm', 'acp', 'other'), 'Virgin - up')
Venn4(set1_down, set2_down, set3_down, set4_down, c('notMcGraw', 'sperm', 'acp', 'other'), 'Courted - up')
mtext(paste(sub_analyse, ' FDR = ', FDR, ' - logFC = ', logFC, sep=""), outer = TRUE, cex = 1)
dev.off()

overlapping_genes <- subset(merged2, merged2[14] < FDR & abs(merged2[10]) > logFC & merged2$category!='NA') 
annotated_overlapping <- merge(overlapping_genes, mcgraw_annotation, by.x='FBID', by.y='FBID', all=F)

de_nomac <- length(as.character(factor(merged2$dpse[merged2[14] < FDR & abs(merged2[10]) > logFC & merged2$category=='NA' ])))
de_mac <- length(as.character(factor(merged2$dpse[merged2[14] < FDR & abs(merged2[10]) > logFC & merged2$category!='NA' ])))
node_nomac <- length(as.character(factor(merged2$dpse[merged2[14] > FDR & abs(merged2[10]) > logFC & merged2$category=='NA' ])))
node_mac <- length(as.character(factor(merged2$dpse[merged2[14] > FDR & merged2$category!='NA' ])))

# FB
de_nomac <- 280
de_mac <- 40
node_nomac <- 10075
node_mac <- 1091

# MH
de_nomac <- 602
de_mac <- 58
node_nomac <- 11262
node_mac <- 1073

total_significant <- de_nomac + de_mac
total_mac <- sum(de_mac + node_mac)
total_nomac <- sum(de_nomac + node_nomac)
exp_mac <- round(total_significant * total_mac/(total_mac + total_nomac))
exp_nomac <- round(total_significant * total_nomac/(total_mac + total_nomac))


Obs <- c(de_mac, de_nomac)
Exp <- c(exp_mac, exp_nomac)
chi_obs_data <- rbind(Obs, Exp)
colnames(chi_obs_data) <- c('McGraw', 'not McGraw')
chi_test <- chisq.test(c(de_mac, de_nomac), p=c(total_mac/(total_mac + total_nomac), total_nomac/(total_mac + total_nomac)))
pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

cbblue <- '#0072B2'
cblblue <- '#56B4E9'

par(mar=c(5,5,4,3), oma=c(0,0,2,0))
mp <- barplot(chi_obs_data, col=c(cbgreen ,cbblue), main = paste('DE ', sub_analyse, '\nFDR=',FDR, ' logFC=', logFC, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
names(vals) <- LETTERS[1:4]
text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)

chi_test
