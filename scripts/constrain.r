library(edgeR)
library(gplots)
library(ggplot2)
library(dynamicTreeCut)

cbred <- 2 # '#D55E00'
cbblue <- '#0072B2'
cborange <- '#E69F00'
cbgreen <- '#009E73'
cbpink <- '#CC79A7'
cblblue <- '#56B4E9'
cdyellow <- '#F0E442'
col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.15)
col1b <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)
col3 <- rgb(red = 0, green = 1, blue = 0, alpha = 0.6)
col4 <- rgb(red = 0, green = 0, blue = 1, alpha = 0.6)
col5 <- rgb(red = 0.7, green = 0, blue = 0.5, alpha = 0.6)

# Base file naming on script input
args <- commandArgs(trailingOnly = TRUE)
ssline1 = paste(args[1])
ssline2 = paste(args[2])
gender = paste(args[3])
courtship = paste(args[4])
tissue = paste(args[5] )
FDR2use = as.numeric(paste(args[6]))

sexbthres <- 1 # threshold logFC for sex biased genes

# ssline1 <- "E"
# ssline2 <- "M"
# gender <- "M"
# courtship <- "V"
# tissue <- "B"
# FDR2use  <- 0.1

courtcomp <- paste('EMB_',gender,tissue,sep="")
compname <- paste(ssline1,ssline2,gender,courtship,tissue,sep="")

datapath <- "~/git/feminisation_direction/input"
scriptpath <- "~/git/feminisation_direction/scripts/"
outpath <- "~/git/feminisation_direction/output/subsets"
outpath2 <- paste(outpath,'/',compname,sep="")
dir.create(file.path(outpath))
dir.create(file.path(outpath2))

load(file.path(datapath, "designtab.Rdata")) # table of design information 
designtab
count <- read.table(file.path(datapath, 'new_count.txt'), header=T)

# courtship analysis
EMB_MH <- c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93)
EMB_MB <- EMB_MH+1
EMB_FH <- EMB_MH+2
EMB_FB <- EMB_MH+3

if (courtcomp == 'EMB_MH') { 
	EMB_count <- data.frame(count[,EMB_MH])
	EMB_designtab <- designtab[EMB_MH,]
}
if (courtcomp == 'EMB_MB') { 
	EMB_count <- data.frame(count[,EMB_MB])
	EMB_designtab <- designtab[EMB_MB,]
}
if (courtcomp == 'EMB_FH') { 
	EMB_count <- data.frame(count[,EMB_FH])
	EMB_designtab <- designtab[EMB_FH,]
}
if (courtcomp == 'EMB_FB') { 
	EMB_count <- data.frame(count[,EMB_FB])
	EMB_designtab <- designtab[EMB_FB,]
}

# load(file.path(datapath, "neworthanno.Rdata"))
# anno <- as.matrix(neworthanno)

annotation <- read.delim(file.path(datapath, "exongeneanno.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

# EMB libraries are used (B is retained) so that the low variance of B helps in dispersion estimation for genes.

model.formula <- as.formula("~0+group")
EMB_dmat <- model.matrix(model.formula, data=as.data.frame(EMB_designtab))

EMB_group <- as.character(EMB_designtab[,5])
EMB_dgl <- DGEList(counts <- as.matrix(EMB_count), group=EMB_group, genes=annotation, remove.zeros=T)

filter_file <- file.path(paste(outpath, '/',courtcomp, 'filtering_info.txt', sep=""))
ave0 <- EMB_dgl[aveLogCPM(EMB_dgl) >= 0,]
ave1 <- EMB_dgl[aveLogCPM(EMB_dgl) >= 1,]
sum2 <- EMB_dgl[rowSums(cpm(EMB_dgl)) >= 2,]
sum10 <- EMB_dgl[rowSums(cpm(EMB_dgl)) >= 10,]

write(paste("Comp\tRemaining\tMin\t1Q\tMedian\tMean\t3Q\tMax"), filter_file)
write(paste("All_"), filter_file, append=T)
write(paste(nrow(EMB_dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(EMB_dgl)/ncol(EMB_dgl))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave0_"), filter_file, append=T)
write(paste(nrow(ave0), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave0)/ncol(ave0))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave1_"), filter_file, append=T)
write(paste(nrow(ave1), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave1)/ncol(ave1))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum2_"), filter_file, append=T)
write(paste(nrow(sum2), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum2)/ncol(sum2))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum10_"), filter_file, append=T)
write(paste(nrow(sum10), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum10)/ncol(sum10))), filter_file, append=T, sep='\t', ncol=6)

# compile nice table with perl -pe 's/_\n/\t/g' #_filtering_info.txt

# "Genes were filtered from the analysis if their average log2 count per million (as computed by edgeR's aveLogCPM function) was negative. This had the effect of keeping genes with an average count of about 50 or more per sample."

# to check effect
# EMB_dgl <- EMB_dgl[aveLogCPM(EMB_dgl) >= 0,]

EMB_dgl <- EMB_dgl[aveLogCPM(EMB_dgl) >= 0,] # keep average cpm > 0 for further analysis

xcpm <- mglmOneGroup(EMB_dgl$counts)        # computing a logCPM for making dispersion plot

# estimate data normlization factors and dispersion
EMB_dgl <- calcNormFactors(EMB_dgl)
EMB_dgl <- estimateGLMCommonDisp(EMB_dgl,EMB_dmat)
EMB_dgl <- estimateGLMTrendedDisp(EMB_dgl,EMB_dmat,min.n=1000)
EMB_dgl <- estimateGLMTagwiseDisp(EMB_dgl,EMB_dmat)

y <- EMB_dgl
colnames(y) <- paste(colnames(y), EMB_group,sep="\n")

## MDS plot
pdf(file.path(outpath,paste('MDS_', courtcomp, '.pdf', sep="")), width=8, height=8) 
par(mar=c(5,5,4,3))
plotMDS(y, cex=0.4, col=as.numeric(y$samples$group), main=paste(courtcomp, "_MDSplot", sep=""))
# legend("bottomleft", as.character(unique(y$samples$group)), col=1:6, pch=20, cex=0.5)
dev.off()

## Dispersion plot
pdf(file.path(outpath, paste('Dispersion_', courtcomp, '.pdf', sep="")), width=8, height=8)
plot(xcpm, EMB_dgl$tagwise.dispersion,pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste(courtcomp," dispersion", sep=""))
if(!is.null(EMB_dgl$trended.dispersion)) points(xcpm ,EMB_dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
points(xcpm ,EMB_dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=EMB_dgl$common.dispersion, col=cblblue, lwd=2)
dev.off()
##  fit the data model ------------------------------------------
fitres <- glmFit(EMB_dgl, EMB_dmat)

dmat <- EMB_dmat
count <- EMB_count

# load yongxiang function
source(file.path(scriptpath, "mkmaincontrast.R"))

## Female-virgin-2Tissues
contrast=list(fct="R", lev = c(ssline1,ssline2))
constraint=list(fct=c("G","S","T"), lev=c(gender, courtship, tissue)) # Gender, courtship Status, Tissue 

c1mat <- mkmaincontrast(dmat, contrast=contrast, constraint=constraint)

allzeros <- which(rowSums(count)==0)

comparison <- glmLRT(fitres, contrast=c1mat)
PV <- comparison$table[,"PValue"]
FDR <- p.adjust(comparison$table[,"PValue"], method="BH")
logFC <- comparison$table[,"logFC"]
xcpm <- comparison$table[,"logCPM"]
resultscomparison <- cbind(rownames(comparison$table), comparison$table, FDR)
colnames(resultscomparison) <- c("gene", "logFC", "logCPM", "LR", "PValue", "FDR")
resultscomparisonannotated <- merge(resultscomparison, annotation, by.x='gene', by.y='gid', all=F)
write.table(resultscomparisonannotated, file=file.path(outpath2,paste(compname,'_all.txt', sep="")), quote=F, row.names=F, sep="\t")

significant_comparison <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR < FDR2use) 
# nrow(significant_comparison)
up_significant_comparison <- subset(resultscomparisonannotated, resultscomparison$FDR < FDR2use & resultscomparisonannotated$logFC > 0) 
# nrow(up_significant_comparison)
down_significant_comparison <- subset(resultscomparisonannotated, resultscomparison$FDR < FDR2use & resultscomparisonannotated$logFC < 0) 
# nrow(down_significant_comparison)

write.table(sort(significant_comparison$gene, decreasing=T), file=file.path(outpath2,paste(compname,'_dpse_gIDs.txt', sep="")), quote=F, row.names=F, sep="\t")

## B only sex biased lists

if (courtship!='All') {
sexbiasdata <- read.table(file.path(datapath, paste('sex_bias_B_', courtship, tissue, '.txt', sep="")), header=T)
sexbiasdata$sexbias <- 'neutral'
sexbiasdata$sexbias[sexbiasdata[2] < -sexbthres & sexbiasdata[3] < 0.05] <- 'female' 
sexbiasdata$sexbias[sexbiasdata[2] > sexbthres & sexbiasdata[3] < 0.05] <- 'male' 
maleDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='male') 
femaleDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='female') 
neutralDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='neutral') 
}

if (courtship=='All') {
sexbiasdata <- read.table(file.path(datapath, paste('sex_bias_B_', tissue, '.txt', sep="")), header=T)
sexbiasdata$sexbias <- 'neutral'
sexbiasdata$sexbias[sexbiasdata[2] < -sexbthres & sexbiasdata[3] < 0.05] <- 'female' 
sexbiasdata$sexbias[sexbiasdata[2] > sexbthres & sexbiasdata[3] < 0.05] <- 'male' 
maleDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='male') 
femaleDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='female') 
neutralDpse <- subset(sexbiasdata, sexbiasdata$sexbias=='neutral') 
}



## export for Hollis comparison, normalised counts
d = estimateCommonDisp(EMB_dgl, verbose=TRUE)
pseu_counts_EMB_dgl <- d$pseudo.counts
pseu_counts_EMB_dgl <- data.frame(row.names(pseu_counts_EMB_dgl), pseu_counts_EMB_dgl)
colnames(pseu_counts_EMB_dgl) <- c('ID',EMB_group)
pseu_counts_EMB_dgl2 = merge(resultscomparison, pseu_counts_EMB_dgl, by.x="gene", by.y="ID", all=F )
pseu_counts_EMB_dgl2$bias <- 'empty'

pseu_counts_EMB_dgl2$bias[match(pseu_counts_EMB_dgl2$gene, femaleDpse$geneID) !='NA'] <- 'female'
pseu_counts_EMB_dgl2$bias[match(pseu_counts_EMB_dgl2$gene, maleDpse$geneID) !='NA'] <- 'male'
pseu_counts_EMB_dgl2$bias[match(pseu_counts_EMB_dgl2$gene, neutralDpse$geneID) !='NA'] <- 'neutral'
pseu_counts_EMB_dgl2$bias <- as.factor(pseu_counts_EMB_dgl2$bias)

write.table(pseu_counts_EMB_dgl2, file=file.path(outpath2, paste(compname,'_pseudo_counts.txt', sep="")), quote=F, row.names=F, sep="\t")



## export for bootstaps
### calculation of moderated log-counts-per-million for use in heatmaps
nc <- cpm(EMB_dgl, prior.count=2, log=T)
nc2 <- data.frame(row.names(nc), nc)

colnames(nc2) <- c('ID', EMB_group)
resultscomparison_mod_logcpm = merge(resultscomparison, nc2, by.x="gene", by.y="ID", all=F )
resultscomparison_mod_logcpm$bias <- 'empty'

resultscomparison_mod_logcpm$bias[match(resultscomparison_mod_logcpm$gene, femaleDpse$geneID) !='NA'] <- 'female'
resultscomparison_mod_logcpm$bias[match(resultscomparison_mod_logcpm$gene, maleDpse$geneID) !='NA'] <- 'male'
resultscomparison_mod_logcpm$bias[match(resultscomparison_mod_logcpm$gene, neutralDpse$geneID) !='NA'] <- 'neutral'
resultscomparison_mod_logcpm$bias <- as.factor(resultscomparison_mod_logcpm$bias)
row.names(resultscomparison_mod_logcpm) <- resultscomparison_mod_logcpm$gene
nrow(subset(resultscomparison_mod_logcpm, resultscomparison_mod_logcpm$bias=='empty'))

write.table(resultscomparison_mod_logcpm, file=file.path(outpath2, paste('Results_long_FDR_', FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep="\t")

female_biased=merge(femaleDpse, significant_comparison, by.x="geneID", by.y="gene", all=F )
# nrow(female_biased)
male_biased=merge(maleDpse, significant_comparison, by.x="geneID", by.y="gene", all=F )
# nrow(male_biased)
unbiased=merge(neutralDpse, significant_comparison, by.x="geneID", by.y="gene", all=F )
# nrow(unbiased)
comparison_numbers <- cbind(nrow(female_biased), nrow(male_biased), nrow(unbiased))
colnames(comparison_numbers) <- c("female", "male", "unbiased")
rownames(comparison_numbers) <- "comparison"
nr_fem_positive <- nrow(data.frame(female_biased$logFC[female_biased$logFC > 0]))
nr_fem_negative <- nrow(data.frame(female_biased$logFC[female_biased$logFC < 0]))
nr_mal_positive <- nrow(data.frame(male_biased$logFC[male_biased$logFC > 0]))
nr_mal_negative <- nrow(data.frame(male_biased$logFC[male_biased$logFC < 0]))
nr_neutr_positive <- nrow(data.frame(unbiased$logFC[unbiased$logFC > 0]))
nr_neutr_negative <- nrow(data.frame(unbiased$logFC[unbiased$logFC < 0]))

write.table(significant_comparison, file=file.path(outpath2, paste('Results_FDR_', FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep="\t")
write.table(female_biased, file=file.path(outpath2, paste('Results_femalebias_FDR_', FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep="\t")
write.table(male_biased, file=file.path(outpath2, paste('Results_malebias_FDR_', FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep="\t")
write.table(unbiased, file=file.path(outpath2, paste('Results_unbias_FDR_', FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep="\t")

rgb_red <- rgb(1,0,0,1/4)
rgb_blue <- rgb(0,0,1,1/4)

p1 <- hist(male_biased$logFC, breaks=seq(-11,11,0.25))
p2 <- hist(female_biased$logFC, breaks=seq(-11,11,0.25))

dev.off()
par(fig=c(0,1,0.1,1))
dev.copy(pdf,file.path(outpath2, paste('Sexualisation_DE_', FDR2use, '_', compname, '_Bonly_sex_bias_tissue_specific.pdf', sep="")), width=8, height=6)
par(fig=c(0,1,0.1,1))
if (nrow(male_biased) > nrow(female_biased)) { # male_more
plot( p1, col=rgb_blue, xlim=c(-9,9), ylim=range(p1$counts)*1.1, main=paste(compname, "DE genes"), xaxt='n', border=rgb_blue, xlab='', ylab="Gene number")
plot( p2, col=rgb_red, xlim=c(-9,9), ylim=range(p2$counts)*1.1, add=T, border=rgb_red)
abline(v=0, lty=2, lw=2)
} else {
plot( p2, col=rgb_red, xlim=c(-9,9),ylim=range(p2$counts)*1.1, main=paste(compname,"DE genes"), xaxt='n', border=rgb_red, xlab='', ylab="Gene number")
plot( p1, col=rgb_blue, xlim=c(-9,9), ylim=range(p1$counts)*1.1, add=T, border=rgb_blue)
abline(v=0, lty=2, lw=3)
} 

par(fig=c(0,1,0,0.40), new=TRUE)
plot(c(-9.5,9.5), c(0.5,3.5), type="n", ylab='', yaxt='n', xlab='logFC', axes=F)
boxplot(female_biased$logFC, unbiased$logFC, male_biased$logFC, col=c(rgb_red, 'white', rgb_blue), xlab="logFC", horizontal=T, boxwex=0.8, add=T, yaxt='n', frame.plot=FALSE)

# var.test(male_biased$logFC, unbiased$logFC)
if ( length(male_biased$logFC) > 0 & length(unbiased$logFC) > 0) {
	male_p <- wilcox.test(male_biased$logFC, unbiased$logFC)$p.value
	} else { 
		male_p <- 1
}
if ( length(female_biased$logFC) > 0 & length(unbiased$logFC) > 0) {
	female_p <- wilcox.test(female_biased$logFC, unbiased$logFC)$p.value
	} else { 
		female_p <- 1
}
if ( length(male_biased$logFC) > 0 & length(female_biased$logFC) > 0) {
	male_female_p <- wilcox.test(male_biased$logFC, female_biased$logFC)$p.value
	} else { 
		male_female_p <- 1
}

if (male_p < 0.001) {text(9.6, 2.5, labels="***", cex=2)} else if (male_p < 0.01) {text(9.6, 2.5, labels="**", cex=2)} else if (male_p < 0.05) {text(9.6, 2.5, labels="*", cex=2)}
if (male_p < 0.05) {text(8.8, 2.5, labels="}", cex=1.2)}

if (female_p < 0.001) {text(9.6, 1.5, labels="***", cex=2)} else if (female_p < 0.01) {text(9.6, 1.5, labels="**", cex=2)} else if (female_p < 0.05) {text(9.6, 1.5, labels="*", cex=2)}
if (female_p < 0.05) {text(8.8, 1.5, labels="}", cex=1.2)}

if (male_female_p < 0.001) {text(-9.4, 2, labels="***", cex=2)} else if (male_female_p < 0.01) {text(-9.4, 2, labels="**", cex=2)} else if (male_female_p < 0.05) {text(-9.4, 2, labels="*", cex=2)}
if (male_female_p < 0.05) {text(-8.4, 2, labels="{", cex=3)}
	
pasteM = paste("♂: +ve", nr_mal_positive, ", -ve", nr_mal_negative, sep=" ")
pasteF = paste("♀: +ve", nr_fem_positive, ", -ve", nr_fem_negative, sep=" ")
pasteN = paste("☯: +ve", nr_neutr_positive, ", -ve", nr_neutr_negative, sep=" ")
mtext(c(pasteM,pasteN,pasteF),side=1,line=2,at=c(-7,0,7), cex=0.9, col=c(rgb(0,0,1,1/2), 'black', rgb(1,0,0,1/2)))
dev.off()



# elina plot
par(mar=c(5,5,4,3))
dev.copy(pdf,file.path(outpath2,paste('Elina_plot_DE_',FDR2use, '_', compname,'.pdf', sep="")), width=8, height=8)
ElinaPlotMergeDE = merge(sexbiasdata, significant_comparison, by.x="geneID", by.y="gene", all=F )
par(fig=c(0,1,0.1,1))
plot(ElinaPlotMergeDE$logFC, ElinaPlotMergeDE[[2]], type="n", xlab="EM logFC", ylab="sex logFC", main=paste(compname, "- DE"), cex.main=1.5, cex.lab=1.2)
points(ElinaPlotMergeDE$logFC[ElinaPlotMergeDE$sexbias=="neutral"], ElinaPlotMergeDE[[2]][ElinaPlotMergeDE$sexbias=="neutral"], pch=16, col=col1)
points(ElinaPlotMergeDE$logFC[ElinaPlotMergeDE$sexbias=="female"], ElinaPlotMergeDE[[2]][ElinaPlotMergeDE$sexbias=="female"], pch=16, col=col2, cex=1)
points(ElinaPlotMergeDE$logFC[ElinaPlotMergeDE$sexbias=="male"], ElinaPlotMergeDE[[2]][ElinaPlotMergeDE$sexbias=="male"], pch=16, col=col4, cex=1)
# legend('topright', inset=0.05, legend=c('Female bias', 'Male bias', 'Unbiased'), pch =c(20,20,1), col=c(2,4,1) )
abline(v=0, lty=2, lw=1)
abline(h=0, lty=2, lw=1)
par(fig=c(0,1,0,0.40), new=TRUE)
pasteM = paste("♂: +ve", nr_mal_positive, ", -ve", nr_mal_negative, sep=" ")
pasteF = paste("♀: +ve", nr_fem_positive, ", -ve", nr_fem_negative, sep=" ")
pasteN = paste("☯: +ve", nr_neutr_positive, ", -ve", nr_neutr_negative, sep=" ")
mtext(c(pasteM,pasteN,pasteF),side=1,line=2,at=c(-4,-0.3,3.5), cex=0.9, col=c(rgb(0,0,1,1/2), 'black', rgb(1,0,0,1/2)))
dev.off()


## Plot all genes
# names(resultscomparison)
female_comparison = merge(femaleDpse, resultscomparison, by.x="geneID", by.y="gene", all=F )
female_available <- nrow(female_comparison)
male_comparison = merge(maleDpse, resultscomparison, by.x="geneID", by.y="gene", all=F )
male_available <- nrow(male_comparison)
unbiased_comparison = merge(neutralDpse, resultscomparison, by.x="geneID", by.y="gene", all=F )
unbiased_available <- nrow(unbiased_comparison)

nr_fem_positive <- nrow(data.frame(female_comparison$logFC[female_comparison$logFC > 0]))
nr_fem_negative <- nrow(data.frame(female_comparison$logFC[female_comparison$logFC < 0]))
nr_mal_positive <- nrow(data.frame(male_comparison$logFC[male_comparison$logFC > 0]))
nr_mal_negative <- nrow(data.frame(male_comparison$logFC[male_comparison$logFC < 0]))
nr_neutr_positive <- nrow(data.frame(unbiased_comparison$logFC[unbiased_comparison$logFC > 0]))
nr_neutr_negative <- nrow(data.frame(unbiased_comparison$logFC[unbiased_comparison$logFC < 0]))

# var.test(male_comparison$logFC, unbiased_comparison$logFC)
male_p <- wilcox.test(male_comparison$logFC, unbiased_comparison$logFC)$p.value
female_p <- wilcox.test(female_comparison$logFC, unbiased_comparison$logFC)$p.value
male_female_p <- wilcox.test(male_comparison$logFC, female_comparison$logFC)$p.value

## plot
p1 <- hist(male_comparison$logFC, breaks=seq(-9.5,9.5,0.25))
p2 <- hist(female_comparison$logFC,  breaks=seq(-9.5,9.5,0.25))

dev.off()
par(fig=c(0,1,0.1,1))
dev.copy(pdf,file.path(outpath2, paste('Sexualisation_all_', compname, '_Bonly_sex_bias_tissue_specific.pdf', sep="")), width=8, height=6)
par(fig=c(0,1,0.1,1))
if (nrow(male_comparison) > nrow(female_comparison)) { # male_more
plot( p1, col=rgb_blue, xlim=c(-4,4), ylim=range(p1$counts)*1.1, main=paste(compname, " - all genes"), xaxt='n', border=rgb_blue, xlab='', ylab="Gene number")
plot( p2, col=rgb_red, xlim=c(-4,4), ylim=range(p2$counts)*1.1, add=T, border=rgb_red)
abline(v=0, lty=2, lw=2)
} else {
par(mar=c(5,5,4,3))
plot( p2, col=rgb_red, xlim=c(-4,4), ylim=range(p2$counts)*1.1, main=paste(compname, " - all genes"), xaxt='n', border=rgb_red, xlab='', ylab="Gene number")
plot( p1, col=rgb_blue, xlim=c(-4,4), ylim=range(p1$counts)*1.1, add=T, border=rgb_blue)
abline(v=0, lty=2, lw=2)
}

par(fig=c(0,1,0,0.40), new=TRUE)
plot(c(-4.5,4.5), c(0.5,3.5), type="n", ylab='', yaxt='n', xlab='logFC', axes=F)
boxplot(female_comparison$logFC, unbiased_comparison$logFC, male_comparison$logFC, col=c(rgb_red, 'white', rgb_blue), xlab="logFC", horizontal=T, boxwex=0.8, add=T, yaxt='n', frame.plot=FALSE)

if (male_p < 0.001) {text(4.4, 2.5, labels="***", cex=1.5)} else if (male_p < 0.01) {text(4.4, 2.5, labels="**", cex=1.5)} else if (male_p < 0.05) {text(4.4, 2.5, labels="*", cex=1.5)}
if (male_p < 0.05) {text(4, 2.5, labels="}", cex=1.2)}

if (female_p < 0.001) {text(4.4, 1.5, labels="***", cex=1.5)} else if (female_p < 0.01) {text(4.4, 1.5, labels="**", cex=1.5)} else if (female_p < 0.05) {text(4.4, 1.5, labels="*", cex=1.5)}
if (female_p < 0.05) {text(4, 1.5, labels="}", cex=1.2)}

if (male_female_p < 0.001) {text(-4.6, 2.5, labels="***", cex=1.5)} else if (male_female_p < 0.01) {text(-4.4, 2, labels="**", cex=1.5)} else if (male_female_p < 0.05) {text(-4.4, 2, labels="*", cex=1.5)}
if (male_female_p < 0.05) {text(-5, 2, labels="{", cex=3)}

pasteM = paste("♂: +ve", nr_mal_positive, ", -ve", nr_mal_negative, '\n', nrow(male_comparison), sep=" ")
pasteF = paste("♀: +ve", nr_fem_positive, ", -ve", nr_fem_negative, '\n', nrow(female_comparison), sep=" ")
pasteN = paste("☯: +ve", nr_neutr_positive, ", -ve", nr_neutr_negative, '\n', nrow(unbiased_comparison), sep=" ")
mtext(c(pasteM,pasteN,pasteF),side=1,line=2,at=c(-3.5,0,3.5), padj=0.5, cex=0.9, col=c(rgb(0,0,1,1/2), 'black', rgb(1,0,0,1/2)))
dev.off()

## Plot all genes - omnibus - flyatlas
tissueDatabase <- c('omnibus', 'flyatlas') 

for(m in 1:length(tissueDatabase)) {
ovary_genes <- read.table(file.path(datapath, paste("ovary_genes_dpse_", tissueDatabase[m], ".txt", sep="")), header=T)
testis_genes <- read.table(file.path(datapath, paste("testis_genes_dpse_", tissueDatabase[m], ".txt", sep="")), header=T)

merged1 <- merge(ovary_genes, resultscomparison, by.x="gene", by.y="gene", all=T )
results_omni_removed <- merged1[is.na(merged1$tissue),]
merged2 <- merge(results_omni_removed, testis_genes, by.x="gene", by.y="gene", all=T)
results_omni_removed2 <- merged2[is.na(merged2$tissue.y),]

female_comparison = merge(femaleDpse, results_omni_removed2, by.x="geneID", by.y="gene", all=F )
# nrow(female_comparison)
male_comparison = merge(maleDpse, results_omni_removed2, by.x="geneID", by.y="gene", all=F )
# nrow(male_comparison)
unbiased_comparison = merge(neutralDpse, results_omni_removed2, by.x="geneID", by.y="gene", all=F )
# nrow(unbiased_comparison)

nr_fem_positive <- nrow(data.frame(female_comparison$logFC[female_comparison$logFC > 0]))
nr_fem_negative <- nrow(data.frame(female_comparison$logFC[female_comparison$logFC < 0]))
nr_mal_positive <- nrow(data.frame(male_comparison$logFC[male_comparison$logFC > 0]))
nr_mal_negative <- nrow(data.frame(male_comparison$logFC[male_comparison$logFC < 0]))
nr_neutr_positive <- nrow(data.frame(unbiased_comparison$logFC[unbiased_comparison$logFC > 0]))
nr_neutr_negative <- nrow(data.frame(unbiased_comparison$logFC[unbiased_comparison$logFC < 0]))

var.test(male_comparison$logFC, unbiased_comparison$logFC)
male_p <- wilcox.test(male_comparison$logFC, unbiased_comparison$logFC)$p.value
female_p <- wilcox.test(female_comparison$logFC, unbiased_comparison$logFC)$p.value
male_female_p <- wilcox.test(male_comparison$logFC, female_comparison$logFC)$p.value

## plot
p1 <- hist(male_comparison$logFC, breaks=seq(-9.5,9.5,0.25))
p2 <- hist(female_comparison$logFC,  breaks=seq(-9.5,9.5,0.25))

dev.off()
par(fig=c(0,1,0.1,1))
dev.copy(pdf,file.path(outpath2, paste('Sexualisation_all_', tissueDatabase[m], '_', compname, '_Bonly_sex_bias_tissue_specific_', '.pdf', sep="")), width=8, height=6)
par(fig=c(0,1,0.1,1))
if (nrow(male_comparison) > nrow(female_comparison) ) { #male_more
plot( p1, col=rgb_blue, xlim=c(-4,4), ylim=range(p1$counts)*1.1, main=paste(compname, " - all genes (",tissueDatabase[m],")", sep="" ), xaxt='n', border=rgb_blue, xlab='', ylab="Gene number")
plot( p2, col=rgb_red, xlim=c(-4,4), ylim=range(p2$counts)*1.1, add=T, border=rgb_red)
abline(v=0, lty=2, lw=2)
} else {
par(mar=c(5,5,4,3))
plot( p2, col=rgb_red, xlim=c(-4,4), ylim=range(p2$counts)*1.1, main=paste(compname, " - all genes (",tissueDatabase[m],")", sep="" ), xaxt='n', border=rgb_red, xlab='', ylab="Gene number")
plot( p1, col=rgb_blue, xlim=c(-4,4), ylim=range(p1$counts)*1.1, add=T, border=rgb_blue)
abline(v=0, lty=2, lw=2)
}

par(fig=c(0,1,0,0.40), new=TRUE)
plot(c(-4.5,4.5), c(0.5,3.5), type="n", ylab='', yaxt='n', xlab='logFC', axes=F)
boxplot(female_comparison$logFC, unbiased_comparison$logFC, male_comparison$logFC, col=c(rgb_red, 'white', rgb_blue), xlab="logFC", horizontal=T, boxwex=0.8, add=T, yaxt='n', frame.plot=FALSE)

if (male_p < 0.001) {text(4.4, 2.5, labels="***", cex=1.5)} else if (male_p < 0.01) {text(4.4, 2.5, labels="**", cex=1.5)} else if (male_p < 0.05) {text(4.4, 2.5, labels="*", cex=1.5)}
if (male_p < 0.05) {text(4, 2.5, labels="}", cex=1.2)}

if (female_p < 0.001) {text(4.4, 1.5, labels="***", cex=1.5)} else if (female_p < 0.01) {text(4.4, 1.5, labels="**", cex=1.5)} else if (female_p < 0.05) {text(4.4, 1.5, labels="*", cex=1.5)}
if (female_p < 0.05) {text(4, 1.5, labels="}", cex=1.2)}

if (male_female_p < 0.001) {text(-4.4, 2, labels="***", cex=1.5)} else if (male_female_p < 0.01) {text(-4.4, 2, labels="**", cex=1.5)} else if (male_female_p < 0.05) {text(-4.4, 2, labels="*", cex=1.5)}
if (male_female_p < 0.05) {text(-4, 2, labels="{", cex=3)}

pasteM = paste("♂: +ve", nr_mal_positive, ", -ve", nr_mal_negative, sep=" ")
pasteF = paste("♀: +ve", nr_fem_positive, ", -ve", nr_fem_negative, sep=" ")
pasteN = paste("☯: +ve", nr_neutr_positive, ", -ve", nr_neutr_negative, sep=" ")
mtext(c(pasteM,pasteN,pasteF),side=1,line=2,at=c(-3.5,0,3.5), cex=0.9, col=c(rgb(0,0,1,1/2), 'black', rgb(1,0,0,1/2)))
dev.off()
}

## Chisq tests {
## Converted the top genes to FBIDs or GeneIDs for Dpse or Dmel in flymine
chi_out <- data.frame(c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0))
colnames(chi_out) <- c('O1', 'E1', 'O1', 'E2', 'statistic', 'p')
row.names(chi_out) <- c('Male vs female All', 'Sex vs unbiased', 'Male vs female Up', 'Male vs female Down')

female_biased_hits <- nrow(female_biased)
male_biased_hits <- nrow(male_biased)
un_biased_hits <- nrow(unbiased)

total_female <- female_available # nrow(femaleDpse)
total_male <- male_available # nrow(maleDpse)
total_neutral <- unbiased_available # nrow(neutralDpse)

total_sex_biased <- total_female + total_male + total_neutral # 16510

total_significant <- nrow(significant_comparison)

exp_female_biased <- round(total_significant * total_female / total_sex_biased)
exp_male_biased <- round(total_significant * total_male / total_sex_biased)
exp_un_biased <- round(total_significant * total_neutral / total_sex_biased)

Obs  <- c(female_biased_hits, male_biased_hits)
Exp <- c(exp_female_biased, exp_male_biased)

RADobsData <- rbind(Obs, Exp)
colnames(RADobsData)<-c('Females', 'Males')

vals <- cbind(RADobsData[1,1], RADobsData[2,1], RADobsData[1,2], RADobsData[2,2])

par(mar=c(5,5,4,3))

dev.copy(pdf,file.path(outpath2, paste('Chisqs_FDR_', FDR2use, '_', compname,'.pdf', sep="")), width=12, height=12)
par(mfrow=c(2,2)) 
if (female_biased_hits > 0 | male_biased_hits > 0) {
	chi_test <- chisq.test(c(female_biased_hits, male_biased_hits),p=c(total_female/(total_female+total_male), total_male/(total_female+total_male)))
	pasteX = paste("chisq =", round(chi_test$p.value, 4), sep=" ")

	mp <- barplot(RADobsData, names=list("Female biased", "Male biased"), col=c(cbgreen,cbblue), main = paste('DE',compname, 'FDR=',FDR2use, '\nMale vs female biased'), beside=T,  cex.main=1.5, cex.lab=1.2, ylab="Gene number")
	if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}
	if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}

	if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
vals <- cbind(RADobsData[1,1], RADobsData[2,1], RADobsData[1,2], RADobsData[2,2])
	names(vals) <- LETTERS[1:4]
	text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
	chi_out$O1[1] <- round(chi_test$observed[1] , digits=2)
	chi_out$O2[1] <- round(chi_test$observed[2] , digits=2)
	chi_out$E1[1] <- round(chi_test$expected[1] , digits=2)
	chi_out$E2[1] <- round(chi_test$expected[2], digits=2)
	chi_out$statistic[1] <- round(chi_test$statistic[1] , digits=2)
	chi_out$p[1] <- round(chi_test$p.value[1] , digits=2)
}

par(mar=c(5,5,4,3))
Obs2 <- c(female_biased_hits + male_biased_hits, un_biased_hits)
Exp2 <- c(exp_female_biased + exp_male_biased, exp_un_biased)
RADobsData2 <- rbind(Obs2, Exp2)
colnames(RADobsData2)<-c('sex_biased', 'un_biased')
vals <- cbind(RADobsData2[1,1], RADobsData2[2,1], RADobsData2[1,2], RADobsData2[2,2])

if (female_biased_hits + male_biased_hits > 0 | un_biased_hits > 0) {
	chi_test <- chisq.test(c(female_biased_hits + male_biased_hits, un_biased_hits),p=c((total_female+total_male)/(total_female+total_male+total_neutral), total_neutral/(total_female+total_male+total_neutral)))
	pasteX = paste("chisq =", round(chi_test$p.value, 3), sep=" ")
	mp <- barplot(RADobsData2, names=list("Sex biased", "Unbiased"), col=c(cbgreen,cbblue),main = paste('DE',compname, 'FDR=',FDR2use, '\nSex biased vs unbiased'), beside=T,  cex.main=1.5, cex.lab=1.2, ylab="Gene number")
	if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}
	if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}

	if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
	vals <- cbind(RADobsData2[1,1], RADobsData2[2,1], RADobsData2[1,2], RADobsData2[2,2])
	names(vals) <- LETTERS[1:4]
	text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
	chi_out$O1[2] <- round(chi_test$observed[1] , digits=2)
	chi_out$O2[2] <- round(chi_test$observed[2] , digits=2)
	chi_out$E1[2] <- round(chi_test$expected[1] , digits=2)
	chi_out$E2[2] <- round(chi_test$expected[2], digits=2)
	chi_out$statistic[2] <- round(chi_test$statistic[1] , digits=2)
	chi_out$p[2] <- round(chi_test$p.value[1] , digits=2)
}

## upregulated genes chisq
female_biased_up = merge(femaleDpse, up_significant_comparison, by.x="geneID", by.y="gene", all=F )
male_biased_up = merge(maleDpse, up_significant_comparison, by.x="geneID", by.y="gene", all=F )
female_biased_up_hits <- nrow(female_biased_up)
male_biased_up_hits <- nrow(male_biased_up)
female_biased_up_hits <- nrow(female_biased_up)
male_biased_up_hits <- nrow(male_biased_up)

total_up_significant <- nrow(up_significant_comparison)
exp_female_biased_up <- round(total_up_significant * total_female / total_sex_biased)
exp_male_biased_up <- round(total_up_significant * total_male / total_sex_biased)

Obs  <- c(female_biased_up_hits, male_biased_up_hits)
Exp <- c(exp_female_biased_up, exp_male_biased_up)

RADobsData <- rbind(Obs, Exp)
colnames(RADobsData)<-c('Females', 'Males')
vals <- cbind(RADobsData[1,1], RADobsData[2,1], RADobsData[1,2], RADobsData[2,2])
par(mar=c(5,5,4,3))

if (female_biased_up_hits > 0 | male_biased_up_hits > 0) {
	chi_test <- chisq.test(c(female_biased_up_hits, male_biased_up_hits),p=c(total_female/(total_female+total_male), total_male/(total_female+total_male)))
	pasteX = paste("chisq =", round(chi_test$p.value, 3), sep=" ")
	mp <- barplot(RADobsData, names=list("Female biased", "Male biased"), col=c(cbgreen,cbblue), beside=T,  cex.main=1.5, cex.lab=1.2, ylab="Gene number", main = paste('DE UP',compname, 'FDR=',FDR2use, '\nMale vs female biased' ))
	if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
	if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}
	if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}

	names(vals) <- LETTERS[1:4]
	text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
	chi_out$O1[3] <- round(chi_test$observed[1] , digits=2)
	chi_out$O2[3] <- round(chi_test$observed[2] , digits=2)
	chi_out$E1[3] <- round(chi_test$expected[1] , digits=2)
	chi_out$E2[3] <- round(chi_test$expected[2], digits=2)
	chi_out$statistic[3] <- round(chi_test$statistic[1] , digits=2)
	chi_out$p[3] <- round(chi_test$p.value[1], digits=2) 
}

## downregulated genes chisq
female_biased_down=merge(femaleDpse, down_significant_comparison, by.x="geneID", by.y="gene", all=F )
male_biased_down=merge(maleDpse, down_significant_comparison, by.x="geneID", by.y="gene", all=F )
female_biased_down_hits <- nrow(female_biased_down)
male_biased_down_hits <- nrow(male_biased_down)
female_biased_down_hits <- nrow(female_biased_down)
male_biased_down_hits <- nrow(male_biased_down)

total_down_significant <- nrow(down_significant_comparison)
exp_female_biased_down <- round(total_down_significant * total_female / total_sex_biased)
exp_male_biased_down <- round(total_down_significant * total_male / total_sex_biased)

Obs  <- c(female_biased_down_hits, male_biased_down_hits)
Exp <- c(exp_female_biased_down, exp_male_biased_down)

RADobsData <- rbind(Obs, Exp)
colnames(RADobsData)<-c('Females', 'Males')
vals <- cbind(RADobsData[1,1], RADobsData[2,1], RADobsData[1,2], RADobsData[2,2])
par(mar=c(5,5,4,3))

if (female_biased_down_hits > 0 | male_biased_down_hits > 0) {
	chi_test <- chisq.test(c(female_biased_down_hits, male_biased_down_hits),p=c(total_female/(total_female+total_male), total_male/(total_female+total_male)))
	pasteX = paste("chisq =", round(chi_test$p.value, 3))

	mp <- barplot(RADobsData, names=list("Female biased", "Male biased"), col=c(cbgreen,cbblue),main = paste('DE DOWN',compname, 'FDR=',FDR2use, '\nMale vs female biased'), beside=T,  cex.main=1.5, cex.lab=1.2, ylab="Gene number")
	if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
	if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}
	if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=c('Obs', 'Exp'), pch =15, col=c(cbgreen,cbblue) )}
	names(vals) <- LETTERS[1:4]
	text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
	chi_out$O1[4] <- round(chi_test$observed[1] , digits=2)
	chi_out$O2[4] <- round(chi_test$observed[2], digits=2) 
	chi_out$E1[4] <- round(chi_test$expected[1], digits=2) 
	chi_out$E2[4] <- round(chi_test$expected[2], digits=2)
	chi_out$statistic[4] <- round(chi_test$statistic[1], digits=2)
	chi_out$p[4] <- round(chi_test$p.value[1] , digits=2)
}
dev.off()
write.table(chi_out, file=file.path(outpath2, paste(compname, "chi_tbl.txt")), quote=F, row.names=T, sep="\t")

## Elina plot all
par(mar=c(5,5,4,3))
dev.copy(pdf,file.path(outpath2,paste('Elina_plot_All_', compname, '.pdf', sep="")), width=6, height=6)
ElinaPlotMergeAll=merge(sexbiasdata, resultscomparison, by.x="geneID", by.y="gene", all=F )
plot(ElinaPlotMergeAll$logFC, ElinaPlotMergeAll[[2]], type="n", xlab="EM logFC", ylab="sex logFC", main=paste(compname, "- All"), cex.main=1.5, cex.lab=1.2)
points(ElinaPlotMergeAll$logFC[ElinaPlotMergeAll$sexbias=="neutral"], ElinaPlotMergeAll[[2]][ElinaPlotMergeAll$sexbias=="neutral"], pch=16, col=col1)
points(ElinaPlotMergeAll$logFC[ElinaPlotMergeAll$sexbias=="female"], ElinaPlotMergeAll[[2]][ElinaPlotMergeAll$sexbias=="female"], pch=16, col=col2, cex=1)
points(ElinaPlotMergeAll$logFC[ElinaPlotMergeAll$sexbias=="male"], ElinaPlotMergeAll[[2]][ElinaPlotMergeAll$sexbias=="male"], pch=16, col=col4, cex=1)
# legend('topright', inset=0.05, legend=c('Female bias', 'Male bias', 'Unbiased'), pch =c(20,20,1), col=c(2,4,1) )
abline(v=0, lty=2, lw=1)
abline(h=0, lty=2, lw=1)
dev.off()

## Load tissue bias data 
tissueBiasSexSpData <- read.table(file.path(datapath, paste('tissue_bias_B_', gender, '.txt', sep="")), header=T)

resultscomparison_all_out <- resultscomparison
merge_comparison <- merge(tissueBiasSexSpData, resultscomparison_mod_logcpm, by.x="geneID", by.y="gene", all=F )
# nrow(merge_comparison)
# str(merge_comparison)

dev.off()
par(mar=c(5,5,4,3))
dev.copy(pdf,file.path(outpath2, paste('contrast_vs_tissue_bias_', compname, '.pdf', sep="")), width=8, height=8)

plot(merge_comparison$logFC, merge_comparison[[2]], type="n", xlab="EM difference", ylab="BH Baseline bias", main=paste(compname, "vs B tissue bias"), cex.main=1.5, cex.lab=1.2)
points(merge_comparison$logFC, merge_comparison[[2]], pch=16, col=col1)
points(merge_comparison$logFC[merge_comparison$FDR < 0.05], merge_comparison[[2]][merge_comparison$FDR < 0.05], pch=16, col=col2, cex=1)

abline(v=0, lty=2, lw=2, col=col4)
abline(h=0, lty=2, lw=2, col=col4)
dev.off()


## COURTSHIP bias

## subset input
Courtship_bias_data <- read.table(file.path(datapath, paste('Courtship_bias_B_', gender, tissue, '.txt', sep="")), header=T)
merge_c_bias <- merge(Courtship_bias_data, sexbiasdata, by.x="geneID", by.y="geneID", all=F )
merge_c_bias_sig=merge(merge_c_bias, resultscomparison_all_out, by.x="geneID", by.y="gene", all=F )

nr_fem_positive_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] > 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="female"]))
nr_fem_negative_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] < 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="female"]))
nr_mal_positive_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] > 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="male"]))
nr_mal_negative_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] < 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="male"]))
nr_neutr_positive_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] > 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="neutral"]))
nr_neutr_negative_cr <- nrow(data.frame(merge_c_bias_sig[,2][merge_c_bias_sig[,2] < 0 & merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="neutral"]))

par(fig=c(0,1,0.1,1))
dev.copy(pdf,file.path(outpath2,paste('Courtship_vs_sexBias_',FDR2use,'_', compname, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(merge_c_bias_sig[,2], merge_c_bias_sig[[4]], type="n", xlab=paste("VC bias (", gender,tissue,")", sep=""), ylab=paste("MF bias (", tissue, ")", sep=""), main=paste(compname, "- Courtship vs sex bias"), cex.main=1.5, cex.lab=1.2)
points(merge_c_bias_sig[,2][merge_c_bias_sig[,3] < FDR2use], merge_c_bias_sig[[4]][merge_c_bias_sig[,3] < FDR2use], pch=16, col=col1)

points(merge_c_bias_sig[,2][merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="male"], merge_c_bias_sig[[4]][merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="male"], pch=16, col=col4, cex=1)	
points(merge_c_bias_sig[,2][merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="female"], merge_c_bias_sig[[4]][merge_c_bias_sig[,3] < FDR2use & merge_c_bias_sig$sexbias=="female"], pch=16, col=col2, cex=1)
abline(v=0, lty=2, lw=0.4, col=1)
abline(h=0, lty=2, lw=0.4, col=1)

par(fig=c(0,1,0,0.40), new=TRUE)
pasteM = paste("Male: +ve", nr_mal_positive_cr, ", -ve", nr_mal_negative_cr, sep=" ")
pasteF = paste("Female: +ve", nr_fem_positive_cr, ", -ve", nr_fem_negative_cr, sep=" ")
pasteN = paste("Neutral: +ve", nr_neutr_positive_cr, ", -ve", nr_neutr_negative_cr, sep=" ")
xrange <- range(merge_c_bias_sig[,2])
mtext(c(pasteM,pasteN,pasteF),side=1,line=2,at=c(-5,0,5), cex=0.9, col=c(rgb(0,0,1,1/2), 'black', rgb(1,0,0,1/2)))
dev.off()
graphics.off()


par(mar=c(5,5,4,3))
dev.copy(pdf,file.path(outpath2, paste('DE_', FDR2use,'_vs_courtship_vs_sexBias_',compname,'.pdf', sep="")), width=8, height=8)
plot(merge_c_bias_sig[,2], merge_c_bias_sig[,4], type="n", xlab=paste("VC bias", gender,tissue), ylab=paste("MF bias", tissue), main=paste(compname, "- DE vs courtship vs sex bias"), cex.main=1.5, cex.lab=1.2)
points(merge_c_bias_sig[,2], merge_c_bias_sig[,4], pch=16, col=col1, cex=0.5)
points(merge_c_bias_sig[,2][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="male"], merge_c_bias_sig[,4][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="male"], pch=16, col=col4, cex=1)
points(merge_c_bias_sig[,2][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="female"], merge_c_bias_sig[,4][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="female"], pch=16, col=col2, cex=1)
points(merge_c_bias_sig[,2][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="neutral"], merge_c_bias_sig[,4][merge_c_bias_sig$FDR < FDR2use & merge_c_bias_sig$sexbias=="neutral"], pch=16, col=col5, cex=1)

abline(v=0, lty=2, lw=0.4, col=1)
abline(h=0, lty=2, lw=0.4, col=1)
dev.off()




## Map  basics
if (!is.numeric(resultscomparisonannotated$start)) {
	resultscomparisonannotated$start <- as.numeric(levels(resultscomparisonannotated$start))[resultscomparisonannotated$start] 
}

annotation$chr <- as.factor(annotation$chr)
annotation$start <- as.numeric(annotation$start)

mapData <- subset(resultscomparisonannotated, resultscomparisonannotated$chr!='NA') 
mapData2 <- subset(mapData, mapData$start!='NA')
nrow(mapData2)


mapLength <- sum(max(mapData2$end[mapData2$chr=='2_'], na.rm=T), max(mapData2$end[mapData2$chr=='3_'], na.rm=T), max(mapData2$end[mapData2$chr=='XL_group1a'], na.rm=T), max(mapData2$end[mapData2$chr=='XL_group1e'], na.rm=T), max(mapData2$end[mapData2$chr=='XL_group3a'], na.rm=T), max(mapData2$end[mapData2$chr=='XL_group3b'], na.rm=T), max(mapData2$end[mapData2$chr=='XR_group3a'], na.rm=T), max(mapData2$end[mapData2$chr=='XR_group5'], na.rm=T), max(mapData2$end[mapData2$chr=='XR_group6'], na.rm=T), max(mapData2$end[mapData2$chr=='XR_group8'], na.rm=T), max(mapData2$end[mapData2$chr=='4_group1'], na.rm=T), max(mapData2$end[mapData2$chr=='4_group2'], na.rm=T), max(mapData2$end[mapData2$chr=='4_group3'], na.rm=T), max(mapData2$end[mapData2$chr=='4_group4'], na.rm=T), max(mapData2$end[mapData2$chr=='4_group5'], na.rm=T))

chr_XL_group1a_end <- max(mapData2$end[mapData2$chr=='XL_group1a'], na.rm=T)
chr_XL_group1e_end <- chr_XL_group1a_end+(max(mapData2$end[mapData2$chr=='XL_group1e'], na.rm=T))
chr_XL_group3a_end <- chr_XL_group1e_end+(max(mapData2$end[mapData2$chr=='XL_group3a'], na.rm=T))
chr_XL_group3b_end <- chr_XL_group3a_end+(max(mapData2$end[mapData2$chr=='XL_group3b'], na.rm=T))

chr_XL_end <- chr_XL_group3b_end 

chr_XR_group3a_end <- max(mapData2$end[mapData2$chr=='XR_group3a'], na.rm=T) + chr_XL_end
chr_XR_group5_end <- chr_XR_group3a_end+(max(mapData2$end[mapData2$chr=='XR_group5'], na.rm=T))
chr_XR_group6_end <- chr_XR_group5_end+(max(mapData2$end[mapData2$chr=='XR_group6'], na.rm=T))
chr_XR_group8_end <- chr_XR_group6_end+(max(mapData2$end[mapData2$chr=='XR_group8'], na.rm=T))
chr_XR_end <- chr_XR_group8_end

chr_2_end <- max(mapData2$end[mapData2$chr=='2_'], na.rm=T) + chr_XR_end
chr_3_end <- max(mapData2$end[mapData2$chr=='3_'], na.rm=T) + chr_2_end

chr_4_group1_end <- max(mapData2$end[mapData2$chr=='4_group1'], na.rm=T) + chr_3_end
chr_4_group2_end <- chr_4_group1_end+(max(mapData2$end[mapData2$chr=='4_group2'], na.rm=T))
chr_4_group3_end <- chr_4_group2_end+(max(mapData2$end[mapData2$chr=='4_group3'], na.rm=T))
chr_4_group4_end <- chr_4_group3_end+(max(mapData2$end[mapData2$chr=='4_group4'], na.rm=T))
chr_4_group5_end <- chr_4_group4_end+(max(mapData2$end[mapData2$chr=='4_group5'], na.rm=T))

list_de <- list()
list_nonde <- list()



for(logFC_use in c(3, 2, 1, 0) ) {

		list_de <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR < FDR2use & abs(resultscomparisonannotated$logFC) > logFC_use)
		list_nonde <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR > FDR2use | abs(resultscomparisonannotated$logFC) < logFC_use)  

# Chisq comparison Chr XL, XR, 4 with others
	for (chrom in c('2_', '3_', '4_', 'XL_', 'XR_', 'X') ) {
		## estimation of loop counts
		de_chrfocus <- sum(list_de$logFC[grepl(chrom, list_de$chr)] != 'na', na.rm=T) 
		de_chrother <- sum(list_de$logFC[!grepl('Unknown*', list_de$chr) & !grepl(chrom, list_de$chr)] != 'na', na.rm=T)
		nonde_chrfocus <- sum(list_nonde$logFC[grepl(chrom, list_nonde$chr)] != 'na', na.rm=T)
		nonde_chrother <- sum(list_nonde$logFC[!grepl('Unknown*', list_nonde$chr) & !grepl(chrom, list_nonde$chr)] !='na')

		up_chrfocus <- sum(list_de$logFC[grepl(chrom, list_de$chr)] > 0, na.rm=T) 
		up_chrother <- sum(list_de$logFC[!grepl('Unknown*', list_de$chr) & !grepl(chrom, list_de$chr)] > 0 , na.rm=T)
		down_chrfocus <- sum(list_de$logFC[grepl(chrom, list_de$chr)] < 0, na.rm=T) 
		down_chrother <- sum(list_de$logFC[!grepl('Unknown*', list_de$chr) & !grepl(chrom, list_de$chr)] < 0 , na.rm=T)

		total_significant <- de_chrfocus + de_chrother
		total_chrfocus <- sum(de_chrfocus + nonde_chrfocus)
		total_chrother <- sum(de_chrother + nonde_chrother)
		exp_chrfocus <- round(total_significant * total_chrfocus/(total_chrfocus+total_chrother))
		exp_chrother <- round(total_significant * total_chrother/(total_chrfocus+total_chrother))
	 
		if (sum(de_chrfocus + de_chrother) > 10 & de_chrfocus !=0 ) { # if 0 error in chisqtest stops the script 
			Obs <- c(de_chrfocus, de_chrother)
			Exp <- c(exp_chrfocus, exp_chrother)
			chi_obs_data <- rbind(Obs, Exp)
			colnames(chi_obs_data) <- c(chrom, 'Chr other')
			chi_test <- chisq.test(c(de_chrfocus, de_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

			pdf(file.path(outpath2, paste('Chisq_Chr_', chrom, 'FDR_', FDR2use, '_logFC_', logFC_use, '_', compname, '.pdf', sep="")), width=16, height=6)
			par(mfrow=c(1,3))
			par(mar=c(5,5,4,3), oma=c(0,0,2,0))
			mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste(compname, 'DE ', '\nFDR=',FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
			if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
			vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
			names(vals) <- LETTERS[1:4]
			text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
			if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {
				legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
			if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {
				legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
		

		if (sum(up_chrfocus + up_chrother) > 10 & up_chrfocus != 0) { 
			total_significant_up <- up_chrfocus + up_chrother
			exp_chrfocus_up <- round(total_significant_up * total_chrfocus/(total_chrfocus+total_chrother))
			exp_chrother_up <- round(total_significant_up * total_chrother/(total_chrfocus+total_chrother))
			Obs_up <- c(up_chrfocus, up_chrother)
			Exp_up <- c(exp_chrfocus_up, exp_chrother_up)
			chi_obs_data_up <- rbind(Obs_up, Exp_up)
			colnames(chi_obs_data_up) <- c(chrom, 'Chr other')
			chi_test_up <- chisq.test(c(up_chrfocus, up_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX_up <- paste('chisq =', round(chi_test_up$p.value, digits=2), sep=" ")
			mp_up <- barplot(chi_obs_data_up, col=c(cbgreen,cbblue), main = paste(compname, 'UP ', '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
			if (chi_test_up$p.value < 0.001) {title(xlab=pasteX_up, font.lab=2)} else {title(xlab=pasteX_up, font.lab=1)}
			vals_up <- cbind(chi_obs_data_up[1,1], chi_obs_data_up[2,1], chi_obs_data_up[1,2], chi_obs_data_up[2,2])
			names(vals_up) <- LETTERS[1:4]
			if (sum(vals_up[1], vals_up[2]) > sum(vals_up[3], vals_up[4])) {
				legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
			if (sum(vals_up[1], vals_up[2]) < sum(vals_up[3], vals_up[4])) {
				legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
			text(mp_up, vals_up, labels = vals_up, pos = 1, col="white", cex=1.5)
		}

		if (sum(down_chrfocus + down_chrother) > 10 & down_chrfocus != 0) { # if 0 error in chisqtest stops the script
			total_significant_down <- down_chrfocus + down_chrother
			exp_chrfocus_down <- round(total_significant_down * total_chrfocus/(total_chrfocus+total_chrother))
			exp_chrother_down <- round(total_significant_down * total_chrother/(total_chrfocus+total_chrother))
			Obs_down <- c(down_chrfocus, down_chrother)
			Exp_down <- c(exp_chrfocus_down, exp_chrother_down)
			chi_obs_data_down <- rbind(Obs_down, Exp_down)
			colnames(chi_obs_data_down) <- c(chrom, 'Chr other')
			chi_test_down <- chisq.test(c(down_chrfocus, down_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX_down <- paste('chisq =', round(chi_test_down$p.value, digits=2), sep=" ")
			mp_down <- barplot(chi_obs_data_down, col=c(cbgreen,cbblue), main = paste(compname, 'DOWN', '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
			if (chi_test_down$p.value < 0.001) {title(xlab=pasteX_down, font.lab=2)} else {title(xlab=pasteX_down, font.lab=1)}
			vals_down <- cbind(chi_obs_data_down[1,1], chi_obs_data_down[2,1], chi_obs_data_down[1,2], chi_obs_data_down[2,2])
			names(vals_down) <- LETTERS[1:4]
			text(mp_down, vals_down, labels = vals_down, pos = 1, col="white", cex=1.5)
			if (sum(vals_down[1], vals_down[2]) > sum(vals_down[3], vals_down[4])) {
				legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
			if (sum(vals_down[1], vals_down[2]) < sum(vals_down[3], vals_down[4])) {
				legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )
			}
		}
		mtext(paste('Chrom', chrom, 'vs rest'), outer = TRUE, cex = 1.5)
		dev.off()
		}
	## estimation of X counts
		de_XL <- sum(list_de$logFC[grepl('XL_', list_de$chr)] != 'na', na.rm=T) 
		nonde_XL <- sum(list_nonde$logFC[grepl('XL_', list_nonde$chr)] != 'na', na.rm=T) 
		de_XR <- sum(list_de$logFC[grepl('XR_', list_de$chr)] != 'na', na.rm=T) 
		nonde_XR <- sum(list_nonde$logFC[grepl('XR_', list_nonde$chr)] != 'na', na.rm=T) 
		de_X <- de_XL + de_XR
		nonde_X <- nonde_XL + nonde_XR

		up_XL <- sum(list_de$logFC[grepl('XL', list_de$chr)] > 0, na.rm=T) 
		up_XR <- sum(list_de$logFC[grepl('XR', list_de$chr)] > 0, na.rm=T) 
		down_XL <- sum(list_de$logFC[grepl('XL', list_de$chr)] < 0, na.rm=T) 
		down_XR <- sum(list_de$logFC[grepl('XR', list_de$chr)] < 0, na.rm=T) 

		if (chrom == 'XL_') {
			total_significant <- de_XL + de_XR
			total_chrfocus <- sum(de_XL + nonde_XL)
			total_chrother <- sum(de_XR + nonde_XR)
			exp_chrfocus <- round(total_significant * total_chrfocus/(total_chrfocus+total_chrother))
			exp_chrother <- round(total_significant * total_chrother/(total_chrfocus+total_chrother))
	 
				if (sum(de_XL + de_XR) > 10 & de_XL !=0 ) { # if 0 error in chisqtest stops the script 
					Obs <- c(de_XL, de_XR)
					Exp <- c(exp_chrfocus, exp_chrother)
					chi_obs_data <- rbind(Obs, Exp)
					colnames(chi_obs_data) <- c(chrom, 'XR')
					chi_test <- chisq.test(c(de_XL, de_XR), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
					pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")
					pdf(file.path(outpath2, paste('Chisq_Chr_', chrom,'vs_XR_FDR_', FDR2use, '_logFC_', logFC_use, '_',compname, '.pdf', sep="")), width=16, height=6)
					par(mfrow=c(1,3))
					par(mar=c(5,5,4,3), oma=c(0,0,2,0))
					mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste(compname, 'DE ', '\nFDR=',FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
					if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
					vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
					names(vals) <- LETTERS[1:4]
				text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
					if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
					if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
				}

				if (sum(up_XL + up_XR) > 10 & up_XL != 0) { 
					total_significant_up <- up_XL + up_XR
					exp_chrfocus_up <- round(total_significant_up * total_chrfocus/(total_chrfocus+total_chrother))
					exp_chrother_up <- round(total_significant_up * total_chrother/(total_chrfocus+total_chrother))
					Obs_up <- c(up_XL, up_XR)
					Exp_up <- c(exp_chrfocus_up, exp_chrother_up)
					chi_obs_data_up <- rbind(Obs_up, Exp_up)
					colnames(chi_obs_data_up) <- c(chrom, 'XR')
					chi_test_up <- chisq.test(c(up_XL, up_XR), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
					pasteX_up <- paste('chisq =', round(chi_test_up$p.value, digits=2), sep=" ")
					mp_up <- barplot(chi_obs_data_up, col=c(cbgreen,cbblue), main = paste(compname, 'UP ', '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
					if (chi_test_up$p.value < 0.001) {title(xlab=pasteX_up, font.lab=2)} else {title(xlab=pasteX_up, font.lab=1)}
					vals_up <- cbind(chi_obs_data_up[1,1], chi_obs_data_up[2,1], chi_obs_data_up[1,2], chi_obs_data_up[2,2])
					names(vals_up) <- LETTERS[1:4]
					if (sum(vals_up[1], vals_up[2]) > sum(vals_up[3], vals_up[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
					if (sum(vals_up[1], vals_up[2]) < sum(vals_up[3], vals_up[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
					text(mp_up, vals_up, labels = vals_up, pos = 1, col="white", cex=1.5)
				}

				if (sum(down_XL + down_XR) > 10 & down_XL != 0) { # if 0 error in chisqtest stops the script
					total_significant_down <- down_XL + down_XR
					exp_chrfocus_down <- round(total_significant_down * total_chrfocus/(total_chrfocus+total_chrother))
					exp_chrother_down <- round(total_significant_down * total_chrother/(total_chrfocus+total_chrother))
					Obs_down <- c(down_XL, down_XR)
					Exp_down <- c(exp_chrfocus_down, exp_chrother_down)
					chi_obs_data_down <- rbind(Obs_down, Exp_down)
					colnames(chi_obs_data_down) <- c(chrom, 'XR')
					chi_test_down <- chisq.test(c(down_XL, down_XR), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
					pasteX_down <- paste('chisq =', round(chi_test_down$p.value, digits=2), sep=" ")
					mp_down <- barplot(chi_obs_data_down, col=c(cbgreen,cbblue), main = paste(compname, 'DOWN', '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
					if (chi_test_down$p.value < 0.001) {title(xlab=pasteX_down, font.lab=2)} else {title(xlab=pasteX_down, font.lab=1)}
					vals_down <- cbind(chi_obs_data_down[1,1], chi_obs_data_down[2,1], chi_obs_data_down[1,2], chi_obs_data_down[2,2])
					names(vals_down) <- LETTERS[1:4]
					text(mp_down, vals_down, labels = vals_down, pos = 1, col="white", cex=1.5)
					if (sum(vals_down[1], vals_down[2]) > sum(vals_down[3], vals_down[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
					if (sum(vals_down[1], vals_down[2]) < sum(vals_down[3], vals_down[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
				}
				if (sum(de_XL + de_XR) > 10 & de_XL !=0 ) {
			mtext(paste('Chrom', chrom, 'vs XR'), outer = TRUE, cex = 1.5)
			dev.off()
			}
		}
	}

## Physical map
		par(mar=c(5,5,4,3))
		pdf(file.path(outpath2, paste('Chr_location_FDR_', FDR2use, '_logFC_', logFC_use, '_',compname, '.pdf', sep="")), width=16, height=8)
		plot(c(0, mapLength), c(min(1-log(resultscomparisonannotated$FDR)), max(1-log(resultscomparisonannotated$FDR))), axes=F, lwd=2, xlab="Chromosome", ylab="1-log(p)", type="null", col="black", main=compname)
		axis(2, c(0, chr_4_group5_end, 5,10,15))

		points(list_nonde$start[list_nonde$chr=='XL_group1a'], 1-log(list_nonde$FDR[list_nonde$chr=='XL_group1a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group1a_end + list_nonde$start[list_nonde$chr=='XL_group1e'], 1-log(list_nonde$FDR[list_nonde$chr=='XL_group1e']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group1e_end + list_nonde$start[list_nonde$chr=='XL_group3a'], 1-log(list_nonde$FDR[list_nonde$chr=='XL_group3a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group3a_end + list_nonde$start[list_nonde$chr=='XL_group3b'], 1-log(list_nonde$FDR[list_nonde$chr=='XL_group3b']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_end + list_nonde$start[list_nonde$chr=='XR_group3a'], 1-log(list_nonde$FDR[list_nonde$chr=='XR_group3a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group3a_end + list_nonde$start[list_nonde$chr=='XR_group5'], 1-log(list_nonde$FDR[list_nonde$chr=='XR_group5']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group5_end + list_nonde$start[list_nonde$chr=='XR_group6'], 1-log(list_nonde$FDR[list_nonde$chr=='XR_group6']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group6_end + list_nonde$start[list_nonde$chr=='XR_group8'], 1-log(list_nonde$FDR[list_nonde$chr=='XR_group8']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_end + list_nonde$start[list_nonde$chr=='2_'], 1-log(list_nonde$FDR[list_nonde$chr=='2_']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_2_end + list_nonde$start[list_nonde$chr=='3_'], 1-log(list_nonde$FDR[list_nonde$chr=='3_']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_3_end + list_nonde$start[list_nonde$chr=='4_group1'], 1-log(list_nonde$FDR[list_nonde$chr=='4_group1']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group1_end + list_nonde$start[list_nonde$chr=='4_group2'], 1-log(list_nonde$FDR[list_nonde$chr=='4_group2']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group2_end + list_nonde$start[list_nonde$chr=='4_group3'], 1-log(list_nonde$FDR[list_nonde$chr=='4_group3']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group3_end + list_nonde$start[list_nonde$chr=='4_group4'], 1-log(list_nonde$FDR[list_nonde$chr=='4_group4']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group4_end + list_nonde$start[list_nonde$chr=='4_group5'], 1-log(list_nonde$FDR[list_nonde$chr=='4_group5']), pch=16, lwd=2, type="p",col=col1, main="")

		abline (v = c(chr_XL_end), lty='dashed')
		abline (v = c(chr_XL_group1a_end), col=col1, lty='dashed')
		abline (v = c(chr_XL_group1e_end), col=col1, lty='dashed')
		abline (v = c(chr_XL_group3a_end), col=col1, lty='dashed')
		abline (v = c(chr_XR_end), lty='dashed')
		abline (v = c(chr_XR_group3a_end), col=col1, lty='dashed')
		abline (v = c(chr_XR_group5_end), col=col1, lty='dashed')
		abline (v = c(chr_XR_group6_end), col=col1, lty='dashed')
		abline (v = c(chr_2_end), lty='dashed')
		abline (v = c(chr_3_end), lty='dashed')
		abline (v = c(chr_4_group1_end), col=col1, lty='dashed')
		abline (v = c(chr_4_group2_end), col=col1, lty='dashed')
		abline (v = c(chr_4_group3_end), col=col1, lty='dashed')
		abline (v = c(chr_4_group4_end), col=col1, lty='dashed')

		chrPosition <- c(chr_XL_group1a_end/2, chr_XL_group1a_end/2 + chr_XL_group1e_end/2, chr_XL_group1e_end/2 + chr_XL_group3a_end/2, chr_XL_group3a_end/2 + chr_XL_group3b_end/2, chr_XL_group3b_end/2 + chr_XR_group3a_end/2, chr_XR_group3a_end/2 + chr_XR_group5_end/2, chr_XR_group5_end/2 + chr_XR_group6_end/2, chr_XR_group6_end/2 + chr_XR_group8_end/2, chr_XR_group8_end/2 + chr_2_end/2, chr_2_end/2 + chr_3_end/2, chr_3_end/2 + chr_4_group1_end/2, chr_4_group1_end/2 + chr_4_group2_end/2, chr_4_group2_end/2 + chr_4_group3_end/2, chr_4_group3_end/2 + chr_4_group4_end/2, chr_4_group4_end/2 + chr_4_group5_end/2) 

		chromosomes <- c("XL_group1a", "XL_group1e", "XL_group3a", "XL_group3b", "XR_group3a", "XR_group5", "XR_group6", "XR_group8", "2", "3", "4_group1", "4_group2", "4_group3", "4_group4", "4_group5")

		axis(side=1, lty=0, at = chrPosition, cex.axis=.5, las=1, labels=chromosomes,  las=2)

		# outliers
		points(list_de$start[list_de$chr=='XL_group1a' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group1a' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="") 
		points(list_de$start[list_de$chr=='XL_group1a' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group1a' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="") 
		points(chr_XL_group1a_end + list_de$start[list_de$chr=='XL_group1e' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group1e' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group1a_end + list_de$start[list_de$chr=='XL_group1e' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group1e' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_group1e_end + list_de$start[list_de$chr=='XL_group3a' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group3a' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group1e_end + list_de$start[list_de$chr=='XL_group3a' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group3a' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_group3a_end + list_de$start[list_de$chr=='XL_group3b' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group3b' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group3a_end + list_de$start[list_de$chr=='XL_group3b' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XL_group3b' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_end + list_de$start[list_de$chr=='XR_group3a' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group3a' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_end + list_de$start[list_de$chr=='XR_group3a' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group3a' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group3a_end + list_de$start[list_de$chr=='XR_group5' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group5' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group3a_end + list_de$start[list_de$chr=='XR_group5' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group5' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group5_end + list_de$start[list_de$chr=='XR_group6' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group6' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group5_end + list_de$start[list_de$chr=='XR_group6' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group6' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group6_end + list_de$start[list_de$chr=='XR_group8' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group8' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group6_end + list_de$start[list_de$chr=='XR_group8' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='XR_group8' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_end + list_de$start[list_de$chr=='2_' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='2_' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_end + list_de$start[list_de$chr=='2_' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='2_' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_2_end + list_de$start[list_de$chr=='3_' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='3_' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_2_end + list_de$start[list_de$chr=='3_' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='3_' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_3_end + list_de$start[list_de$chr=='4_group1' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group1' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_3_end + list_de$start[list_de$chr=='4_group1' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group1' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
	points(chr_4_group1_end + list_de$start[list_de$chr=='4_group2' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group2' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group1_end + list_de$start[list_de$chr=='4_group2' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group2' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group2_end + list_de$start[list_de$chr=='4_group3' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group3' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group2_end + list_de$start[list_de$chr=='4_group3' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group3' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group3_end + list_de$start[list_de$chr=='4_group4' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group4' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group3_end + list_de$start[list_de$chr=='4_group4' & list_de$logFC < -logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group4' & list_de$logFC < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group4_end + list_de$start[list_de$chr=='4_group5' & list_de$logFC > logFC_use], 1-log(list_de$FDR[list_de$chr=='4_group5' & list_de$logFC > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
	legend('topleft', inset=0.05, legend=c('E up', 'M up'), pch =c(25,25), col=c(col2,col4), cex=0.5 )
		dev.off()
}


### Output subsets FDR05

FDR05 <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR < FDR2use)
FDR05_unbiased <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR > FDR2use)
UP <- subset(FDR05, FDR05$logFC > 0)
DOWN <- subset(FDR05, FDR05$logFC < 0)

write.table(list_de$gene, file=file.path(outpath2, paste('gene_',FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep='\t')

write.table(data.frame(gene=FDR05_unbiased$gene, logFC=FDR05_unbiased$logFC, FDR=FDR05_unbiased$FDR), file=file.path(outpath2, paste('unbiased_',FDR2use, '_', compname, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=UP$gene, logFC=UP$logFC, FDR=UP$FDR), file=file.path(outpath2, paste('UP_',FDR2use, '_', compname, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=DOWN$gene, logFC=DOWN$logFC, FDR=DOWN$FDR), file=file.path(outpath2, paste('DOWN_',FDR2use, '_', compname, '.txt', sep="")), quote=F, sep='\t')

gopath <- file.path(outpath2, paste('GO_',FDR2use, '_', compname, sep=""))
dir.create(gopath)

write.table(data.frame(gene=FDR05$gene), file=file.path(gopath,'FDR_genes.txt'), quote=F, row.names=F, sep='\t')
write.table(data.frame(gene=UP$gene, logFC=UP$logFC), file=file.path(gopath, '_E_biased.txt'), quote=F, row.names=F, sep="\t")
write.table(data.frame(gene=DOWN$gene, logFC=DOWN$logFC), file=file.path(gopath, '_M_biased.txt'), quote=F, row.names=F, sep="\t")

GO_pvalues <- data.frame(gene=as.character(resultscomparisonannotated$gene), FDR=resultscomparisonannotated$FDR, logFC=resultscomparisonannotated$logFC)
write.table(GO_pvalues, file=file.path(gopath, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_UP <- data.frame(gene=as.character(resultscomparisonannotated$gene), FDR=resultscomparisonannotated$FDR)
GO_pvalues_UP[,2][resultscomparisonannotated$logFC < 0 & resultscomparisonannotated$FDR < FDR2use] <- 1
write.table(GO_pvalues_UP, file=file.path(gopath, "GO_pvalues_UP.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_DOWN <- data.frame(gene=as.character(resultscomparisonannotated$gene), FDR=resultscomparisonannotated$FDR)
GO_pvalues_DOWN[,2][resultscomparisonannotated$logFC > 0 & resultscomparisonannotated$FDR < FDR2use] <- 1
write.table(GO_pvalues_DOWN, file=file.path(gopath, "GO_pvalues_DOWN.txt"), quote=F, row.names=F, sep="\t")


## Violin plot
par(mar=c(5,5,4,3))
### all genes
violinData <- data.frame(subset(resultscomparisonannotated, !grepl("Unknown*", chr))['chr'], subset(resultscomparisonannotated, !grepl("Unknown*", chr))['logFC'])

violinData$chrCombined <- NULL
violinData$chrCombined[grepl("XR_", violinData$chr)] <- 'XR'
violinData$chrCombined[grepl("XL_", violinData$chr)] <- 'XL'
violinData$chrCombined[grepl("2_", violinData$chr)] <- '2'
violinData$chrCombined[grepl("3_", violinData$chr)] <- '3'
violinData$chrCombined[grepl("4_", violinData$chr)] <- '4'
violinData$chrCombined <- as.factor(violinData$chrCombined)

p = ggplot(violinData, aes_string(factor(violinData$chrCombined), 'logFC'), ylim)
p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC', compname, 'All'), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
ggsave(file.path(outpath2, paste('Violin_per_chr_all_', compname, '.pdf', sep="")), width=8, height=8)

# tests on chromosomal effects of logFC
# library(FSA)
 mod1 <- lm(logFC ~ chrCombined, data=violinData)
 anova(mod1)
## summary(mod1)
## aov(logFC~chrCombined, data=violinData)
 TukeyHSD(aov(logFC~chrCombined, data=violinData))
# bartlett.test(logFC~chrCombined, data=violinData)
# kruskal.test(logFC~chrCombined, data=violinData)
# PT <- dunnTest(logFC~chrCombined, data=violinData)
# PT
# compname
tapply(violinData$logFC, violinData$chrCombined, mean)

# write.table(violinData, file=file.path(outpath, paste("aov_data", compname,".txt", sep="")), quote=F, row.names=F, sep="\t")


### DE genes
restab_logCPM_subset <- subset(resultscomparisonannotated, resultscomparisonannotated$FDR < FDR2use) 
violinData <- data.frame(subset(restab_logCPM_subset, !grepl("Unknown*", chr))['chr'], subset(restab_logCPM_subset, !grepl("Unknown*", chr))['logFC'])
violinData$chrCombined <- NULL
violinData$chrCombined[grepl("XR_", violinData$chr)] <- 'XR'
violinData$chrCombined[grepl("XL_", violinData$chr)] <- 'XL'
violinData$chrCombined[grepl("2_", violinData$chr)] <- '2'
violinData$chrCombined[grepl("3_", violinData$chr)] <- '3'
violinData$chrCombined[grepl("4_", violinData$chr)] <- '4'
violinData$chrCombined <- as.factor(violinData$chrCombined)

p = ggplot(violinData, aes_string(factor(violinData$chrCombined), 'logFC'), ylim)
p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC',compname, 'DE',FDR2use), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
ggsave(file.path(outpath2, paste('Violin_per_chr_DE_',FDR2use, '_', compname, '.pdf', sep="")), width=8, height=8)


## Heatmaps FDR
rownames(resultscomparison_mod_logcpm) <- resultscomparison_mod_logcpm$gene
DE_counts <- subset(resultscomparison_mod_logcpm, resultscomparison_mod_logcpm$FDR < FDR2use)
if (courtship!='All') {
keepcolumns1 <- paste(ssline1, gender,courtship, tissue, sep="_")
keepcolumns2 <- paste(ssline2, gender,courtship, tissue, sep="_")
DE_counts_relevant <- DE_counts[ , grepl( paste(c(keepcolumns1,keepcolumns2),collapse="|") , names( DE_counts ) ) ]
colnames(DE_counts_relevant) <- c('M1','M2','M3','M4','E1','E2','E3','E4')
}
if (courtship=='All') {
	if (gender=='F') {
		DE_counts_relevant <- DE_counts[ , grepl( paste(c("E_F","M_F"),collapse="|") , names( DE_counts ) ) ]
	}
	if (gender=='M') {
		DE_counts_relevant <- DE_counts[ , grepl( paste(c("E_M","M_M"),collapse="|") , names( DE_counts ) ) ]
	}
colnames(DE_counts_relevant) <- c('M1c','M2c','M3c','M4c','M1v','M2v','M3v','M4v','E1c','E2c', 'E3c', 'E4c','E1v','E2v', 'E3v','E4v')
}
DE_colour <- DE_counts$bias
DE_colour <- gsub("female", "2", DE_colour)
DE_colour <- gsub("male", "4", DE_colour)
DE_colour <- gsub("empty", "1", DE_colour)
DE_colour <- gsub("neutral", "1", DE_colour)


if (nrow(DE_counts_relevant) > 2) {
pdf(file.path(outpath2, paste('Heatmap_robustDE_', compname, '_', FDR2use, '.pdf', sep="")), width=8, height=8)
cmethods <- c('median', 'average', 'ward.D') 
for(m in 1:length(cmethods)) {
d <- as.matrix(DE_counts_relevant)
myheatcol <- colorpanel(100, cbred,'white', cbblue) # choose a color palette for the heat map
distmatrix <- as.dist(1-cor(t(d), method="pearson"))
hr <- hclust(distmatrix, method=cmethods[m])  # plot(hr)
mycl <- cutree(hr, k=2) # dynamic tree cut

heatmap.2(d, col=myheatcol,  Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.9, cexRow=0.5, main=paste(compname, cmethods[m], 'FDR', FDR2use), srtCol=45, key.title=NA,colRow=DE_colour) 
}
dev.off()
}



if (nrow(DE_counts_relevant) > 2) {

pdf(file.path(outpath2, paste('Heatmap_DE_', compname, '_', FDR2use, '.pdf', sep="")), width=8, height=8)
# heatmap.2(as.matrix(DE_counts_relevant), col=colorpanel(100, cbred,'white', cbblue), scale="row", key=T, keysize=1.3, density.info="density", trace="none", cexCol=0.5, cexRow=0.6, main=paste(sub_analyse, strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use), srtCol=0, key.title=NA)

cmethods <- c('median', 'average', 'ward.D') # complete ward.D average median ward.D2 mcquitty centroid single

for(m in 1:length(cmethods)) {
d <- as.matrix(DE_counts_relevant)
myheatcol <- colorpanel(100, cbred,'white', cbblue) # choose a color palette for the heat map
distmatrix <- as.dist(1-cor(t(d), method="pearson"))
hr <- hclust(distmatrix, method=cmethods[m])  # plot(hr)
mycl <- cutreeDynamic(hr, distM=as.matrix(distmatrix)) # dynamic tree cut
clusterCols <- rainbow(length(unique(mycl))) # get a color palette equal to the number of clusters
# myClusterSideBar <- clusterCols[mycl] # create vector of colors for side bar, old code
myClusterSideBar <- as.character(as.numeric(mycl))
heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.9, cexRow=0.5, RowSideColors = myClusterSideBar, main=paste(compname, cmethods[m], 'FDR', FDR2use), srtCol=45, key.title=NA,colRow=DE_colour)
# heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.9, cexRow=0.5, main=paste(compname, cmethods[m], 'FDR', FDR2use), srtCol=45, key.title=NA,colRow=DE_colour) # if tree cut does not work

legend(x=0.15,y=1.12, legend = unique(mycl), col = unique(as.numeric(mycl)), lty= 1, lwd = 3, cex=.5, title="clusters", xpd=T)
DE_counts[[ cmethods[m] ]] = as.factor(mycl)
gopath_method <- file.path(outpath2, paste('GO_',FDR2use,'_',cmethods[m], sep=""))
dir.create(gopath_method)

for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
FDR05 <- data.frame(gene=DE_counts$gene, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts$FDR )
FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1
write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')
}
}
dev.off()

print(paste('groups:  median', length(levels(DE_counts$median)), '  ward', length(levels(DE_counts$ward)), '  complete', length(levels(DE_counts$complete)), '  average', length(levels(DE_counts$average)), '  ward2', length(levels(DE_counts$ward.D2)), '  mcquitty', length(levels(DE_counts$mcquitty)), '  centroid', length(levels(DE_counts$centroid)), '  centroid', length(levels(DE_counts$centroid)),sep=" "))

write.table(DE_counts, file=file.path(outpath2, paste('clusters_',FDR2use, '_', compname, '.txt', sep="")), quote=F, row.names=F, sep='\t')

}


