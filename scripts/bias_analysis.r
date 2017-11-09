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

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col1b <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)
col3 <- rgb(red = 0, green = 1, blue = 0, alpha = 0.6)
col4 <- rgb(red = 0, green = 0, blue = 1, alpha = 0.6)
col5 <- rgb(red = 0.7, green = 0, blue = 0.5, alpha = 0.6)

# Base file naming on script input
args <- commandArgs(trailingOnly = TRUE)
sub_analyse = paste(args[1])
FDR2use = as.numeric(paste(args[2]))
bias_analysis = paste(args[3])

# For a manual run
# sub_analyse <- 'B_CH'
# FDR2use <- 0.05
# bias_analysis <- 'sex' # sex courtship tissue

datapath <- "~/git/feminisation_direction/input"
outpath <- paste('~/git/feminisation_direction/output/', bias_analysis, '_bias_analysis', sep="")
dir.create(file.path(outpath))
outpath2 <- paste(outpath, sub_analyse, sep="/")
dir.create(file.path(outpath2))

annotation <- read.delim(file.path(datapath, "exongeneanno.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

count <- read.table(file.path(datapath, paste(sub_analyse,'_count.txt', sep="")), header=T, row.names=1)
design <- read.table(file.path(datapath, paste(sub_analyse,'_design.txt', sep="")), header=T)
model.formula <- as.formula("~0+group")
dmat <- model.matrix(model.formula, data=design)

dgl <- DGEList(counts=count, group=design$group, genes=annotation)

filter_file <- file.path(paste(outpath, '/',sub_analyse, 'filtering_info.txt', sep=""))
ave0 <- dgl[aveLogCPM(dgl) >= 0,]
ave1 <- dgl[aveLogCPM(dgl) >= 1,]
sum2 <- dgl[rowSums(cpm(dgl)) >= 2,]
sum10 <- dgl[rowSums(cpm(dgl)) >= 10,]

write(paste("Comp\tRemaining\tMin\t1Q\tMedian\tMean\t3Q\tMax"), filter_file)
write(paste("All_"), filter_file, append=T)
write(paste(nrow(dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(dgl)/ncol(dgl))), filter_file, append=T, sep='\t', ncol=6)
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

## New to try effect of
# dgl <- dgl[aveLogCPM(dgl) >= 0,]

# dgl <- dgl[rowSums(cpm(dgl)>2) > 5,] # Keep gees expressed in more than half libraries 5 for B, 14 for EMB

dgl <- dgl[aveLogCPM(dgl) >= 0,] # keep average cpm > 0 for further analysis

y <- dgl
colnames(y) <- paste(colnames(y), design$group, sep="\n")

# MDS
pdf(file.path(outpath2,paste('MDS_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plotMDS(y, cex=0.4, col=as.numeric(y$samples$group), main=paste(sub_analyse,"_MDS plot", sep=""))
dev.off()

# estimate data normlization factors and dispersion
xcpm <- mglmOneGroup(dgl$counts) # computing a logCPM for making dispersion plot
dgl <- calcNormFactors(dgl)
dgl <- estimateGLMCommonDisp(dgl, dmat)
dgl <- estimateGLMTrendedDisp(dgl, dmat, min.n=1000)
dgl <- estimateGLMTagwiseDisp(dgl, dmat)

## dispersion plot
pdf(file.path(outpath2, paste('Dispersion_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste(sub_analyse," dispersion", sep=""))
if(!is.null(dgl$trended.dispersion)) points(xcpm ,dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=dgl$common.dispersion, col=cblblue, lwd=2)
legend("topright", c("Common","Trended","Tagwise"), pch=16, col=c(cblblue, cbgreen, "black"), title="Dispersion")
dev.off()


##  fit the data model ------------------------------------------
fitres <- glmFit(dgl, dmat)
x <- read.delim(file.path(datapath, paste(sub_analyse,'_matrix.txt', sep="")), header=T, sep="\t")
sortedX <- data.frame(x[order(x$model_coefficients, decreasing=F),])
cmat <- as.matrix(sortedX[,-1])
colnames(cmat)[1] <- colnames(sortedX[2])
rownames(cmat) <- as.character(sortedX[,1])
# cmat # check cmat loaded correctly

######## ------- obtain contrast fit and test results
lrtres <- list()
for(k in 1:ncol(cmat)) lrtres[[k]] <- glmLRT(fitres, contrast=cmat[,k])

logFC <- NULL
PV <- NULL
FDR <- NULL
for(k in 1:ncol(cmat)) {PV <- cbind(PV, lrtres[[k]]$table[,"PValue"])
	FDR= cbind(FDR,p.adjust(PV[,k],method="BH"))
	logFC= cbind(logFC,lrtres[[k]]$table[,"logFC"])
	}

xcpm <- lrtres[[1]]$table[,"logCPM"]
allzeros <- which(rowSums(cpm(dgl)) < 3,)
allused <- which(rowSums(cpm(dgl)) >= 3,)

cname <-colnames(cmat)
colnames(logFC) <- paste("logFC",cname,sep=".")
colnames(PV) <- paste("PV",cname,sep=".")
colnames(FDR) <- paste("FDR",cname,sep=".")
rownames(logFC) <- rownames(PV) <- rownames(FDR) <- rownames(fitres$coefficients)




## Make results table
idxzeros <- allzeros
restab <- data.frame(rownames(dgl$counts),dgl$genes[,6], dgl$genes[,2], dgl$genes[,3],dgl$genes[,4],logCPM=xcpm,logFC,PV,FDR)
colnames(restab)[1] <- 'gene'
colnames(restab)[2] <- 'gname'
colnames(restab)[3] <- 'chr'
colnames(restab)[4] <- 'start'
colnames(restab)[5] <- 'end'

## pseudo counts export
d = estimateCommonDisp(dgl, verbose=TRUE)
pseu_counts_dgl <- d$pseudo.counts
pseu_counts_dgl <- data.frame(row.names(pseu_counts_dgl), pseu_counts_dgl)
colnames(pseu_counts_dgl) <- c('ID',as.character(y$samples$group))
pseu_counts_dgl2 = merge(restab, pseu_counts_dgl, by.x="gene", by.y="ID", all=F )
pseu_counts_dgl2$bias <- 'empty'

# pseu_counts_dgl2$bias[match(pseu_counts_dgl2$gene, femaleDpse$geneID) !='NA'] <- 'female'
# pseu_counts_dgl2$bias[match(pseu_counts_dgl2$gene, maleDpse$geneID) !='NA'] <- 'male'
# pseu_counts_dgl2$bias[match(pseu_counts_dgl2$gene, neutralDpse$geneID) !='NA'] <- 'neutral'
# pseu_counts_dgl2$bias <- as.factor(pseu_counts_dgl2$bias)

write.table(pseu_counts_dgl2, file=file.path(outpath2, paste(sub_analyse,'_pseudo_counts.txt', sep="")), quote=F, row.names=F, sep="\t")



## pvalue histogram
pdf(file.path(outpath2, paste('p_hist_',FDR2use, '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
npanel <- ncol(logFC)
np <- ceiling(sqrt(npanel))
if(np*(np-1)>= npanel) mfcol <- c(np-1,np) else
{if((np-1)*(np+1)>= npanel) mfcol <- c(np+1,np-1) else mfcol <- c(np,np)}
par(mfcol=mfcol)
for(k in 1:npanel) {
hist(PV[,k],n=100,xlab="P-value",main=colnames(logFC)[k])
}
dev.off()

#----- whinin group pairwise scatter plot
wx <- dgl$counts
wg <- as.character(dgl$samples$group)
wug <- unique(wg)
wn <- length(wug)
for (k in 1:wn) {
	ix <- wg %in% wug[k]
	xmat <- log2(wx[,ix])
	pdf(file.path(outpath2, paste('pairwise_raw_count_', wug[k], '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
if (sum(ix) > 1 ) { 
	pairs(xmat,pch=16,cex=0.4,main=wug[k]) 
	}
#	dev.copy(pdf,file.path(outpath2,'courtship_out',paste(sub_analyse,wug[k],'_pairwise_raw_count.pdf', sep="")), width=8, height=8)
	dev.off()
}	

# write.table(restab, file=file.path(outpath2,paste(sub_analyse,'_all.txt', sep="")), quote=F, row.names=F, sep='\t')

restab_frame <- as.data.frame(restab)

for(logFC_use in c(3, 2, 1, 0) ) {
de.yes.no <- FDR < FDR2use & abs(logFC) > logFC_use
if (ncol(cmat) == 1) {
	de4 <- which((FDR[,1] < FDR2use) == 0 & abs(logFC[,1] > logFC_use) != 0) } else {
de4 <- which(rowSums(FDR[,c(1:ncol(FDR))] < FDR2use) == 0 & rowSums(abs(logFC[,c(1:ncol(FDR))]) > logFC_use) != 0 )
}

de.yes.no[de4,] <- FALSE
deidx <- ii <- rowSums(de.yes.no) > 0
delabel <- (sign(logFC)*de.yes.no)[ii,]
delabel[is.na(delabel)] <- 0
combinedp <- 1; for(k in 1:ncol(PV)) combinedp <- combinedp*PV[,k]
deidx <- ii

wtable_1 <- restab[deidx,] # this is the result table of DE genes.
demat <- as.matrix(logFC[deidx,])
write.table(wtable_1, file=file.path(outpath, paste('de_',FDR2use, '_',logFC_use, '_', sub_analyse,'.txt', sep="")), quote=F, row.names=F, sep='\t')

if (nrow(wtable_1) > 1) {
	if (ncol(cmat) == 1) {
		wtable_3 <- rbind(NUM_DE=sum(abs(delabel)), NUM_UP_DE =sum(delabel>0), NUM_DOWN_DE =sum(delabel<0)) } else {
		wtable_3 <- rbind(NUM_DE=colSums(abs(delabel)), NUM_UP_DE =colSums(delabel>0), NUM_DOWN_DE =colSums(delabel<0))
	}
	wtable_3 <- cbind(DE_numbers=rownames(wtable_3), wtable_3)
	write.table(wtable_3, file=file.path(outpath, paste('Number_de_',FDR2use, '_',logFC_use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
	}
}


fc <- logFC
fc[(fc)>10] <- 10
fc[ fc < -10] <- -10
par(mfcol=mfcol)

for(k in 1:ncol(cmat)) {
	pdf(file.path(outpath2, paste('FC-CPMplot_',FDR2use, '_', sub_analyse,'_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=8, height=8)
	ylab <- colnames(logFC)[k]
	deix <- which(de.yes.no[,k])
	maPlot(x=NULL,y=NULL,logAbundance= xcpm, logFC = fc[,k], xlab = bquote(paste(log^2, CPM)), ylab = paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' - ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), de.tags= deix,pch = 19, cex = 0.3, smearWidth = 0.5, panel.first = grid(), smooth.scatter = FALSE, lowess = FALSE, na.rm =TRUE, main = paste('logFC plot ', strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""))
	dev.off()
}
		
# generation of subsets
# for(k in 2:ncol(restab)) restab[,k] <- as.numeric(levels(restab[,k]))[restab[,k]]
# str(restab)

if (!is.numeric(restab$start)) {
	restab$start <- as.numeric(levels(restab$start))[restab$start] 
}

annotation$chr <- as.factor(annotation$chr)
annotation$start <- as.numeric(annotation$start)

mapData <- subset(restab, restab$chr!='NA') 
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

	for(k in 1:ncol(cmat)) {
		list_de[[k]] <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use & abs(restab[[paste('logFC.',colnames(cmat)[k], sep="")]]) > logFC_use)
		list_nonde[[k]] <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] > FDR2use | abs(restab[[paste('logFC.',colnames(cmat)[k], sep="")]]) < logFC_use)  
	}

# Chisq comparison Chr XL, XR, 4 with others
	for (chrom in c('2_', '3_', '4_', 'XL_', 'XR_', 'X') ) {
		for(k in 1:ncol(cmat)) {

		## estimation of loop counts
		de_chrfocus <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl(chrom, list_de[[k]]$chr)] != 'na', na.rm=T) 
		de_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('Unknown*', list_de[[k]]$chr) & !grepl(chrom, list_de[[k]]$chr)] != 'na', na.rm=T)
		nonde_chrfocus <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl(chrom, list_nonde[[k]]$chr)] != 'na', na.rm=T)
		nonde_chrother <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('Unknown*', list_nonde[[k]]$chr) & !grepl(chrom, list_nonde[[k]]$chr)] !='na')

		up_chrfocus <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl(chrom, list_de[[k]]$chr)] > 0, na.rm=T) 
		up_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('Unknown*', list_de[[k]]$chr) & !grepl(chrom, list_de[[k]]$chr)] > 0 , na.rm=T)
		down_chrfocus <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl(chrom, list_de[[k]]$chr)] < 0, na.rm=T) 
		down_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('Unknown*', list_de[[k]]$chr) & !grepl(chrom, list_de[[k]]$chr)] < 0 , na.rm=T)

		total_significant <- de_chrfocus + de_chrother
		total_chrfocus <- sum(de_chrfocus + nonde_chrfocus)
		total_chrother <- sum(de_chrother + nonde_chrother)
		exp_chrfocus <- round(total_significant * total_chrfocus/(total_chrfocus+total_chrother))
		exp_chrother <- round(total_significant * total_chrother/(total_chrfocus+total_chrother))
	 
		if (sum(de_chrfocus + de_chrother) > 10 & de_chrfocus !=0 ) { # if 0 error in chisqtest stops the script 
			Obs <- c(de_chrfocus, de_chrother)
			Exp <- c(exp_chrfocus, exp_chrother)
			chi_obs_data <- rbind(Obs, Exp)
			colnames(chi_obs_data) <- c('Chr 2', 'Chr other')
			chi_test <- chisq.test(c(de_chrfocus, de_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

			pdf(file.path(outpath2, paste('Chisq_Chr_', chrom, 'FDR_', FDR2use, '_logFC_', logFC_use, '_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=16, height=6)
			par(mfrow=c(1,3))
			par(mar=c(5,5,4,3), oma=c(0,0,2,0))
			mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste('DE ',(colnames(cmat)[k]), '\nFDR=',FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
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
			colnames(chi_obs_data_up) <- c('Chr 2', 'Chr other')
			chi_test_up <- chisq.test(c(up_chrfocus, up_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX_up <- paste('chisq =', round(chi_test_up$p.value, digits=2), sep=" ")
			mp_up <- barplot(chi_obs_data_up, col=c(cbgreen,cbblue), main = paste('Up ',colnames(cmat)[k], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
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
			colnames(chi_obs_data_down) <- c('chrfocus', 'Chr other')
			chi_test_down <- chisq.test(c(down_chrfocus, down_chrother), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
			pasteX_down <- paste('chisq =', round(chi_test_down$p.value, digits=2), sep=" ")
			mp_down <- barplot(chi_obs_data_down, col=c(cbgreen,cbblue), main = paste('Down ',colnames(cmat)[k], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
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
		de_XL <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XL_', list_de[[k]]$chr)] != 'na', na.rm=T) 
		nonde_XL <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XL_', list_nonde[[k]]$chr)] != 'na', na.rm=T) 
		de_XR <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XR_', list_de[[k]]$chr)] != 'na', na.rm=T) 
		nonde_XR <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XR_', list_nonde[[k]]$chr)] != 'na', na.rm=T) 
		de_X <- de_XL + de_XR
		nonde_X <- nonde_XL + nonde_XR

		up_XL <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XL_', list_de[[k]]$chr)] > 0, na.rm=T) 
		up_XR <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XR_', list_de[[k]]$chr)] > 0, na.rm=T) 
		down_XL <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XL_', list_de[[k]]$chr)] < 0, na.rm=T) 
		down_XR <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][grepl('XR_', list_de[[k]]$chr)] < 0, na.rm=T) 

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
					pdf(file.path(outpath2, paste('Chisq_Chr_', chrom,'vs_XR_FDR_', FDR2use, '_logFC_', logFC_use, '_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=16, height=6)
					par(mfrow=c(1,3))
					par(mar=c(5,5,4,3), oma=c(0,0,2,0))
					mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste('DE ',(colnames(cmat)[k]), '\nFDR=',FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
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
					Obs_up <- c(up_chrfocus, up_chrother)
					Exp_up <- c(exp_chrfocus_up, exp_chrother_up)
					chi_obs_data_up <- rbind(Obs_up, Exp_up)
					colnames(chi_obs_data_up) <- c(chrom, 'XR')
					chi_test_up <- chisq.test(c(up_XL, up_XR), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
					pasteX_up <- paste('chisq =', round(chi_test_up$p.value, digits=2), sep=" ")
					mp_up <- barplot(chi_obs_data_up, col=c(cbgreen,cbblue), main = paste('Up ',colnames(cmat)[k], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
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
					Obs_down <- c(down_chrfocus, down_chrother)
					Exp_down <- c(exp_chrfocus_down, exp_chrother_down)
					chi_obs_data_down <- rbind(Obs_down, Exp_down)
					colnames(chi_obs_data_down) <- c(chrom, 'XR')
					chi_test_down <- chisq.test(c(down_XL, down_XR), p=c(total_chrfocus/(total_chrfocus+total_chrother), total_chrother/(total_chrfocus+total_chrother)))
					pasteX_down <- paste('chisq =', round(chi_test_down$p.value, digits=2), sep=" ")
					mp_down <- barplot(chi_obs_data_down, col=c(cbgreen,cbblue), main = paste('Down ',colnames(cmat)[k], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.2, cex.lab=1, ylab="Gene number")
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
	}

## Physical map
	for(k in 1:ncol(cmat)) {
		par(mar=c(5,5,4,3))
		pdf(file.path(outpath2, paste('Chr_location_FDR_', FDR2use, '_logFC_', logFC_use, '_',colnames(cmat)[k], '.pdf', sep="")), width=16, height=8)
		plot(c(0, mapLength), c(min(1-log(restab[[paste('FDR.',colnames(cmat)[k], sep="")]])), max(1-log(restab[[paste('FDR.',colnames(cmat)[k], sep="")]]))), axes=F, lwd=2, xlab="Chromosome", ylab="1-log(p)", type="null", col="black", main=paste(colnames(cmat)[k]))
		axis(2, c(0, chr_4_group5_end, 5,10,15))

		points(list_nonde[[k]]$start[list_nonde[[k]]$chr=='XL_group1a'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XL_group1a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group1a_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XL_group1e'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XL_group1e']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group1e_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XL_group3a'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XL_group3a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_group3a_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XL_group3b'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XL_group3b']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XL_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XR_group3a'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XR_group3a']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group3a_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XR_group5'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XR_group5']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group5_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XR_group6'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XR_group6']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_group6_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='XR_group8'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='XR_group8']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_XR_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='2_'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='2_']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_2_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='3_'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='3_']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_3_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='4_group1'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='4_group1']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group1_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='4_group2'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='4_group2']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group2_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='4_group3'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='4_group3']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group3_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='4_group4'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='4_group4']), pch=16, lwd=2, type="p",col=col1, main="")
		points(chr_4_group4_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='4_group5'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='4_group5']), pch=16, lwd=2, type="p",col=col1, main="")

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
		points(list_de[[k]]$start[list_de[[k]]$chr=='XL_group1a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group1a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="") 
				points(list_de[[k]]$start[list_de[[k]]$chr=='XL_group1a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group1a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="") 
		points(chr_XL_group1a_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group1e' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group1e' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group1a_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group1e' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group1e' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_group1e_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group1e_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_group3a_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group3b' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group3b' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_group3a_end + list_de[[k]]$start[list_de[[k]]$chr=='XL_group3b' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XL_group3b' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XL_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XL_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group3a' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group3a_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group3a_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group5_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group6' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group6' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group5_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group6' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group6' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_group6_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group8' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group8' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_group6_end + list_de[[k]]$start[list_de[[k]]$chr=='XR_group8' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='XR_group8' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_XR_end + list_de[[k]]$start[list_de[[k]]$chr=='2_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='2_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_XR_end + list_de[[k]]$start[list_de[[k]]$chr=='2_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='2_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_2_end + list_de[[k]]$start[list_de[[k]]$chr=='3_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='3_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_2_end + list_de[[k]]$start[list_de[[k]]$chr=='3_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='3_' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_3_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group1' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group1' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_3_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group1' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group1' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
	points(chr_4_group1_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group2' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group2' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group1_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group2' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group2' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group2_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group3' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group3' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group2_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group3' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group3' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group3_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group4' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group4' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		points(chr_4_group3_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group4' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group4' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, main="")
		points(chr_4_group4_end + list_de[[k]]$start[list_de[[k]]$chr=='4_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='4_group5' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, main="")
		legend('topleft', inset=0.05, legend=c('E up', 'M up'), pch =c(25,25), col=c(col2,col4), cex=0.5 )
		dev.off()
		}
}


### Output subsets FDR05
for(k in 1:ncol(cmat)) {
	FDR05 <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
	FDR05_unbiased <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] > FDR2use)
	UP <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0)
	DOWN <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0)

	write.table(list_de[[k]]$gene, file=file.path(outpath2, paste('gene_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

	write.table(data.frame(gene=FDR05_unbiased$gene, logFC=FDR05_unbiased[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=FDR05_unbiased[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath2, paste('unbiased_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
	write.table(data.frame(gene=UP$gene, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=UP[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath2, paste('UP_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
	write.table(data.frame(gene=DOWN$gene, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=DOWN[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath2, paste('DOWN_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')

	gopath <- file.path(outpath2, paste('GO_',FDR2use, '_', colnames(cmat)[k], sep=""))
	dir.create(gopath)

	write.table(data.frame(gene=FDR05$gene), file=file.path(gopath,'FDR_genes.txt'), quote=F, row.names=F, sep='\t')
	write.table(data.frame(gene=UP$gene, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")
	write.table(data.frame(gene=DOWN$gene, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][2], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")

	GO_pvalues <- data.frame(gene=as.character(restab$gene), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]], logFC=restab[[paste('logFC.',colnames(cmat)[k], sep="")]])
	write.table(GO_pvalues, file=file.path(gopath, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

	GO_pvalues_UP <- data.frame(gene=as.character(restab$gene), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
	GO_pvalues_UP[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
	write.table(GO_pvalues_UP, file=file.path(gopath, "GO_pvalues_UP.txt"), quote=F, row.names=F, sep="\t")

	GO_pvalues_DOWN <- data.frame(gene=as.character(restab$gene), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
	GO_pvalues_DOWN[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
	write.table(GO_pvalues_DOWN, file=file.path(gopath, "GO_pvalues_DOWN.txt"), quote=F, row.names=F, sep="\t")
}

## export results and normalised counts
nc <- cpm(dgl, prior.count=2, log=T)
nc2 <- data.frame(row.names(nc), nc)
colnames(nc2) <- c('ID', paste(design$group, design$line, sep="_" ))
restab_logCPM = merge(restab, nc2, by.x="gene", by.y="ID", all=F )
write.table(restab_logCPM, file=file.path(outpath2, paste('LogCPM_',FDR2use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

## EXPORT bias data for use in E-M contrasts
### Courtship: Decided to use courtship bias as calculated from Baseline data only
if (bias_analysis == 'courtship') {
courtship_export <- data.frame(restab_logCPM$gene, restab_logCPM[paste('logFC.S_V.C_', sub_analyse, '_B', sep="")], restab_logCPM[paste('FDR.S_V.C_', sub_analyse, '_B', sep="")])
colnames(courtship_export) <- c('geneID', paste('logFC.', sub_analyse, '_B', sep=""), paste('FDR.', sub_analyse, '_B', sep=""))
write.table(courtship_export, file=file.path(datapath, paste('Courtship_bias_' ,sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
}

### sex bias: used from Baseline only, E and M libraries are not involved in calculating dispersion
if (bias_analysis == 'sex' ) {
sexbias_export <- data.frame(restab_logCPM$gene, restab_logCPM[paste('logFC.G_M.F_', sub_analyse, '_B', sep="")], restab_logCPM[paste('FDR.G_M.F_', sub_analyse, '_B', sep="")])
colnames(sexbias_export) <- c('geneID', paste('logFC.', sub_analyse, '_B', sep=""), paste('FDR.', sub_analyse, '_B', sep=""))
write.table(sexbias_export, file=file.path(datapath, paste('sex_bias_' ,sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
}

### tissue bias: used from baseline only, but with dispersion from E and M libraries 
if (bias_analysis == 'tissue' ) {
tissuebias_export <- data.frame(restab_logCPM$gene, restab_logCPM[paste('logFC.T_H.B_', sub_analyse, '_B', sep="")], restab_logCPM[paste('FDR.T_H.B_', sub_analyse, '_B', sep="")])
colnames(tissuebias_export) <- c('geneID', paste('logFC.', sub_analyse, '_B', sep=""), paste('FDR.', sub_analyse, '_B', sep=""))
write.table(tissuebias_export, file=file.path(datapath, paste('tissue_bias_' ,sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
}

## violin plots
for(k in 1:ncol(cmat)) {
	par(mar=c(5,5,4,3))

	### all genes
	violinData <- data.frame(subset(restab_logCPM, !grepl("Unknown*", chr))["chr"], subset(restab_logCPM, !grepl("Unknown*", chr))[paste('logFC.',colnames(cmat)[k], sep="")])
	violinData$chrCombined <- NULL
	violinData$chrCombined[grepl("XR_", violinData$chr)] <- 'XR'
	violinData$chrCombined[grepl("XL_", violinData$chr)] <- 'XL'
	violinData$chrCombined[grepl("2_", violinData$chr)] <- '2'
	violinData$chrCombined[grepl("3_", violinData$chr)] <- '3'
	violinData$chrCombined[grepl("4_", violinData$chr)] <- '4'
	violinData$chrCombined <- as.factor(violinData$chrCombined)

	p = ggplot(violinData, aes_string(factor(violinData$chrCombined), paste('logFC.', colnames(cmat)[k], sep="")), ylim)
	p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC.',colnames(cmat)[k], sep=""), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
	ggsave(file.path(outpath2, paste('Violin_per_chr_all_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)

	### DE genes
	restab_logCPM_subset <- subset(restab_logCPM, restab_logCPM[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use) 
	violinData <- data.frame(subset(restab_logCPM_subset, !grepl("Unknown*", chr))["chr"], subset(restab_logCPM_subset, !grepl("Unknown*", chr))[paste('logFC.',colnames(cmat)[k], sep="")])
	violinData$chrCombined <- NULL
	violinData$chrCombined[grepl("XR_", violinData$chr)] <- 'XR'
	violinData$chrCombined[grepl("XL_", violinData$chr)] <- 'XL'
	violinData$chrCombined[grepl("2_", violinData$chr)] <- '2'
	violinData$chrCombined[grepl("3_", violinData$chr)] <- '3'
	violinData$chrCombined[grepl("4_", violinData$chr)] <- '4'
	violinData$chrCombined <- as.factor(violinData$chrCombined)

	p = ggplot(violinData, aes_string(factor(violinData$chrCombined), paste('logFC.', colnames(cmat)[k], sep="")), ylim)
	p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC.',colnames(cmat)[k], sep=""), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
	ggsave(file.path(outpath2, paste('Violin_per_chr_DE_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
}


# heatmaps
for(k in 1:ncol(cmat)) {
	if (ncol(cmat) == 1) {
		cmat_subset <- cmat } else {
		cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]]!=0) 
		}
	design_subset <- design[design$group %in% row.names(cmat_subset),]
	design_subset$sample <- paste(design_subset$group, design_subset$line, sep="_")
	length(as.character(design_subset$sample))
	rownames(restab_logCPM) <- restab_logCPM$gene
	DE_counts <- subset(restab_logCPM, restab_logCPM[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
	DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))
	# row.names(DE_counts_relevant) <- DE_counts$gname
	# colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group) 
	if (nrow(DE_counts_relevant) > 2) {

	pdf(file.path(outpath2, paste('Heatmap_', FDR2use, '_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
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
	heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.6, cexRow=0.6, RowSideColors = myClusterSideBar, main=paste(sub_analyse, cmethods[m], colnames(cmat)[k], FDR2use), srtCol=45, key.title=NA)
	legend(x=-1,y=1.12, legend = unique(mycl), col = unique(as.numeric(mycl)), lty= 1, lwd = 3, cex=.5, title="clusters", xpd=T)
	DE_counts[[ cmethods[m] ]] = as.factor(mycl)
	gopath_method <- file.path(outpath2, paste('GO_',FDR2use, '_', colnames(cmat)[k],'/',cmethods[m], sep=""))
	dir.create(gopath_method)

	for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
	FDR05 <- data.frame(gene=DE_counts$gene, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts[[paste('FDR.',colnames(cmat)[k], sep="")]] )
	FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1
	write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')
	}
	}
	dev.off()

	print(paste(colnames(cmat)[k], 'groups:  median', length(levels(DE_counts$median)), '  ward', length(levels(DE_counts$ward)), '  complete', length(levels(DE_counts$complete)), '  average', length(levels(DE_counts$average)), '  ward2', length(levels(DE_counts$ward.D2)), '  mcquitty', length(levels(DE_counts$mcquitty)), '  centroid', length(levels(DE_counts$centroid)), '  centroid', length(levels(DE_counts$centroid)),sep=" "))
	write.table(DE_counts, file=file.path(outpath2, paste('clusters',FDR2use, '_', sub_analyse, '_',colnames(cmat)[k], '.txt', sep="")), quote=F, row.names=F, sep='\t')
	}
}
