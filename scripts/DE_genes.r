require(limma)
require(lmodel2)

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

# args <- commandArgs(trailingOnly = TRUE)
# FDR2use = as.numeric(paste(args[1]))

FDR2use <- 0.1
plottype <- 'sex' # 'EM' 'sex'

annotation <- read.delim(file.path("~/git/feminisation_direction/input/exongeneanno.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

dpse_dmel_data <- read.table('/Users/zabameos/git/feminisation_direction/input/dpseDmel.txt', header=T)
str(dpse_dmel_data)
gerrard_data <- read.table("~/git/feminisation_direction/input/gerrard_genes.txt", header=T)
str(gerrard_data)

if (plottype == 'EM') {
datapath <- "~/git/feminisation_direction/output/subsets/genes"
outpath <- paste("~/git/feminisation_direction/output/subsets/", sep="")
}

if (plottype == 'sex') {
datapath <- "/Users/zabameos/git/feminisation_direction/output/sex_bias_analysis/genes"
outpath <- paste("~/git/feminisation_direction/output/sex_bias_analysis/", sep="")
}

list_results <- list()
list_de <- list()
list_de_up <- list()
list_de_down <- list()
list_de_sex <- list()

## manual step
for(contrast in c('B_VB', 'B_CB', 'B_VH', 'B_CH')) { # sex bias
# for(contrast in c('EMFVB', 'EMMVB', 'EMFVH', 'EMMVH', 'EMFCB', 'EMMCB', 'EMFCH', 'EMMCH')) { # EM bias

	list_results[[contrast]] <- read.table(file.path(datapath, paste(contrast, '_pseudo_counts.txt', sep="")), header=T)
	list_de[[contrast]] <- subset(list_results[[contrast]], list_results[[contrast]]$FDR < FDR2use) 
	list_de_up[[contrast]] <- subset(list_results[[contrast]], list_results[[contrast]]$FDR < FDR2use & list_results[[contrast]][7] > 0 & abs(list_results[[contrast]][7]) > 1) 
	list_de_down[[contrast]] <- subset(list_results[[contrast]], list_results[[contrast]]$FDR < FDR2use & list_results[[contrast]][7] < 0 & abs(list_results[[contrast]][7]) > 1) 
	list_de_sex[[contrast]] <- subset(list_results[[contrast]], list_results[[contrast]]$FDR < FDR2use & abs(list_results[[contrast]][7]) > 1) 
}

merged1to2 <- merge(list_de[[1]][1], list_de[[2]][1], by.x='gene', by.y='gene', all=T)
merged1to3 <- merge(merged1to2, list_de[[3]][1], by.x='gene', by.y='gene', all=T)
merged1to4 <- merge(merged1to3, list_de[[4]][1], by.x='gene', by.y='gene', all=T)

merged1to2_up <- merge(list_de_up[[1]][1], list_de_up[[2]][1], by.x='gene', by.y='gene', all=T)
merged1to3_up <- merge(merged1to2_up, list_de_up[[3]][1], by.x='gene', by.y='gene', all=T)
merged1to4_up <- merge(merged1to3_up, list_de_up[[4]][1], by.x='gene', by.y='gene', all=T)

merged1to2_down <- merge(list_de_down[[1]][1], list_de_down[[2]][1], by.x='gene', by.y='gene', all=T)
merged1to3_down <- merge(merged1to2_down, list_de_down[[3]][1], by.x='gene', by.y='gene', all=T)
merged1to4_down <- merge(merged1to3_down, list_de_down[[4]][1], by.x='gene', by.y='gene', all=T)


merged1to2_sex <- merge(data.frame(list_de_sex[[1]][1], list_de_sex[[1]][7]), data.frame(list_de_sex[[2]][1], list_de_sex[[2]][7]), by.x='gene', by.y='gene', all=T)
merged1to3_sex <- merge(merged1to2_sex, data.frame(list_de_sex[[3]][1], list_de_sex[[3]][7]), by.x='gene', by.y='gene', all=T)
merged1to4_sex <- merge(merged1to3_sex, data.frame(list_de[[4]][1], list_de[[4]][7]), by.x='gene', by.y='gene', all=T)


merged1to5 <- merge(merged1to4, list_de[[5]][1], by.x='gene', by.y='gene', all=T)
merged1to6 <- merge(merged1to5, list_de[[6]][1], by.x='gene', by.y='gene', all=T)
merged1to7 <- merge(merged1to6, list_de[[7]][1], by.x='gene', by.y='gene', all=T)
merged1to8 <- merge(merged1to7, list_de[[8]][1], by.x='gene', by.y='gene', all=T)

annotated_de <- merge(merged1to4, annotation, by.x='gene', by.y='gid', all=F) # use for courthsip, 

# write.table(annotated_de$gene, file=file.path("~/Desktop/out.txt"), quote=F, row.names=F, sep="\t")

if (plottype == 'EM') {
annotated_de <- merge(merged1to8, annotation, by.x='gene', by.y='gid', all=F) # fix for 8 contrasts
merged1to8$DE <- 'yes'
annotated2_de <- merge(merged1to8, annotation, by.x='gene', by.y='gid', all=T) # fix for 8 contrasts
annotated2_de$DE <- as.character(annotated2_de$DE) 
annotated2_de$DE[is.na(annotated2_de$DE)] <- 'NA' 
annotated2_de$DE <- as.factor(annotated2_de$DE)

annotated3 <- merge(annotated2_de, dpse_dmel_data, by.x='gene', by.y='dpse', all=F)
annotated3$dmel <- as.character(annotated3$dmel) 
annotated3$dmel[is.na(annotated3$dmel)] <- 'NA' 
annotated3$dmel <- as.factor(annotated3$dmel)

annotated4 <- merge(annotated3, gerrard_data, by.x='dmel', by.y='id', all=T)
annotated4$tissue <- as.character(annotated4$tissue ) 
annotated4$tissue [is.na(annotated4$tissue )] <- 'NA' 
annotated4$tissue  <- as.factor(annotated4$tissue )

annotated4$updown <- as.character(annotated4$updown ) 
annotated4$updown [is.na(annotated4$updown )] <- 'NA' 
annotated4$updown  <- as.factor(annotated4$updown )

annotated4$gene <- as.character(annotated4$gene  ) 
annotated4$gene  [is.na(annotated4$gene  )] <- 'NA' 
annotated4$gene   <- as.factor(annotated4$gene  )

set1 <- unique(as.character(subset(annotated4, annotated4$gene!='NA')$gene) )
set2 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$DE=='yes')$gene) )
set3 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$dmel=='NA')$gene) )
set4 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$tissue!='NA')$gene) )

source('~/git/feminisation_direction/scripts/vennFunctions.r')

pdf(file.path(outpath, paste('Venn_4.pdf', sep="")), width=16, height=6)
plot.new()
par(mfrow=c(1,3)) 
par(mar=c(5,5,4,3), oma=c(0,0,4,0))
Venn3(set1, set2, set3, c('All', 'DE', 'Dpse only'), 'DE - 8 contrasts')
Venn3(set1, set2, set4, c('All', 'DE', 'Gerrard'), 'DE - 8 contrasts')
Venn4(set1, set2, set3, set4, c('All', 'DE', 'noDmel', 'Gerrard'), 'DE - 8 contrasts')
mtext('All EM genes', outer = TRUE, cex = 1)
dev.off()

# these numbers are correct and are more than table 4 of chromosomal location chisq tests because they include genes not assigned to major Muller elements.
de_dpsedmel <- 390
de_dpse <- 60
node_dpsedmel <- 11616
node_dpse <- 4283

de_dpsedmel <- 71
de_dpse <- 5
node_dpsedmel <- 11545
node_dpse <- 385



total_significant <- de_dpsedmel + de_dpse
total_dmel <- sum(de_dpse + node_dpse)
total_nodmel <- sum(de_dpsedmel + node_dpsedmel)
exp_dmel <- round(total_significant * total_dmel/(total_dmel + total_nodmel))
exp_nodmel <- round(total_significant * total_nodmel/(total_dmel + total_nodmel))


Obs <- c(de_dpse, de_dpsedmel)
Exp <- c(exp_dmel, exp_nodmel)
chi_obs_data <- rbind(Obs, Exp)
colnames(chi_obs_data) <- c('Dpse only', 'Dpse and Dmel')
chi_test <- chisq.test(c(de_dpse, de_dpsedmel), p=c(total_dmel/(total_dmel + total_nodmel), total_nodmel/(total_dmel + total_nodmel)))
pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

cbblue <- '#0072B2'
cblblue <- '#56B4E9'

par(mar=c(5,5,4,3), oma=c(0,0,2,0))
mp <- barplot(chi_obs_data, col=c(cbgreen ,cbblue), main = 'DE - 8 contrasts', beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
names(vals) <- LETTERS[1:4]
text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)

chi_test
}

if (plottype == 'sex') {
annotated_de <- merge(merged1to4, annotation, by.x='gene', by.y='gid', all=F)
merged1to4_sex$DE <- 'yes'
annotated2_de <- merge(merged1to4_sex, annotation, by.x='gene', by.y='gid', all=T) # fix for 8 contrasts
annotated2_de$DE <- as.character(annotated2_de$DE) 
annotated2_de$DE[is.na(annotated2_de$DE)] <- 'NA' 
annotated2_de$DE <- as.factor(annotated2_de$DE)

annotated3 <- merge(annotated2_de, dpse_dmel_data, by.x='gene', by.y='dpse', all=F)
annotated3$dmel <- as.character(annotated3$dmel) 
annotated3$dmel[is.na(annotated3$dmel)] <- 'NA' 
annotated3$dmel <- as.factor(annotated3$dmel)

annotated4 <- merge(annotated3, gerrard_data, by.x='dmel', by.y='id', all=T)
annotated4$tissue <- as.character(annotated4$tissue ) 
annotated4$tissue [is.na(annotated4$tissue )] <- 'NA' 
annotated4$tissue  <- as.factor(annotated4$tissue )

annotated4$updown <- as.character(annotated4$updown ) 
annotated4$updown [is.na(annotated4$updown )] <- 'NA' 
annotated4$updown  <- as.factor(annotated4$updown )

annotated4$gene <- as.character(annotated4$gene  ) 
annotated4$gene  [is.na(annotated4$gene  )] <- 'NA' 
annotated4$gene   <- as.factor(annotated4$gene  )




annotated4_fem1 <- subset(annotated4, annotated4[3] < 0 | is.na(annotated4[3]) ) 
annotated4_fem2 <- subset(annotated4_fem1, annotated4_fem1[4] < 0 | is.na(annotated4_fem1[4]) ) 
annotated4_fem3 <- subset(annotated4_fem2, annotated4_fem2[5] < 0 | is.na(annotated4_fem2[5]) ) 
annotated4_fem4 <- subset(annotated4_fem3, annotated4_fem3[6] < 0 | is.na(annotated4_fem3[6]) ) 
annotated4_fem5 <- subset(annotated4_fem4, annotated4_fem4[6] < 0 | annotated4_fem4[5] < 0 | annotated4_fem4[4] < 0 | annotated4_fem4[3] < 0  ) 
annotated4_fem5$bias <- 'female'

annotated4_male1 <- subset(annotated4, annotated4[3] > 0 | is.na(annotated4[3]) ) 
annotated4_male2 <- subset(annotated4_male1, annotated4_male1[4] > 0 | is.na(annotated4_male1[4]) ) 
annotated4_male3 <- subset(annotated4_male2, annotated4_male2[5] > 0 | is.na(annotated4_male2[5]) ) 
annotated4_male4 <- subset(annotated4_male3, annotated4_male3[6] > 0 | is.na(annotated4_male3[6]) ) 
annotated4_male5 <- subset(annotated4_male4, annotated4_male4[6] > 0 | annotated4_male4[5] > 0 | annotated4_male4[4] > 0 | annotated4_male4[3] > 0  ) 
annotated4_male5$bias <- 'male'

annotated4_ub <- subset(annotated4_fem4, is.na(annotated4_fem4[6]) & is.na(annotated4_fem4[5]) & is.na(annotated4_fem4[4]) & is.na(annotated4_fem4[3]) ) 
annotated4_ub$bias <- 'ub'

annotated5 <- rbind(annotated4_fem5, annotated4_male5, annotated4_ub)
annotated5$bias <- factor(annotated5$bias)

set1 <- unique(as.character(subset(annotated4, annotated4$gene!='NA')$gene) )
set2 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$DE=='yes')$gene) )
set3 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$dmel=='NA')$gene) )
set4 <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$tissue!='NA')$gene) )

set1s <- unique(as.character(subset(annotated5, annotated5$gene!='NA')$gene) )
set2s <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$DE=='yes')$gene) )
set3s <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$dmel=='NA')$gene) )
set4s <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$tissue!='NA')$gene) )



merged1to4_up$DE <- 'yes'
annotated2_de <- merge(merged1to4_up, annotation, by.x='gene', by.y='gid', all=T) # fix for 8 contrasts
annotated2_de$DE <- as.character(annotated2_de$DE) 
annotated2_de$DE[is.na(annotated2_de$DE)] <- 'NA' 
annotated2_de$DE <- as.factor(annotated2_de$DE)

annotated3 <- merge(annotated2_de, dpse_dmel_data, by.x='gene', by.y='dpse', all=F)
annotated3$dmel <- as.character(annotated3$dmel) 
annotated3$dmel[is.na(annotated3$dmel)] <- 'NA' 
annotated3$dmel <- as.factor(annotated3$dmel)

annotated4 <- merge(annotated3, gerrard_data, by.x='dmel', by.y='id', all=T)
annotated4$tissue <- as.character(annotated4$tissue ) 
annotated4$tissue [is.na(annotated4$tissue )] <- 'NA' 
annotated4$tissue  <- as.factor(annotated4$tissue )

annotated4$updown <- as.character(annotated4$updown ) 
annotated4$updown [is.na(annotated4$updown )] <- 'NA' 
annotated4$updown  <- as.factor(annotated4$updown )

annotated4$gene <- as.character(annotated4$gene  ) 
annotated4$gene  [is.na(annotated4$gene  )] <- 'NA' 
annotated4$gene   <- as.factor(annotated4$gene  )

set1_up <- unique(as.character(subset(annotated4, annotated4$gene!='NA')$gene) )
set2_up <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$DE=='yes')$gene) )
set3_up <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$dmel=='NA')$gene) )
set4_up <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$tissue!='NA')$gene) )

set1_ups <- unique(as.character(subset(annotated5, annotated5$gene!='NA')$gene) )
set2_ups <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$DE=='yes')$gene) )
set3_ups <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$dmel=='NA')$gene) )
set4_ups <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$tissue!='NA')$gene) )



merged1to4_down$DE <- 'yes'
annotated2_de <- merge(merged1to4_down, annotation, by.x='gene', by.y='gid', all=T) # fix for 8 contrasts
annotated2_de$DE <- as.character(annotated2_de$DE) 
annotated2_de$DE[is.na(annotated2_de$DE)] <- 'NA' 
annotated2_de$DE <- as.factor(annotated2_de$DE)

annotated3 <- merge(annotated2_de, dpse_dmel_data, by.x='gene', by.y='dpse', all=F)
annotated3$dmel <- as.character(annotated3$dmel) 
annotated3$dmel[is.na(annotated3$dmel)] <- 'NA' 
annotated3$dmel <- as.factor(annotated3$dmel)

annotated4 <- merge(annotated3, gerrard_data, by.x='dmel', by.y='id', all=T)
annotated4$tissue <- as.character(annotated4$tissue ) 
annotated4$tissue [is.na(annotated4$tissue )] <- 'NA' 
annotated4$tissue  <- as.factor(annotated4$tissue )

annotated4$updown <- as.character(annotated4$updown ) 
annotated4$updown [is.na(annotated4$updown )] <- 'NA' 
annotated4$updown  <- as.factor(annotated4$updown )

annotated4$gene <- as.character(annotated4$gene  ) 
annotated4$gene  [is.na(annotated4$gene  )] <- 'NA' 
annotated4$gene   <- as.factor(annotated4$gene  )

set1_down <- unique(as.character(subset(annotated4, annotated4$gene!='NA')$gene) )
set2_down <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$DE=='yes')$gene) )
set3_down <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$dmel=='NA')$gene) )
set4_down <- unique(as.character(subset(annotated4, annotated4$gene!='NA' & annotated4$tissue!='NA')$gene) )

set1_downs <- unique(as.character(subset(annotated5, annotated5$gene!='NA')$gene) )
set2_downs <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$DE=='yes')$gene) )
set3_downs <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$dmel=='NA')$gene) )
set4_downs <- unique(as.character(subset(annotated5, annotated5$gene!='NA' & annotated5$tissue!='NA')$gene) )


source('~/git/feminisation_direction/scripts/vennFunctions.r')


pdf(file.path(outpath, paste('Venn_4.pdf', sep="")), width=16, height=12)
par(mfrow=c(2,3)) 
par(mar=c(5,5,4,3), oma=c(0,0,4,0))
Venn3(set1, set2, set3, c('All', 'DE', 'Dpse only'), 'Sex biased')
Venn3(set1_up, set2_up, set3_up, c('All', 'DE', 'Dpse only'), 'DE - male')
Venn3(set1_down, set2_down, set3_down, c('All', 'DE', 'Dpse only'), 'DE - female')

Venn3(set1s, set2s, set3s, c('All', 'DE', 'Dpse only'), 'Sex biased strict')
Venn3(set1_ups, set2_ups, set3_ups, c('All', 'DE', 'Dpse only'), 'DE - male strict')
Venn3(set1_downs, set2_downs, set3_downs, c('All', 'DE', 'Dpse only'), 'DE - female strict')

# Venn4(set1, set2, set3, set4, c('All', 'DE', 'noDmel', 'Gerrard'), 'DE - ntrasts')
# Venn4(set1_up, set2_up, set3_up, set4_up, c('All', 'DE', 'noDmel', 'Gerrard'), 'DE - ntrasts')
# Venn4(set1_down, set2_down, set3_down, set4_down, c('All', 'DE', 'noDmel', 'Gerrard'), 'DE - ntrasts')

mtext('All sex biased genes', outer = TRUE, cex = 1)
dev.off()

# All
de_dpsedmel <- 8770
de_dpse <- 2180
node_dpsedmel <- 3236
node_dpse <- 2163
de_dpse/(de_dpsedmel+de_dpse)
(de_dpse+node_dpse)/(de_dpsedmel+de_dpse+node_dpsedmel+node_dpse)
# X-squared = 218.45, df = 1, p-value < 2.2e-16 , 0.20 vs 0.26

# all strict
de_dpsedmel <- 8686
de_dpse <- 2115
node_dpsedmel <- 3236
node_dpse <- 2163
de_dpse/(de_dpsedmel+de_dpse)
(de_dpse+node_dpse)/(de_dpsedmel+de_dpse+node_dpsedmel+node_dpse)
# X-squared = 258.95, df = 1, p-value < 2.2e-16, 0.195 vs 0.26


# male
de_dpsedmel <- 4530
de_dpse <- 1691
node_dpsedmel <- 7476
node_dpse <- 2652
de_dpse/(de_dpsedmel+de_dpse)
(de_dpse+node_dpse)/(de_dpsedmel+de_dpse+node_dpsedmel+node_dpse)
# X-squared = 1.2172, df = 1, p-value = 0.2699, 0.27 vs 0.27

# female
de_dpsedmel <- 4171
de_dpse <- 532
node_dpsedmel <- 7835
node_dpse <- 3811
de_dpse/(de_dpsedmel+de_dpse)
(de_dpse+node_dpse)/(de_dpsedmel+de_dpse+node_dpsedmel+node_dpse)
# X-squared = 560.85, df = 1, p-value < 2.2e-16, 0.11 vs 0.27


total_significant <- de_dpsedmel + de_dpse
total_dmel <- sum(de_dpse + node_dpse)
total_nodmel <- sum(de_dpsedmel + node_dpsedmel)
exp_dmel <- round(total_significant * total_dmel/(total_dmel + total_nodmel))
exp_nodmel <- round(total_significant * total_nodmel/(total_dmel + total_nodmel))


Obs <- c(de_dpse, de_dpsedmel)
Exp <- c(exp_dmel, exp_nodmel)
chi_obs_data <- rbind(Obs, Exp)
colnames(chi_obs_data) <- c('Dpse only', 'Dpse and Dmel')
chi_test <- chisq.test(c(de_dpse, de_dpsedmel), p=c(total_dmel/(total_dmel + total_nodmel), total_nodmel/(total_dmel + total_nodmel)))
pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

cbblue <- '#0072B2'
cblblue <- '#56B4E9'

par(mar=c(5,5,4,3), oma=c(0,0,2,0))
mp <- barplot(chi_obs_data, col=c(cbgreen ,cbblue), main = 'DE - 8 contrasts', beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
names(vals) <- LETTERS[1:4]
text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)

annotated5$nohomolog <- factor(annotated5$dmel=='NA')
annotated6 <- subset(annotated5, annotated5$gene!='NA')
cbind(annotated6$DE=='yes', annotated6$DE!='yes')

mod1 <- glm(cbind(DE=='yes', DE!='yes') ~ nohomolog + bias + nohomolog:bias, data=annotated6, family=binomial )
summary(mod1)
anova(mod1, test='Chisq')

# numde <- c(

}




mapdata_de <- subset(annotated_de, annotated_de$chr!='NA')

mapdata_all <- subset(annotation, annotation$chr!='NA' )


mapLength <- sum(max(mapdata_all$end[mapdata_all$chr=='2_'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='3_'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XL_group1a'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XL_group1e'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XL_group3a'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XL_group3b'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XR_group3a'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XR_group5'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XR_group6'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='XR_group8'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='4_group1'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='4_group2'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='4_group3'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='4_group4'], na.rm=T), max(mapdata_all$end[mapdata_all$chr=='4_group5'], na.rm=T))

chr_XL_group1a_end <- max(mapdata_all$end[mapdata_all$chr=='XL_group1a'], na.rm=T)
chr_XL_group1e_end <- chr_XL_group1a_end+(max(mapdata_all$end[mapdata_all$chr=='XL_group1e'], na.rm=T))
chr_XL_group3a_end <- chr_XL_group1e_end+(max(mapdata_all$end[mapdata_all$chr=='XL_group3a'], na.rm=T))
chr_XL_group3b_end <- chr_XL_group3a_end+(max(mapdata_all$end[mapdata_all$chr=='XL_group3b'], na.rm=T))

chr_XL_end <- chr_XL_group3b_end 

chr_XR_group3a_end <- max(mapdata_all$end[mapdata_all$chr=='XR_group3a'], na.rm=T) + chr_XL_end
chr_XR_group5_end <- chr_XR_group3a_end+(max(mapdata_all$end[mapdata_all$chr=='XR_group5'], na.rm=T))
chr_XR_group6_end <- chr_XR_group5_end+(max(mapdata_all$end[mapdata_all$chr=='XR_group6'], na.rm=T))
chr_XR_group8_end <- chr_XR_group6_end+(max(mapdata_all$end[mapdata_all$chr=='XR_group8'], na.rm=T))
chr_XR_end <- chr_XR_group8_end

chr_2_end <- max(mapdata_all$end[mapdata_all$chr=='2_'], na.rm=T) + chr_XR_end
chr_3_end <- max(mapdata_all$end[mapdata_all$chr=='3_'], na.rm=T) + chr_2_end

chr_4_group1_end <- max(mapdata_all$end[mapdata_all$chr=='4_group1'], na.rm=T) + chr_3_end
chr_4_group2_end <- chr_4_group1_end+(max(mapdata_all$end[mapdata_all$chr=='4_group2'], na.rm=T))
chr_4_group3_end <- chr_4_group2_end+(max(mapdata_all$end[mapdata_all$chr=='4_group3'], na.rm=T))
chr_4_group4_end <- chr_4_group3_end+(max(mapdata_all$end[mapdata_all$chr=='4_group4'], na.rm=T))
chr_4_group5_end <- chr_4_group4_end+(max(mapdata_all$end[mapdata_all$chr=='4_group5'], na.rm=T))


# Chisq comparison Chr XL, XR, 4 with others
chi_out <- data.frame(c(0,0,0,0,0,0), c(0,0,0,0,0,0), c(0,0,0,0,0,0), c(0,0,0,0,0,0), c(0,0,0,0,0,0), c(0,0,0,0,0,0))
colnames(chi_out) <- c('O1', 'E1', 'O2', 'E2', 'statistic', 'p')
row.names(chi_out) <- c('2_', '3_', '4_', 'X', 'XR_', 'XL_')
count <- 0
	for (chrom in c('2_', '3_', '4_', 'X', 'XR_', 'XL_') ) {
		## estimation of loop counts
		count <- count+1
		de_chrfocus <- sum(annotated_de$chr[grepl(chrom, annotated_de$chr)] != 'na', na.rm=T) 
		de_chrother <- sum(annotated_de$chr[!grepl('Unknown*', annotated_de$chr) & !grepl(chrom, annotated_de$chr)] != 'na', na.rm=T)
		nonde_chrfocus <- sum(mapdata_all$chr[grepl(chrom, mapdata_all$chr)] != 'na', na.rm=T)
		nonde_chrother <- sum(mapdata_all$chr[!grepl('Unknown*', mapdata_all$chr) & !grepl(chrom, mapdata_all$chr)] !='na')

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
			chi_out$O1[count] <- round(chi_test$observed[1] , digits=2)
			chi_out$O2[count] <- round(chi_test$observed[2] , digits=2)
			chi_out$E1[count] <- round(chi_test$expected[1] , digits=2)
			chi_out$E2[count] <- round(chi_test$expected[2], digits=2)
			chi_out$statistic[count] <- round(chi_test$statistic[1] , digits=2)
			chi_out$p[count] <- round(chi_test$p.value[1] , digits=2)

			pdf(file.path(datapath, paste('Chisq_Chr_', chrom, 'FDR_', FDR2use, '_all_', plottype, '.pdf', sep="")), width=9, height=9)
			par(mar=c(5,5,4,3), oma=c(0,0,2,0))
			mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste('All DE ', '\nFDR=',FDR2use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
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

		mtext(paste('Chrom', chrom, 'vs rest'), outer = TRUE, cex = 1.5)
		dev.off()
		}
	## estimation of X counts
		de_XL <- sum(annotated_de$chr[grepl('XL_', annotated_de$chr)] != 'na', na.rm=T) 
		nonde_XL <- sum(mapdata_all$chr[grepl('XL_', mapdata_all$chr)] != 'na', na.rm=T) 
		de_XR <- sum(annotated_de$chr[grepl('XR_', annotated_de$chr)] != 'na', na.rm=T) 
		nonde_XR <- sum(mapdata_all$chr[grepl('XR_', mapdata_all$chr)] != 'na', na.rm=T) 
		de_X <- de_XL + de_XR
		nonde_X <- nonde_XL + nonde_XR

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
					pdf(file.path(datapath, paste('Chisq_Chr_', chrom,'vs_XR_FDR_', FDR2use, '_all_', plottype, '.pdf', sep="")), width=9, height=9)
# 					par(mfrow=c(1,3))
					par(mar=c(5,5,4,3), oma=c(0,0,2,0))
					mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste('All DE ', '\nFDR=',FDR2use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
					if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
					vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
					names(vals) <- LETTERS[1:4]
				text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
					if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
					if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
				}

				if (sum(de_XL + de_XR) > 10 & de_XL !=0 ) {
			mtext(paste('Chrom', chrom, 'vs XR'), outer = TRUE, cex = 1.5)
			dev.off()
			}
		}
		write.table(chi_out, file=file.path(datapath, paste(plottype,"_FDR_",FDR2use, "_chi_tbl.txt", sep='')), quote=F, row.names=T, sep="\t")
	}

