## replace EMMVH

datapath <- "~/git/feminisation_direction/input"
outpath <- "~/git/feminisation_direction/output/correlations"
subsetpath <- "~/git/feminisation_direction/output/subsets"
dir.create(file.path(outpath))

HollData <- read.table(file.path(datapath, "GSE50915_normalized_counts.txt"), header=T)
str(HollData)

mdata <- read.table(file.path(datapath, 'GSE39742_series_matrix.txt'), header=T)
str(mdata) 
fdata <- read.table(file.path(datapath, 'fbgn_fbtr.txt'), header=T)
str(fdata)

mfdata <- merge(mdata, fdata, by.x="ID_REF", by.y="FBtr", all=F )

dpseDmel <- read.table(file.path(datapath,"dpseDmel.txt"), header=T)
str(dpseDmel)

EMMVB <- read.table(file.path(subsetpath, "EMMVB/EMMVB_pseudo_counts.txt"), header=T)
EMFVB <- read.table(file.path(subsetpath, "EMFVB/EMFVB_pseudo_counts.txt"), header=T)
str(EMMVB)
str(EMFVB)

merged1_m <- merge(dpseDmel, EMMVB, by.x="dpse", by.y="gene", all=F)
merged2_m <- merge(merged1_m, mfdata, by.x="dmel", by.y="FBgn", all=F)

merged1_f <- merge(dpseDmel, EMFVB, by.x="dpse", by.y="gene", all=F)
merged2_f <- merge(merged1_f, mfdata, by.x="dmel", by.y="FBgn", all=F)

### for females
EM_Vel_all_f <- (merged2_f$E_F_V_B + merged2_f$E_F_V_B.1 + merged2_f$E_F_V_B.2 + merged2_f$E_F_V_B.3)/4 - (merged2_f$M_F_V_B + merged2_f$M_F_V_B.1 + merged2_f$M_F_V_B.2 + merged2_f$M_F_V_B.3)/4

EM_Holl_all_f <- (merged2_f$E1F+ merged2_f$E2F + merged2_f$E3F)/3 - (merged2_f$M1F + merged2_f$M2F + merged2_f$M3F)/3

### for males
EM_Vel_all_m <- (merged2_m$E_M_V_B + merged2_m$E_M_V_B.1 + merged2_m$E_M_V_B.2 + merged2_m$E_M_V_B.3)/4 - (merged2_m$M_M_V_B + merged2_m$M_M_V_B.1 + merged2_m$M_M_V_B.2 + merged2_m$M_M_V_B.3)/4

EM_Holl_all_m <- (merged2_m$E1M + merged2_m$E2M + merged2_m$E3M)/3 - (merged2_m$M1M + merged2_m$M2M + merged2_m$M3M)/3

cor(EM_Vel_all_m, EM_Holl_all_m)
cor.test(EM_Vel_all_m, EM_Holl_all_m)$p.value
paste_all_m = paste("cor =", round(cor(EM_Vel_all_m, EM_Holl_all_m),digits=2), ", p = ", round(cor.test(EM_Vel_all_m, EM_Holl_all_m)$p.value, digits=2), sep=" ")

cor(EM_Vel_all_f, EM_Holl_all_f)
cor.test(EM_Vel_all_f, EM_Holl_all_f)$p.value
paste_all_f = paste("cor =", round(cor(EM_Vel_all_f, EM_Holl_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_all_f, EM_Holl_all_f)$p.value, digits=2), sep=" ")

merged3_m <- subset(merged2_m, merged2_m$FDR < 0.1) 
merged3_f <- subset(merged2_f, merged2_f$FDR < 0.1) 


### for females
EM_Vel_DE_f <- (merged3_f$E_F_V_B + merged3_f$E_F_V_B.1 + merged3_f$E_F_V_B.2 + merged3_f$E_F_V_B.3)/4 - (merged3_f$M_F_V_B + merged3_f$M_F_V_B.1 + merged3_f$M_F_V_B.2 + merged3_f$M_F_V_B.3)/4

EM_Holl_DE_f <- (merged3_f$E1F + merged3_f$E2F + merged3_f$E3F)/3 - (merged3_f$M1F + merged3_f$M2F + merged3_f$M3F)/3

### for males
EM_Vel_DE_m <- (merged3_m$E_M_V_B + merged3_m$E_M_V_B.1 + merged3_m$E_M_V_B.2 + merged3_m$E_M_V_B.3)/4 - (merged3_m$M_M_V_B + merged3_m$M_M_V_B.1 + merged3_m$M_M_V_B.2 + merged3_m$M_M_V_B.3)/4

EM_Holl_DE_m <- (merged3_m$E1M + merged3_m$E2M + merged3_m$E3M)/3 - (merged3_m$M1M + merged3_m$M2M + merged3_m$M3M)/3

cor(EM_Vel_DE_f, EM_Holl_DE_f)
cor.test(EM_Vel_DE_f, EM_Holl_DE_f)$p.value
cor(EM_Vel_DE_m, EM_Holl_DE_m)
cor.test(EM_Vel_DE_m, EM_Holl_DE_m)$p.value

paste_DE_f = paste("cor =", round(cor(EM_Vel_DE_f, EM_Holl_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_DE_f, EM_Holl_DE_f)$p.value, digits=2), sep=" ")
paste_DE_m = paste("cor =", round(cor(EM_Vel_DE_m, EM_Holl_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_DE_m, EM_Holl_DE_m)$p.value, digits=2), sep=" ")

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_all_f, EM_Holl_all_f, xlim=c(-9000,10000), ylim=c(-5,5), xlab="Mean E - M (Dpse)", ylab="Mean E - M (Dmel)", main="All genes E - M (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_all_f, cex=0.8)
plot(EM_Vel_DE_f, EM_Holl_DE_f, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes E - M (females)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_f, cex=0.8)

plot(EM_Vel_all_m, EM_Holl_all_m, xlab="Mean E - M (Dpse)", ylim=c(-5,5), ylab="Mean E - M (Dmel)", main="All genes E - M (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_m, cex=0.8)
plot(EM_Vel_DE_m, EM_Holl_DE_m, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes E - M (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_m, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VB_dpse_vs_dmel_EminusM.pdf', sep="")), width=12, height=12)
dev.off()




### males E DE only
EM_Vel_E_DE_m <- (merged3_m$E_M_V_B + merged3_m$E_M_V_B.1 + merged3_m$E_M_V_B.2 + merged3_m$E_M_V_B.3)/4
EM_Holl_E_DE_m <- (merged3_m$E1M + merged3_m$E2M + merged3_m$E3M)/3

cor(EM_Vel_E_DE_m, EM_Holl_E_DE_m)
cor.test(EM_Vel_E_DE_m, EM_Holl_E_DE_m)$p.value

paste_E_DE_m = paste("cor =", round(cor(EM_Vel_E_DE_m, EM_Holl_E_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_E_DE_m, EM_Holl_E_DE_m)$p.value, digits=2), sep=" ")

EM_Vel_E_all_m <- (merged2_m$E_M_V_B + merged2_m$E_M_V_B.1 + merged2_m$E_M_V_B.2 + merged2_m$E_M_V_B.3)/4
EM_Holl_E_all_m <- (merged2_m$E1M + merged2_m$E2M + merged2_m$E3M)/3

cor(EM_Vel_E_all_m, EM_Holl_E_all_m)
cor.test(EM_Vel_E_all_m, EM_Holl_E_all_m)$p.value

paste_E_all_m = paste("cor =", round(cor(EM_Vel_E_all_m, EM_Holl_E_all_m),digits=2), ", p = ", round(cor.test(EM_Vel_E_all_m, EM_Holl_E_all_m)$p.value, digits=2), sep=" ")

### males M DE only
EM_Vel_M_DE_m <- (merged3_m$M_M_V_B + merged3_m$M_M_V_B.1 + merged3_m$M_M_V_B.2 + merged3_m$M_M_V_B.3)/4
EM_Holl_M_DE_m <- (merged3_m$M1M + merged3_m$M2M + merged3_m$M3M)/3

cor(EM_Vel_M_DE_m, EM_Holl_M_DE_m)
cor.test(EM_Vel_M_DE_m, EM_Holl_M_DE_m)$p.value

paste_M_DE_m = paste("cor =", round(cor(EM_Vel_M_DE_m, EM_Holl_M_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_M_DE_m, EM_Holl_M_DE_m)$p.value, digits=2), sep=" ")

EM_Vel_M_all_m <- (merged2_m$M_M_V_B + merged2_m$M_M_V_B.1 + merged2_m$M_M_V_B.2 + merged2_m$M_M_V_B.3)/4
EM_Holl_M_all_m <- (merged2_m$M1M + merged2_m$M2M + merged2_m$M3M)/3

cor(EM_Vel_M_all_m, EM_Holl_M_all_m)
cor.test(EM_Vel_M_all_m, EM_Holl_M_all_m)$p.value

paste_M_all_m = paste("cor =", round(cor(EM_Vel_M_all_m, EM_Holl_M_all_m),digits=2), ", p = ", round(cor.test(EM_Vel_M_all_m, EM_Holl_M_all_m)$p.value, digits=2), sep=" ")

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_E_all_m, EM_Holl_E_all_m, xlim=c(0,50000), ylim=c(0,20), xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="All genes E (males)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_E_all_m, cex=0.8)
plot(EM_Vel_E_DE_m, EM_Holl_E_DE_m, xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="DE E genes (males)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_E_DE_m, cex=0.8)

plot(EM_Vel_M_all_m, EM_Holl_M_all_m, xlim=c(0,50000), ylim=c(0,20), xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="All M genes (males)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_M_all_m, cex=0.8)
plot(EM_Vel_M_DE_m, EM_Holl_M_DE_m, xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="DE M genes (males)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_M_DE_m, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VB_dpse_vs_dmel_males.pdf', sep="")), width=12, height=12)
dev.off()





### females E DE only
EM_Vel_E_DE_f <- (merged3_f$E_F_V_B + merged3_f$E_F_V_B.1 + merged3_f$E_F_V_B.2 + merged3_f$E_F_V_B.3)/4
EM_Holl_E_DE_f <- (merged3_f$E1M + merged3_f$E2M + merged3_f$E3M)/3

cor(EM_Vel_E_DE_f, EM_Holl_E_DE_f)
cor.test(EM_Vel_E_DE_f, EM_Holl_E_DE_f)$p.value

paste_E_DE_f = paste("cor =", round(cor(EM_Vel_E_DE_f, EM_Holl_E_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_E_DE_f, EM_Holl_E_DE_f)$p.value, digits=2), sep=" ")

EM_Vel_E_all_f <- (merged2_f$E_F_V_B + merged2_f$E_F_V_B.1 + merged2_f$E_F_V_B.2 + merged2_f$E_F_V_B.3)/4
EM_Holl_E_all_f <- (merged2_f$E1M + merged2_f$E2M + merged2_f$E3M)/3

cor(EM_Vel_E_all_f, EM_Holl_E_all_f)
cor.test(EM_Vel_E_all_f, EM_Holl_E_all_f)$p.value

paste_E_all_f = paste("cor =", round(cor(EM_Vel_E_all_f, EM_Holl_E_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_E_all_f, EM_Holl_E_all_f)$p.value, digits=2), sep=" ")

### females M DE only
EM_Vel_f_DE_f <- (merged3_f$M_F_V_B + merged3_f$M_F_V_B.1 + merged3_f$M_F_V_B.2 + merged3_f$M_F_V_B.3)/4
EM_Holl_f_DE_f <- (merged3_f$M1M + merged3_f$M2M + merged3_f$M3M)/3

cor(EM_Vel_f_DE_f, EM_Holl_f_DE_f)
cor.test(EM_Vel_f_DE_f, EM_Holl_f_DE_f)$p.value

paste_f_DE_f = paste("cor =", round(cor(EM_Vel_f_DE_f, EM_Holl_f_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_f_DE_f, EM_Holl_f_DE_f)$p.value, digits=2), sep=" ")

EM_Vel_M_all_f <- (merged2_f$M_F_V_B + merged2_f$M_F_V_B.1 + merged2_f$M_F_V_B.2 + merged2_f$M_F_V_B.3)/4
EM_Holl_f_all_f <- (merged2_f$M1M + merged2_f$M2M + merged2_f$M3M)/3

cor(EM_Vel_M_all_f, EM_Holl_f_all_f)
cor.test(EM_Vel_M_all_f, EM_Holl_f_all_f)$p.value

paste_f_all_f = paste("cor =", round(cor(EM_Vel_M_all_f, EM_Holl_f_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_M_all_f, EM_Holl_f_all_f)$p.value, digits=2), sep=" ")

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_E_all_f, EM_Holl_E_all_f, xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="All E genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_E_all_f, cex=0.8)
plot(EM_Vel_E_DE_f, EM_Holl_E_DE_f, xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="DE E genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_E_DE_f, cex=0.8)

plot(EM_Vel_M_all_f, EM_Holl_f_all_f, xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="All M genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_f_all_f, cex=0.8)
plot(EM_Vel_f_DE_f, EM_Holl_f_DE_f, xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="DE M genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_f_DE_f, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VB_dpse_vs_dmel_females.pdf', sep="")), width=12, height=12)
dev.off()













# VIRGIN HEADS

EMMVH <- read.table(file.path(subsetpath, "EMMVH/EMMVH_pseudo_counts.txt"), header=T)
EMFVH <- read.table(file.path(subsetpath, "EMFVH/EMFVH_pseudo_counts.txt"), header=T)
str(EMMVH)
str(EMFVH)


merged1_m <- merge(dpseDmel, EMMVH, by.x="dpse", by.y="gene", all=F)
merged2_m <- merge(merged1_m, HollData, by.x="dmel", by.y="gene", all=F)

merged1_f <- merge(dpseDmel, EMFVH, by.x="dpse", by.y="gene", all=F)
merged2_f <- merge(merged1_f, HollData, by.x="dmel", by.y="gene", all=F)

### for females
EM_Vel_all_f <- (merged2_f$E_F_V_H + merged2_f$E_F_V_H.1 + merged2_f$E_F_V_H.2 + merged2_f$E_F_V_H.3)/4 - (merged2_f$M_F_V_H + merged2_f$M_F_V_H.1 + merged2_f$M_F_V_H.2 + merged2_f$M_F_V_H.3)/4

EM_Holl_all_f <- (merged2_f$E1F+ merged2_f$E2F + merged2_f$E3F)/3 - (merged2_f$M1F + merged2_f$M2F + merged2_f$M3F)/3

### for males
EM_Vel_all_m <- (merged2_m$E_M_V_H + merged2_m$E_M_V_H.1 + merged2_m$E_M_V_H.2 + merged2_m$E_M_V_H.3)/4 - (merged2_m$M_M_V_H + merged2_m$M_M_V_H.1 + merged2_m$M_M_V_H.2 + merged2_m$M_M_V_H.3)/4

EM_Holl_all_m <- (merged2_m$E1M + merged2_m$E2M + merged2_m$E3M)/3 - (merged2_m$M1M + merged2_m$M2M + merged2_m$M3M)/3

#### remove outliers
EM_Vel_all_m <- data.frame(EM_Vel_all_m)
EM_Holl_all_m <- data.frame(EM_Holl_all_m)

xvel_m <- quantile(EM_Vel_all_m$EM_Vel_all_m, c(0.005,0.995) )
xhol_m  <- quantile(EM_Holl_all_m$EM_Holl_all_m, c(0.005,0.995) )

xvel_m_clean <- EM_Vel_all_m[EM_Vel_all_m$EM_Vel_all_m >=xvel_m[1] & EM_Vel_all_m$EM_Vel_all_m<=xvel_m[2] & EM_Holl_all_m$EM_Holl_all_m>=xhol_m[1] & EM_Holl_all_m$EM_Holl_all_m<=xhol_m[2],]

xhol_m_clean <- EM_Holl_all_m[EM_Vel_all_m$EM_Vel_all_m >=xvel_m[1] & EM_Vel_all_m$EM_Vel_all_m<=xvel_m[2] & EM_Holl_all_m$EM_Holl_all_m>=xhol_m[1] & EM_Holl_all_m$EM_Holl_all_m<=xhol_m[2],]


EM_Vel_all_f <- data.frame(EM_Vel_all_f)
EM_Holl_all_f <- data.frame(EM_Holl_all_f)

xvel_f <- quantile(EM_Vel_all_f$EM_Vel_all_f, c(0.005,0.995) )
xhol_f  <- quantile(EM_Holl_all_f$EM_Holl_all_f, c(0.005,0.995) )

xvel_f_clean <- EM_Vel_all_f[EM_Vel_all_f$EM_Vel_all_f >=xvel_f[1] & EM_Vel_all_f$EM_Vel_all_f<=xvel_f[2] & EM_Holl_all_f$EM_Holl_all_f>=xhol_f[1] & EM_Holl_all_f$EM_Holl_all_f<=xhol_f[2],]

xhol_f_clean <- EM_Holl_all_f[EM_Vel_all_f$EM_Vel_all_f >=xvel_f[1] & EM_Vel_all_f$EM_Vel_all_f<=xvel_f[2] & EM_Holl_all_f$EM_Holl_all_f>=xhol_f[1] & EM_Holl_all_f$EM_Holl_all_f<=xhol_f[2],]

## with outliers
cor(EM_Vel_all_m, EM_Holl_all_m)
cor.test(EM_Vel_all_m, EM_Holl_all_m)$p.value
# paste_all_m_outliers = paste("cor =", round(cor(xvel_m_clean, xhol_m_clean),digits=2), ", p = ", cor.test(xvel_m_clean, xhol_m_clean)$p.value, sep=" ")
paste_all_m_outliers = paste("cor =", round(cor(xvel_m_clean, xhol_m_clean),digits=2), ", p < 0.0001", sep=" ")

paste_all_m = paste("cor =", round(cor(EM_Vel_all_m, EM_Holl_all_m),digits=2), ", p = ", cor.test(EM_Vel_all_m, EM_Holl_all_m)$p.value, sep=" ")

cor(EM_Vel_all_f, EM_Holl_all_f)
cor.test(EM_Vel_all_f, EM_Holl_all_f)$p.value
# paste_all_f_outliers = paste("cor =", round(cor(xvel_f_clean, xhol_f_clean),digits=2), ", p = ", round(cor.test(xvel_f_clean, xhol_f_clean)$p.value, digits=2), sep=" ")
paste_all_f_outliers = paste("cor =", round(cor(xvel_f_clean, xhol_f_clean),digits=2), ", p < 0.0001", sep=" ")

paste_all_f = paste("cor =", round(cor(EM_Vel_all_f, EM_Holl_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_all_f, EM_Holl_all_f)$p.value, digits=2), sep=" ")


merged3_m <- subset(merged2_m, merged2_m$FDR < 0.1) 
merged3_f <- subset(merged2_f, merged2_f$FDR < 0.1) 



### for females
EM_Vel_DE_f <- (merged3_f$E_F_V_H + merged3_f$E_F_V_H.1 + merged3_f$E_F_V_H.2 + merged3_f$E_F_V_H.3)/4 - (merged3_f$M_F_V_H + merged3_f$M_F_V_H.1 + merged3_f$M_F_V_H.2 + merged3_f$M_F_V_H.3)/4

EM_Holl_DE_f <- (merged3_f$E1F + merged3_f$E2F + merged3_f$E3F)/3 - (merged3_f$M1F + merged3_f$M2F + merged3_f$M3F)/3

### for males
EM_Vel_DE_m <- (merged3_m$E_M_V_H + merged3_m$E_M_V_H.1 + merged3_m$E_M_V_H.2 + merged3_m$E_M_V_H.3)/4 - (merged3_m$M_M_V_H + merged3_m$M_M_V_H.1 + merged3_m$M_M_V_H.2 + merged3_m$M_M_V_H.3)/4

EM_Holl_DE_m <- (merged3_m$E1M + merged3_m$E2M + merged3_m$E3M)/3 - (merged3_m$M1M + merged3_m$M2M + merged3_m$M3M)/3

cor(EM_Vel_DE_f, EM_Holl_DE_f)
cor.test(EM_Vel_DE_f, EM_Holl_DE_f)$p.value
cor(EM_Vel_DE_m, EM_Holl_DE_m)
cor.test(EM_Vel_DE_m, EM_Holl_DE_m)$p.value


paste_DE_f = paste("cor =", round(cor(EM_Vel_DE_f, EM_Holl_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_DE_f, EM_Holl_DE_f)$p.value, digits=2), sep=" ")
paste_DE_m = paste("cor =", round(cor(EM_Vel_DE_m, EM_Holl_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_DE_m, EM_Holl_DE_m)$p.value, digits=2), sep=" ")

# only DE genes
# par(mfrow=c(1,2)) 
# par(mar=c(5,5,4,3))
# plot(EM_Vel_DE_f, EM_Holl_DE_f, xlab="Log2FC D. pseudoobscura", ylab="Log2FC D. melanogaster", main="DE genes high vs low sexual selection (females)", cex.main=1.4, cex.lab=1.2)
# legend('topleft', inset=0.05, legend=paste_DE_f, cex=1)

# plot(EM_Vel_DE_m, EM_Holl_DE_m, xlab="Log2FC D. pseudoobscura", ylab="Log2FC D. melanogaster", main="DE genes high vs low sexual selection (males)", cex.main=1.4, cex.lab=1.2)
# legend('topleft', inset=0.05, legend=paste_DE_m, cex=1)

# dev.copy(pdf,file.path(outpath, paste('VH_dpse_vs_dmel_EminusM.pdf', sep="")), width=12, height=6)
# dev.off()


# without outliers
par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(xvel_f_clean, xhol_f_clean, xlab="Mean E - M (Dpse)", ylab="Mean E - M (Dmel)", main="All genes high vs low sexual selection (females)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_f_outliers, cex=0.8)
plot(EM_Vel_DE_f, EM_Holl_DE_f, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes high vs low sexual selection (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomleft', inset=0.05, legend=paste_DE_f, cex=0.8)

plot(xvel_m_clean, xhol_m_clean, xlab="Mean E - M (Dpse)", ylab="Mean E - M (Dmel)", main="All genes high vs low sexual selection (males)", cex.main=1.5, cex.lab=1.2)
legend('bottomleft', inset=0.05, legend=paste_all_m_outliers, cex=0.8)
plot(EM_Vel_DE_m, EM_Holl_DE_m, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes high vs low sexual selection (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_m, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VH_dpse_vs_dmel_EminusM_wo_outliers.pdf', sep="")), width=12, height=12)
dev.off()


# with outliers
par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_all_f, EM_Holl_all_f, xlim=c(-14000,14000), ylim=c(-29000,25000), xlab="Mean E - M (Dpse)", ylab="Mean E - M (Dmel)", main="All genes high vs low sexual selection (females)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_f, cex=0.8)
plot(EM_Vel_DE_f, EM_Holl_DE_f, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes high vs low sexual selection (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomleft', inset=0.05, legend=paste_DE_f, cex=0.8)

plot(EM_Vel_all_m, EM_Holl_all_m, xlim=c(-21000,10000), ylim=c(-30000,13000), xlab="Mean E - M (Dpse)", ylab="Mean E - M (Dmel)", main="All genes high vs low sexual selection (males)", cex.main=1.5, cex.lab=1.2)
legend('bottomleft', inset=0.05, legend=paste_all_m, cex=0.8)
plot(EM_Vel_DE_m, EM_Holl_DE_m, xlab="E - M (Dpse)", ylab="E - M (Dmel)", main="DE genes high vs low sexual selection (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_m, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VH_dpse_vs_dmel_EminusM.pdf', sep="")), width=12, height=12)
dev.off()





### males E DE only
EM_Vel_E_DE_m <- (merged3_m$E_M_V_H + merged3_m$E_M_V_H.1 + merged3_m$E_M_V_H.2 + merged3_m$E_M_V_H.3)/4
EM_Holl_E_DE_m <- (merged3_m$E1M + merged3_m$E2M + merged3_m$E3M)/3

cor(EM_Vel_E_DE_m, EM_Holl_E_DE_m)
cor.test(EM_Vel_E_DE_m, EM_Holl_E_DE_m)$p.value

paste_E_DE_m = paste("cor =", round(cor(EM_Vel_E_DE_m, EM_Holl_E_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_E_DE_m, EM_Holl_E_DE_m)$p.value, digits=2), sep=" ")

EM_Vel_E_all_m <- (merged2_m$E_M_V_H + merged2_m$E_M_V_H.1 + merged2_m$E_M_V_H.2 + merged2_m$E_M_V_H.3)/4
EM_Holl_E_all_m <- (merged2_m$E1M + merged2_m$E2M + merged2_m$E3M)/3

cor(EM_Vel_E_all_m, EM_Holl_E_all_m)
cor.test(EM_Vel_E_all_m, EM_Holl_E_all_m)$p.value

paste_E_all_m = paste("cor =", round(cor(EM_Vel_E_all_m, EM_Holl_E_all_m),digits=2), ", p = ", round(cor.test(EM_Vel_E_all_m, EM_Holl_E_all_m)$p.value, digits=2), sep=" ")

### males M DE only
EM_Vel_M_DE_m <- (merged3_m$M_M_V_H + merged3_m$M_M_V_H.1 + merged3_m$M_M_V_H.2 + merged3_m$M_M_V_H.3)/4
EM_Holl_M_DE_m <- (merged3_m$M1M + merged3_m$M2M + merged3_m$M3M)/3

cor(EM_Vel_M_DE_m, EM_Holl_M_DE_m)
cor.test(EM_Vel_M_DE_m, EM_Holl_M_DE_m)$p.value

paste_M_DE_m = paste("cor =", round(cor(EM_Vel_M_DE_m, EM_Holl_M_DE_m),digits=2), ", p = ", round(cor.test(EM_Vel_M_DE_m, EM_Holl_M_DE_m)$p.value, digits=2), sep=" ")

EM_Vel_M_all_m <- (merged2_m$M_M_V_H + merged2_m$M_M_V_H.1 + merged2_m$M_M_V_H.2 + merged2_m$M_M_V_H.3)/4
EM_Holl_M_all_m <- (merged2_m$M1M + merged2_m$M2M + merged2_m$M3M)/3

cor(EM_Vel_M_all_m, EM_Holl_M_all_m)
cor.test(EM_Vel_M_all_m, EM_Holl_M_all_m)$p.value

paste_M_all_m = paste("cor =", round(cor(EM_Vel_M_all_m, EM_Holl_M_all_m),digits=2), ", p = ", round(cor.test(EM_Vel_M_all_m, EM_Holl_M_all_m)$p.value, digits=2), sep=" ")

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_E_all_m, EM_Holl_E_all_m, xlim=c(0,70000), ylim=c(0,70000), xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="All genes E (males)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_E_all_m, cex=0.8)
plot(EM_Vel_E_DE_m, EM_Holl_E_DE_m, xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="DE E genes (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_E_DE_m, cex=0.8)

plot(EM_Vel_M_all_m, EM_Holl_M_all_m, xlim=c(0,70000), ylim=c(0,70000), xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="All M genes (males)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_M_all_m, cex=0.8)
plot(EM_Vel_M_DE_m, EM_Holl_M_DE_m, xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="DE M genes (males)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_M_DE_m, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VH_dpse_vs_dmel_males.pdf', sep="")), width=12, height=12)
dev.off()



### males E DE only
EM_Vel_E_DE_f <- (merged3_f$E_F_V_H + merged3_f$E_F_V_H.1 + merged3_f$E_F_V_H.2 + merged3_f$E_F_V_H.3)/4
EM_Holl_E_DE_f <- (merged3_f$E1M + merged3_f$E2M + merged3_f$E3M)/3

cor(EM_Vel_E_DE_f, EM_Holl_E_DE_f)
cor.test(EM_Vel_E_DE_f, EM_Holl_E_DE_f)$p.value

paste_E_DE_f = paste("cor =", round(cor(EM_Vel_E_DE_f, EM_Holl_E_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_E_DE_f, EM_Holl_E_DE_f)$p.value, digits=2), sep=" ")

EM_Vel_E_all_f <- (merged2_f$E_F_V_H + merged2_f$E_F_V_H.1 + merged2_f$E_F_V_H.2 + merged2_f$E_F_V_H.3)/4
EM_Holl_E_all_f <- (merged2_f$E1M + merged2_f$E2M + merged2_f$E3M)/3

cor(EM_Vel_E_all_f, EM_Holl_E_all_f)
cor.test(EM_Vel_E_all_f, EM_Holl_E_all_f)$p.value

paste_E_all_f = paste("cor =", round(cor(EM_Vel_E_all_f, EM_Holl_E_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_E_all_f, EM_Holl_E_all_f)$p.value, digits=2), sep=" ")

### males M DE only
EM_Vel_f_DE_f <- (merged3_f$M_F_V_H + merged3_f$M_F_V_H.1 + merged3_f$M_F_V_H.2 + merged3_f$M_F_V_H.3)/4
EM_Holl_f_DE_f <- (merged3_f$M1M + merged3_f$M2M + merged3_f$M3M)/3

cor(EM_Vel_f_DE_f, EM_Holl_f_DE_f)
cor.test(EM_Vel_f_DE_f, EM_Holl_f_DE_f)$p.value

paste_f_DE_f = paste("cor =", round(cor(EM_Vel_f_DE_f, EM_Holl_f_DE_f),digits=2), ", p = ", round(cor.test(EM_Vel_f_DE_f, EM_Holl_f_DE_f)$p.value, digits=2), sep=" ")

EM_Vel_M_all_f <- (merged2_f$M_F_V_H + merged2_f$M_F_V_H.1 + merged2_f$M_F_V_H.2 + merged2_f$M_F_V_H.3)/4
EM_Holl_f_all_f <- (merged2_f$M1M + merged2_f$M2M + merged2_f$M3M)/3

cor(EM_Vel_M_all_f, EM_Holl_f_all_f)
cor.test(EM_Vel_M_all_f, EM_Holl_f_all_f)$p.value

paste_f_all_f = paste("cor =", round(cor(EM_Vel_M_all_f, EM_Holl_f_all_f),digits=2), ", p = ", round(cor.test(EM_Vel_M_all_f, EM_Holl_f_all_f)$p.value, digits=2), sep=" ")

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(EM_Vel_E_all_f, EM_Holl_E_all_f, xlim=c(0,70000), ylim=c(0,70000), xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="All E genes (females)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_E_all_f, cex=0.8)
plot(EM_Vel_E_DE_f, EM_Holl_E_DE_f, xlab="Mean E (Dpse)", ylab="Mean E (Dmel)", main="DE E genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_E_DE_f, cex=0.8)

plot(EM_Vel_M_all_f, EM_Holl_f_all_f, xlim=c(0,70000), ylim=c(0,70000), xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="All M genes (females)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_f_all_f, cex=0.8)
plot(EM_Vel_f_DE_f, EM_Holl_f_DE_f, xlab="Mean M (Dpse)", ylab="Mean M (Dmel)", main="DE M genes (females)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_f_DE_f, cex=0.8)

dev.copy(pdf,file.path(outpath, paste('VH_dpse_vs_dmel_females.pdf', sep="")), width=12, height=12)
dev.off()






# Virgin Heads Males vs Females

EMMVH <- read.table(file.path(subsetpath, "EMMVH/EMMVH_pseudo_counts.txt"), header=T)
EMFVH <- read.table(file.path(subsetpath, "EMFVH/EMFVH_pseudo_counts.txt"), header=T)
str(EMMVH)
str(EMFVH)

### for females
EMFVH$EminM <- (EMFVH$E_F_V_H + EMFVH$E_F_V_H.1 + EMFVH$E_F_V_H.2 + EMFVH$E_F_V_H.3)/4 - (EMFVH$M_F_V_H + EMFVH$M_F_V_H.1 + EMFVH$M_F_V_H.2 + EMFVH$M_F_V_H.3)/4
EMFVH$E <- (EMFVH$E_F_V_H + EMFVH$E_F_V_H.1 + EMFVH$E_F_V_H.2 + EMFVH$E_F_V_H.3)/4
EMFVH$M <- (EMFVH$M_F_V_H + EMFVH$M_F_V_H.1 + EMFVH$M_F_V_H.2 + EMFVH$M_F_V_H.3)/4

EMMVH$EminM <- (EMMVH$E_M_V_H + EMMVH$E_M_V_H.1 + EMMVH$E_M_V_H.2 + EMMVH$E_M_V_H.3)/4 - (EMMVH$M_M_V_H + EMMVH$M_M_V_H.1 + EMMVH$M_M_V_H.2 + EMMVH$M_M_V_H.3)/4
EMMVH$E <- (EMMVH$E_M_V_H + EMMVH$E_M_V_H.1 + EMMVH$E_M_V_H.2 + EMMVH$E_M_V_H.3)/4
EMMVH$M <- (EMMVH$M_M_V_H + EMMVH$M_M_V_H.1 + EMMVH$M_M_V_H.2 + EMMVH$M_M_V_H.3)/4

mergedVH <- merge(EMMVH, EMFVH, by.x='gene', by.y='gene', all=F)

cor(mergedVH$EminM.x, mergedVH$EminM.y)
cor.test(mergedVH$EminM.x, mergedVH$EminM.y)$p.value

paste_all_mergedVH_EminM = paste("cor =", round(cor(mergedVH$EminM.x, mergedVH$EminM.y),digits=2), ", p = ", round(cor.test(mergedVH$EminM.x, mergedVH$EminM.y)$p.value, digits=2), sep=" ")
paste_all_mergedVH_E = paste("cor =", round(cor(mergedVH$E.x, mergedVH$E.y),digits=2), ", p = ", round(cor.test(mergedVH$E.x, mergedVH$E.y)$p.value, digits=2), sep=" ")
paste_all_mergedVH_M = paste("cor =", round(cor(mergedVH$M.x, mergedVH$M.y),digits=2), ", p = ", round(cor.test(mergedVH$M.x, mergedVH$M.y)$p.value, digits=2), sep=" ")

paste_DE_mergedVH_EminM = paste("cor =", round(cor(mergedVH$EminM.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$EminM.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVH$EminM.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$EminM.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedVH_E = paste("cor =", round(cor(mergedVH$E.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$E.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVH$E.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$E.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedVH_M = paste("cor =", round(cor(mergedVH$M.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$M.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVH$M.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$M.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1])$p.value, digits=2), sep=" ")

# Virgin Bodies

EMMVB <- read.table(file.path(subsetpath, "EMMVB/EMMVB_pseudo_counts.txt"), header=T)
EMFVB <- read.table(file.path(subsetpath, "EMFVB/EMFVB_pseudo_counts.txt"), header=T)
str(EMMVB)
str(EMFVB)

### for females
EMFVB$EminM <- (EMFVB$E_F_V_B + EMFVB$E_F_V_B.1 + EMFVB$E_F_V_B.2 + EMFVB$E_F_V_B.3)/4 - (EMFVB$M_F_V_B + EMFVB$M_F_V_B.1 + EMFVB$M_F_V_B.2 + EMFVB$M_F_V_B.3)/4
EMFVB$E <- (EMFVB$E_F_V_B + EMFVB$E_F_V_B.1 + EMFVB$E_F_V_B.2 + EMFVB$E_F_V_B.3)/4
EMFVB$M <- (EMFVB$M_F_V_B + EMFVB$M_F_V_B.1 + EMFVB$M_F_V_B.2 + EMFVB$M_F_V_B.3)/4

EMMVB$EminM <- (EMMVB$E_M_V_B + EMMVB$E_M_V_B.1 + EMMVB$E_M_V_B.2 + EMMVB$E_M_V_B.3)/4 - (EMMVB$M_M_V_B + EMMVB$M_M_V_B.1 + EMMVB$M_M_V_B.2 + EMMVB$M_M_V_B.3)/4
EMMVB$E <- (EMMVB$E_M_V_B + EMMVB$E_M_V_B.1 + EMMVB$E_M_V_B.2 + EMMVB$E_M_V_B.3)/4
EMMVB$M <- (EMMVB$M_M_V_B + EMMVB$M_M_V_B.1 + EMMVB$M_M_V_B.2 + EMMVB$M_M_V_B.3)/4

mergedVB <- merge(EMMVB, EMFVB, by.x='gene', by.y='gene', all=F)

cor(mergedVB$EminM.x, mergedVB$EminM.y)
cor.test(mergedVB$EminM.x, mergedVB$EminM.y)$p.value

paste_all_mergedVB_EminM = paste("cor =", round(cor(mergedVB$EminM.x, mergedVB$EminM.y),digits=2), ", p = ", round(cor.test(mergedVB$EminM.x, mergedVB$EminM.y)$p.value, digits=2), sep=" ")
paste_all_mergedVB_E = paste("cor =", round(cor(mergedVB$E.x, mergedVB$E.y),digits=2), ", p = ", round(cor.test(mergedVB$E.x, mergedVB$E.y)$p.value, digits=2), sep=" ")
paste_all_mergedVB_M = paste("cor =", round(cor(mergedVB$M.x, mergedVB$M.y),digits=2), ", p = ", round(cor.test(mergedVB$M.x, mergedVB$M.y)$p.value, digits=2), sep=" ")

paste_DE_mergedVB_EminM = paste("cor =", round(cor(mergedVB$EminM.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$EminM.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVB$EminM.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$EminM.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedVB_E = paste("cor =", round(cor(mergedVB$E.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$E.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVB$E.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$E.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedVB_M = paste("cor =", round(cor(mergedVB$M.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$M.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedVB$M.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$M.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1])$p.value, digits=2), sep=" ")

# Courted Heads

EMMCH <- read.table(file.path(subsetpath, "EMMCH/EMMCH_pseudo_counts.txt"), header=T)
EMFCH <- read.table(file.path(subsetpath, "EMFCH/EMFCH_pseudo_counts.txt"), header=T)
str(EMMCH)
str(EMFCH)

### for females
EMFCH$EminM <- (EMFCH$E_F_C_H + EMFCH$E_F_C_H.1 + EMFCH$E_F_C_H.2 + EMFCH$E_F_C_H.3)/4 - (EMFCH$M_F_C_H + EMFCH$M_F_C_H.1 + EMFCH$M_F_C_H.2 + EMFCH$M_F_C_H.3)/4
EMFCH$E <- (EMFCH$E_F_C_H + EMFCH$E_F_C_H.1 + EMFCH$E_F_C_H.2 + EMFCH$E_F_C_H.3)/4
EMFCH$M <- (EMFCH$M_F_C_H + EMFCH$M_F_C_H.1 + EMFCH$M_F_C_H.2 + EMFCH$M_F_C_H.3)/4

EMMCH$EminM <- (EMMCH$E_M_C_H + EMMCH$E_M_C_H.1 + EMMCH$E_M_C_H.2 + EMMCH$E_M_C_H.3)/4 - (EMMCH$M_M_C_H + EMMCH$M_M_C_H.1 + EMMCH$M_M_C_H.2 + EMMCH$M_M_C_H.3)/4
EMMCH$E <- (EMMCH$E_M_C_H + EMMCH$E_M_C_H.1 + EMMCH$E_M_C_H.2 + EMMCH$E_M_C_H.3)/4
EMMCH$M <- (EMMCH$M_M_C_H + EMMCH$M_M_C_H.1 + EMMCH$M_M_C_H.2 + EMMCH$M_M_C_H.3)/4

mergedCH <- merge(EMMCH, EMFCH, by.x='gene', by.y='gene', all=F)

cor(mergedCH$EminM.x, mergedCH$EminM.y)
cor.test(mergedCH$EminM.x, mergedCH$EminM.y)$p.value

paste_all_mergedCH_EminM = paste("cor =", round(cor(mergedCH$EminM.x, mergedCH$EminM.y),digits=2), ", p = ", round(cor.test(mergedCH$EminM.x, mergedCH$EminM.y)$p.value, digits=2), sep=" ")
paste_all_mergedCH_E = paste("cor =", round(cor(mergedCH$E.x, mergedCH$E.y),digits=2), ", p = ", round(cor.test(mergedCH$E.x, mergedCH$E.y)$p.value, digits=2), sep=" ")
paste_all_mergedCH_M = paste("cor =", round(cor(mergedCH$M.x, mergedCH$M.y),digits=2), ", p = ", round(cor.test(mergedCH$M.x, mergedCH$M.y)$p.value, digits=2), sep=" ")

paste_DE_mergedCH_EminM = paste("cor =", round(cor(mergedCH$EminM.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$EminM.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCH$EminM.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$EminM.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedCH_E = paste("cor =", round(cor(mergedCH$E.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$E.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCH$E.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$E.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedCH_M = paste("cor =", round(cor(mergedCH$M.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$M.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCH$M.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$M.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1])$p.value, digits=2), sep=" ")


# Courted Bodies

EMMCB <- read.table(file.path(subsetpath, "EMMCB/EMMCB_pseudo_counts.txt"), header=T)
EMFCB <- read.table(file.path(subsetpath, "EMFCB/EMFCB_pseudo_counts.txt"), header=T)
str(EMMCB)
str(EMFCB)

### for females
EMFCB$EminM <- (EMFCB$E_F_C_B + EMFCB$E_F_C_B.1 + EMFCB$E_F_C_B.2 + EMFCB$E_F_C_B.3)/4 - (EMFCB$M_F_C_B + EMFCB$M_F_C_B.1 + EMFCB$M_F_C_B.2 + EMFCB$M_F_C_B.3)/4
EMFCB$E <- (EMFCB$E_F_C_B + EMFCB$E_F_C_B.1 + EMFCB$E_F_C_B.2 + EMFCB$E_F_C_B.3)/4
EMFCB$M <- (EMFCB$M_F_C_B + EMFCB$M_F_C_B.1 + EMFCB$M_F_C_B.2 + EMFCB$M_F_C_B.3)/4

EMMCB$EminM <- (EMMCB$E_M_C_B + EMMCB$E_M_C_B.1 + EMMCB$E_M_C_B.2 + EMMCB$E_M_C_B.3)/4 - (EMMCB$M_M_C_B + EMMCB$M_M_C_B.1 + EMMCB$M_M_C_B.2 + EMMCB$M_M_C_B.3)/4
EMMCB$E <- (EMMCB$E_M_C_B + EMMCB$E_M_C_B.1 + EMMCB$E_M_C_B.2 + EMMCB$E_M_C_B.3)/4
EMMCB$M <- (EMMCB$M_M_C_B + EMMCB$M_M_C_B.1 + EMMCB$M_M_C_B.2 + EMMCB$M_M_C_B.3)/4

mergedCB <- merge(EMMCB, EMFCB, by.x='gene', by.y='gene', all=F)

cor(mergedCB$EminM.x, mergedCB$EminM.y)
cor.test(mergedCB$EminM.x, mergedCB$EminM.y)$p.value

paste_all_mergedCB_EminM = paste("cor =", round(cor(mergedCB$EminM.x, mergedCB$EminM.y),digits=2), ", p = ", round(cor.test(mergedCB$EminM.x, mergedCB$EminM.y)$p.value, digits=2), sep=" ")
paste_all_mergedCB_E = paste("cor =", round(cor(mergedCB$E.x, mergedCB$E.y),digits=2), ", p = ", round(cor.test(mergedCB$E.x, mergedCB$E.y)$p.value, digits=2), sep=" ")
paste_all_mergedCB_M = paste("cor =", round(cor(mergedCB$M.x, mergedCB$M.y),digits=2), ", p = ", round(cor.test(mergedCB$M.x, mergedCB$M.y)$p.value, digits=2), sep=" ")

paste_DE_mergedCB_EminM = paste("cor =", round(cor(mergedCB$EminM.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$EminM.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCB$EminM.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$EminM.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedCB_E = paste("cor =", round(cor(mergedCB$E.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$E.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCB$E.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$E.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1])$p.value, digits=2), sep=" ")
paste_DE_mergedCB_M = paste("cor =", round(cor(mergedCB$M.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$M.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1]),digits=2), ", p = ", round(cor.test(mergedCB$M.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$M.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1])$p.value, digits=2), sep=" ")


## linear models
library(sjPlot)
library(lmodel2)
mergedVB_DE <- subset(mergedVB, mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1) 
mergedCB_DE <- subset(mergedCB, mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1) 
mergedVH_DE <- subset(mergedVH, mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1) 
mergedCH_DE <- subset(mergedCH, mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1) 

VB_male <- data.frame(c(mergedVB_DE$E.x, mergedVB_DE$M.x))
VB_male$sline <- 'M'
for (i in 1:(nrow(VB_male)/2)) VB_male$sline[i] <- 'E'
VB_female <- data.frame(c(mergedVB_DE$E.y, mergedVB_DE$M.y))
VB_frame <- data.frame(VB_male, VB_female)
colnames(VB_frame) <- c('male', 'sline', 'female')
VB_frame$sline <- factor(VB_frame$sline)
mod_VB <- lm(VB_frame$male ~ VB_frame$female * VB_frame$sline)
summary(mod_VB)
anova(mod_VB)
sjt.lm(mod_VB)


CB_male <- data.frame(c(mergedCB_DE$E.x, mergedCB_DE$M.x))
CB_male$sline <- 'M'
for (i in 1:(nrow(CB_male)/2)) CB_male$sline[i] <- 'E'
CB_female <- data.frame(c(mergedCB_DE$E.y, mergedCB_DE$M.y))
CB_frame <- data.frame(CB_male, CB_female)
colnames(CB_frame) <- c('male', 'sline', 'female')
CB_frame$sline <- factor(CB_frame$sline)
mod_CB <- lm(CB_frame$male ~ CB_frame$female * CB_frame$sline)
summary(mod_CB)
anova(mod_CB)
sjt.lm(mod_CB)

VH_male <- data.frame(c(mergedVH_DE$E.x, mergedVH_DE$M.x))
VH_male$sline <- 'M'
for (i in 1:(nrow(VH_male)/2)) VH_male$sline[i] <- 'E'
VH_female <- data.frame(c(mergedVH_DE$E.y, mergedVH_DE$M.y))
VH_frame <- data.frame(VH_male, VH_female)
colnames(VH_frame) <- c('male', 'sline', 'female')
VH_frame$sline <- factor(VH_frame$sline)
mod_VH <- lm(VH_frame$male ~ VH_frame$female * VH_frame$sline)
summary(mod_VH)
anova(mod_VH)
sjt.lm(mod_VH)

CH_male <- data.frame(c(mergedCH_DE$E.x, mergedCH_DE$M.x))
CH_male$sline <- 'M'
for (i in 1:(nrow(CH_male)/2)) CH_male$sline[i] <- 'E'
CH_female <- data.frame(c(mergedCH_DE$E.y, mergedCH_DE$M.y))
CH_frame <- data.frame(CH_male, CH_female)
colnames(CH_frame) <- c('male', 'sline', 'female')
CH_frame$sline <- factor(CH_frame$sline)
mod_CH <- lm(CH_frame$male ~ CH_frame$female * CH_frame$sline)
summary(mod_CH)
anova(mod_CH)
sjt.lm(mod_CH)

sjt.lm(mod_VB, mod_CB, mod_VH, mod_CH, emph.p=T, p.numeric=T, separate.ci.col=F, show.header=T, show.aic=T, file=file.path(outpath, 'lm_DE.html'), sep.column=F)


# lmodel2 (does not do interactions)

mod2 <- lmodel2(male ~ female + sline, data=CH_frame, nperm=99)
plot(mod2, "OLS")

# log ratios

head(VB_frame)

VB_frame$logratio <- log2((1 + VB_frame$male) / (1 + VB_frame$female))
mod_VB <- lm(logratio ~ sline, data=VB_frame)

summary(mod_VB)

anova(mod_VB)

VH_frame$logratio <- log2((1 + VH_frame$male) / (1 + VH_frame$female))
mod_VH <- lm(logratio ~ sline, data=VH_frame)
summary(mod_VH)
anova(mod_VH)

CB_frame$logratio <- log2((1 + CB_frame$male) / (1 + CB_frame$female))
mod_CB <- lm(logratio ~ sline, data=CB_frame)
summary(mod_CB)
anova(mod_CB)

CH_frame$logratio <- log2((1 + CH_frame$male) / (1 + CH_frame$female))
mod_CH <- lm(logratio ~ sline, data=CH_frame)
summary(mod_CH)
anova(mod_CH)


boxplot(VB_frame$logratio[VB_frame$sline=='E'], VB_frame$logratio[VB_frame$sline=='M'],
		VH_frame$logratio[VH_frame$sline=='E'], VH_frame$logratio[VH_frame$sline=='M'],
		CB_frame$logratio[CB_frame$sline=='E'], CB_frame$logratio[CB_frame$sline=='M'],
		CH_frame$logratio[CH_frame$sline=='E'], CH_frame$logratio[CH_frame$sline=='M'],
		names=list("E-VB", "M-VB","E-VH", "M-VH","E-CB", "M-CB","E-CH", "M-CH"), col=c(2,4), main ='Boxplots')

sjt.lm(mod_VB, mod_CB, mod_VH, mod_CH, emph.p=T, p.numeric=T, separate.ci.col=F, show.header=T, show.aic=T, file=file.path(outpath, 'lm_DE.html'), sep.column=F)


# Export files for Mike
write.table(CH_frame, file=file.path(outpath, "CH_out.txt"), quote=F, row.names=F, sep="\t") 
write.table(CB_frame, file=file.path(outpath, "CB_out.txt"), quote=F, row.names=F, sep="\t") 
write.table(VH_frame, file=file.path(outpath, "VH_out.txt"), quote=F, row.names=F, sep="\t") 
write.table(VB_frame, file=file.path(outpath, "VB_out.txt"), quote=F, row.names=F, sep="\t") 


par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$EminM.x, mergedVB$EminM.y, xlim=c(-15000,10000), ylim=c(-15000,10000), xlab="Males", ylab="Females", main="All Virgin Abdomens (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedVB_EminM, cex=0.8)
plot(mergedVH$EminM.x, mergedVH$EminM.y, xlim=c(-40000,40000), ylim=c(-40000,40000), xlab="Males", ylab="Females", main="All Virgin Heads (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedVH_EminM, cex=0.8)
plot(mergedCB$EminM.x, mergedCB$EminM.y, xlim=c(-20000,20000), ylim=c(-20000,20000), xlab="Males", ylab="Females", main="All Courted Abdomens (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedCB_EminM, cex=0.8)
plot(mergedCH$EminM.x, mergedCH$EminM.y, xlim=c(-20000,20000), ylim=c(-20000,20000), xlab="Males", ylab="Females", main="All Courted Heads (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedCH_EminM, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('All_males_females_EminM.pdf', sep="")), width=12, height=12)
dev.off()

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$EminM.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$EminM.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Abdomens (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedVB_EminM, cex=0.8)
plot(mergedVH$EminM.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$EminM.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Heads (E-M)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_DE_mergedVH_EminM, cex=0.8)
plot(mergedCB$EminM.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$EminM.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Abdomens (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedCB_EminM, cex=0.8)
plot(mergedCH$EminM.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$EminM.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Heads (E-M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedCH_EminM, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('DE_males_females_EminM.pdf', sep="")), width=12, height=12)
dev.off()

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$E.x, mergedVB$E.y,xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Virgin Abdomens (E)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_all_mergedVB_E, cex=0.8)
plot(mergedVH$E.x, mergedVH$E.y, xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Virgin Heads (E)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedVH_E, cex=0.8)
plot(mergedCB$E.x, mergedCB$E.y, xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Courted Abdomens (E)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_all_mergedCB_E, cex=0.8)
plot(mergedCH$E.x, mergedCH$E.y, xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Courted Heads (E)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedCH_E, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('All_males_females_E.pdf', sep="")), width=12, height=12)
dev.off()

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$E.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$E.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Abdomens (E)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedVB_E, cex=0.8)
plot(mergedVH$E.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$E.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Heads (E)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_DE_mergedVH_E, cex=0.8)
plot(mergedCB$E.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$E.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Abdomens (E)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedCB_E, cex=0.8)
plot(mergedCH$E.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$E.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Heads (E)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedCH_E, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('DE_males_females_E.pdf', sep="")), width=12, height=12)
dev.off()



par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$M.x, mergedVB$M.y,xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Virgin Abdomens (M)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_all_mergedVB_M, cex=0.8)
plot(mergedVH$M.x, mergedVH$M.y,xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Virgin Heads (M)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_all_mergedVH_M, cex=0.8)
plot(mergedCB$M.x, mergedCB$M.y, xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Courted Abdomens (M)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_all_mergedCB_M, cex=0.8)
plot(mergedCH$M.x, mergedCH$M.y, xlim=c(0,100000), ylim=c(0,100000), xlab="Males", ylab="Females", main="All Courted Heads (M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_all_mergedCH_M, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('All_males_females_M.pdf', sep="")), width=12, height=12)
dev.off()

par(mfrow=c(2,2)) 
par(mar=c(5,5,4,3))
plot(mergedVB$M.x[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], mergedVB$M.y[mergedVB$FDR.x < 0.1 | mergedVB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Abdomens (M)", cex.main=1.5, cex.lab=1.2)
legend('topright', inset=0.05, legend=paste_DE_mergedVB_M, cex=0.8)
plot(mergedVH$M.x[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], mergedVH$M.y[mergedVH$FDR.x < 0.1 | mergedVH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Virgin Heads (M)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_DE_mergedVH_M, cex=0.8)
plot(mergedCB$M.x[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], mergedCB$M.y[mergedCB$FDR.x < 0.1 | mergedCB$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Abdomens (M)", cex.main=1.5, cex.lab=1.2)
legend('topleft', inset=0.05, legend=paste_DE_mergedCB_M, cex=0.8)
plot(mergedCH$M.x[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], mergedCH$M.y[mergedCH$FDR.x < 0.1 | mergedCH$FDR.y < 0.1], xlab="Males", ylab="Females", main="DE Courted Heads (M)", cex.main=1.5, cex.lab=1.2)
legend('bottomright', inset=0.05, legend=paste_DE_mergedCH_M, cex=0.8)
dev.copy(pdf,file.path(outpath, paste('DE_males_females_M.pdf', sep="")), width=12, height=12)
dev.off()



