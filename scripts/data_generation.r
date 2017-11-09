library(edgeR)

cbgreen <- '#009E73'
cblblue <- '#56B4E9'

datapath <- "~/git/feminisation_direction/input"
load(file.path(datapath, "designtab.Rdata")) # table of design information 
count <- read.table(file.path(datapath, 'new_count.txt'), header=T)
outpath <- '~/git/feminisation_direction/output/all'
dir.create(file.path(outpath))
nrow(count)
	
annotation <- read.delim(file.path(datapath, "exongeneanno.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
model.formula <- as.formula("~0+group")
design <- data.frame(designtab)
dmat <- model.matrix(model.formula, data=design)
dgl <- DGEList(counts=count, group=design$group, genes=annotation)

dgl <- dgl[rowSums(cpm(dgl)) >= 10,] # keep cpm > 2 for further analysis

y <- dgl
colnames(y) <- paste(colnames(y), design$group, sep="\n")


## MDS
pdf(file.path(outpath, 'MDS_all.pdf'), width=8, height=8)
par(mar=c(5,5,4,3))
plotMDS(y, cex=0.4, col=as.numeric(y$samples$group), main="all_MDS plot")
dev.off()

xcpm <- mglmOneGroup(dgl$counts) # computing a logCPM for making dispersion plot
dgl <- calcNormFactors(dgl)
dgl <- estimateGLMCommonDisp(dgl, dmat)
dgl <- estimateGLMTrendedDisp(dgl, dmat, min.n=1000)
dgl <- estimateGLMTagwiseDisp(dgl, dmat)

## dispersion plot
pdf(file.path(outpath, paste('Dispersion.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste("dispersion", sep=""))
if(!is.null(dgl$trended.dispersion)) points(xcpm ,dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=dgl$common.dispersion, col=cblblue, lwd=2)
legend("topright", c("Common","Trended","Tagwise"), pch=16, col=c(cblblue, cbgreen, "black"), title="Dispersion")
dev.off()

## export results and normalised counts
nc <- cpm(dgl, prior.count=2, log=T)
nc2 <- data.frame(row.names(nc), nc)
colnames(nc2) <- c('ID', paste(design$group, design$line, sep="_" ))
 # restab_logCPM = merge(restab, nc2, by.x="gene", by.y="ID", all=F )
write.table(nc2, file=file.path(outpath, paste('normCounts.txt', sep="")), quote=F, row.names=F, sep='\t')


# courtship analysis input design and counts generation for EMB, B, E, M contrasts
EMB_MH <- c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93)
EMB_MB <- EMB_MH+1
EMB_FH <- EMB_MH+2
EMB_FB <- EMB_MH+3

EMB_MH_count <- data.frame(count[,EMB_MH])
EMB_MB_count <- data.frame(count[,EMB_MB])
EMB_FH_count <- data.frame(count[,EMB_FH])
EMB_FB_count <- data.frame(count[,EMB_FB])

EMB_MH_design <- designtab[EMB_MH,]
EMB_MB_design <- designtab[EMB_MB,]
EMB_FH_design <- designtab[EMB_FH,]
EMB_FB_design <- designtab[EMB_FB,]

write.table(EMB_MH_count, file=file.path(datapath, 'EMB_MH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_MH_design, file=file.path(datapath, 'EMB_MH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_MB_count, file=file.path(datapath, 'EMB_MB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_MB_design, file=file.path(datapath, 'EMB_MB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_FH_count, file=file.path(datapath, 'EMB_FH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_FH_design, file=file.path(datapath, 'EMB_FH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_FB_count, file=file.path(datapath, 'EMB_FB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_FB_design, file=file.path(datapath, 'EMB_FB_design.txt'), quote=F, row.names=F, sep='\t')


B_MH <- c(1,5,9,13,17,21,25,29)
B_MB <- B_MH+1
B_FH <- B_MH+2
B_FB <- B_MH+3

B_MH_count <- data.frame(count[,B_MH])
B_MB_count <- data.frame(count[,B_MB])
B_FH_count <- data.frame(count[,B_FH])
B_FB_count <- data.frame(count[,B_FB])

B_MH_design <- designtab[B_MH,]
B_MB_design <- designtab[B_MB,]
B_FH_design <- designtab[B_FH,]
B_FB_design <- designtab[B_FB,]

write.table(B_MH_count, file=file.path(datapath, 'B_MH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_MH_design, file=file.path(datapath, 'B_MH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_MB_count, file=file.path(datapath, 'B_MB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_MB_design, file=file.path(datapath, 'B_MB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_FH_count, file=file.path(datapath, 'B_FH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_FH_design, file=file.path(datapath, 'B_FH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_FB_count, file=file.path(datapath, 'B_FB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_FB_design, file=file.path(datapath, 'B_FB_design.txt'), quote=F, row.names=F, sep='\t')

M_MH <- B_MH + 32
M_MB <- B_MB + 32
M_FH <- B_FH + 32
M_FB <- B_FB + 32

M_MH_count <- data.frame(count[,M_MH])
M_MB_count <- data.frame(count[,M_MB])
M_FH_count <- data.frame(count[,M_FH])
M_FB_count <- data.frame(count[,M_FB])

M_MH_design <- designtab[M_MH,]
M_MB_design <- designtab[M_MB,]
M_FH_design <- designtab[M_FH,]
M_FB_design <- designtab[M_FB,]

write.table(M_MH_count, file=file.path(datapath, 'M_MH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_MH_design, file=file.path(datapath, 'M_MH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_MB_count, file=file.path(datapath, 'M_MB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_MB_design, file=file.path(datapath, 'M_MB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_FH_count, file=file.path(datapath, 'M_FH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_FH_design, file=file.path(datapath, 'M_FH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_FB_count, file=file.path(datapath, 'M_FB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_FB_design, file=file.path(datapath, 'M_FB_design.txt'), quote=F, row.names=F, sep='\t')

E_MH <- B_MH + 64
E_MB <- B_MB + 64
E_FH <- B_FH + 64
E_FB <- B_FB + 64

E_MH_count <- data.frame(count[,E_MH])
E_MB_count <- data.frame(count[,E_MB])
E_FH_count <- data.frame(count[,E_FH])
E_FB_count <- data.frame(count[,E_FB])

E_MH_design <- designtab[E_MH,]
E_MB_design <- designtab[E_MB,]
E_FH_design <- designtab[E_FH,]
E_FB_design <- designtab[E_FB,]


write.table(E_MH_count, file=file.path(datapath, 'E_MH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_MH_design, file=file.path(datapath, 'E_MH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_MB_count, file=file.path(datapath, 'E_MB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_MB_design, file=file.path(datapath, 'E_MB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_FH_count, file=file.path(datapath, 'E_FH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_FH_design, file=file.path(datapath, 'E_FH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_FB_count, file=file.path(datapath, 'E_FB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_FB_design, file=file.path(datapath, 'E_FB_design.txt'), quote=F, row.names=F, sep='\t')


# sex analysis input design and counts generation for EMB, B, E, M contrasts
EMB_CH <- c(1,3,5,7,9,11,13,15,33,35,37,39,41,43,45,47,65,67,69,71,73,75,77,79)
EMB_CB <- c(2,4,6,8,10,12,14,16,34,36,38,40,42,44,46,48,66,68,70,72,74,76,78,80)
EMB_VH <- c(17,19,21,23,25,27,29,31,49,51,53,55,57,59,61,63,81,83,85,87,89,91,93,95)
EMB_VB <- c(18,20,22,24,26,28,30,32,50,52,54,56,58,60,62,64,82,84,86,88,90,92,94,96)

EMB_CH_count <- data.frame(count[,EMB_CH])
EMB_CB_count <- data.frame(count[,EMB_CB])
EMB_VH_count <- data.frame(count[,EMB_VH])
EMB_VB_count <- data.frame(count[,EMB_VB])

EMB_CH_design <- designtab[EMB_CH,]
EMB_CB_design <- designtab[EMB_CB,]
EMB_VH_design <- designtab[EMB_VH,]
EMB_VB_design <- designtab[EMB_VB,]


write.table(EMB_CH_count, file=file.path(datapath, 'EMB_CH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_CH_design, file=file.path(datapath, 'EMB_CH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_CB_count, file=file.path(datapath, 'EMB_CB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_CB_design, file=file.path(datapath, 'EMB_CB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_VH_count, file=file.path(datapath, 'EMB_VH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_VH_design, file=file.path(datapath, 'EMB_VH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_VB_count, file=file.path(datapath, 'EMB_VB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_VB_design, file=file.path(datapath, 'EMB_VB_design.txt'), quote=F, row.names=F, sep='\t')


B_CH <- c(1,3,5,7,9,11,13,15)
B_CB <- c(2,4,6,8,10,12,14,16)
B_VH <- c(17,19,21,23,25,27,29,31)
B_VB <- c(18,20,22,24,26,28,30,32)

B_CH_count <- data.frame(count[,B_CH])
B_CB_count <- data.frame(count[,B_CB])
B_VH_count <- data.frame(count[,B_VH])
B_VB_count <- data.frame(count[,B_VB])

B_CH_design <- designtab[B_CH,]
B_CB_design <- designtab[B_CB,]
B_VH_design <- designtab[B_VH,]
B_VB_design <- designtab[B_VB,]


write.table(B_CH_count, file=file.path(datapath, 'B_CH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_CH_design, file=file.path(datapath, 'B_CH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_CB_count, file=file.path(datapath, 'B_CB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_CB_design, file=file.path(datapath, 'B_CB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_VH_count, file=file.path(datapath, 'B_VH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_VH_design, file=file.path(datapath, 'B_VH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_VB_count, file=file.path(datapath, 'B_VB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_VB_design, file=file.path(datapath, 'B_VB_design.txt'), quote=F, row.names=F, sep='\t')


M_CH <- B_CH + 32
M_CB <- B_CB + 32
M_VH <- B_VH + 32
M_VB <- B_VB + 32

M_CH_count <- data.frame(count[,M_CH])
M_CB_count <- data.frame(count[,M_CB])
M_VH_count <- data.frame(count[,M_VH])
M_VB_count <- data.frame(count[,M_VB])

M_CH_design <- designtab[M_CH,]
M_CB_design <- designtab[M_CB,]
M_VH_design <- designtab[M_VH,]
M_VB_design <- designtab[M_VB,]

write.table(M_CH_count, file=file.path(datapath, 'M_CH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_CH_design, file=file.path(datapath, 'M_CH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_CB_count, file=file.path(datapath, 'M_CB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_CB_design, file=file.path(datapath, 'M_CB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_VH_count, file=file.path(datapath, 'M_VH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_VH_design, file=file.path(datapath, 'M_VH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_VB_count, file=file.path(datapath, 'M_VB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_VB_design, file=file.path(datapath, 'M_VB_design.txt'), quote=F, row.names=F, sep='\t')


E_CH <- B_CH + 64
E_CB <- B_CB + 64
E_VH <- B_VH + 64
E_VB <- B_VB + 64

E_CH_count <- data.frame(count[,E_CH])
E_CB_count <- data.frame(count[,E_CB])
E_VH_count <- data.frame(count[,E_VH])
E_VB_count <- data.frame(count[,E_VB])

E_CH_design <- designtab[E_CH,]
E_CB_design <- designtab[E_CB,]
E_VH_design <- designtab[E_VH,]
E_VB_design <- designtab[E_VB,]

write.table(E_CH_count, file=file.path(datapath, 'E_CH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_CH_design, file=file.path(datapath, 'E_CH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_CB_count, file=file.path(datapath, 'E_CB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_CB_design, file=file.path(datapath, 'E_CB_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_VH_count, file=file.path(datapath, 'E_VH_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_VH_design, file=file.path(datapath, 'E_VH_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_VB_count, file=file.path(datapath, 'E_VB_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_VB_design, file=file.path(datapath, 'E_VB_design.txt'), quote=F, row.names=F, sep='\t')


# sex analysis input design and counts generation for EMB, B, E, M contrasts - combined by courtship status
EMB_H <- c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95)
EMB_B <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)

EMB_H_count <- data.frame(count[,EMB_H])
EMB_B_count <- data.frame(count[,EMB_B])

EMB_H_design <- designtab[EMB_H,]
EMB_B_design <- designtab[EMB_B,]

write.table(EMB_H_count, file=file.path(datapath, 'EMB_H_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_H_design, file=file.path(datapath, 'EMB_H_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_B_count, file=file.path(datapath, 'EMB_B_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_B_design, file=file.path(datapath, 'EMB_B_design.txt'), quote=F, row.names=F, sep='\t')


B_H <- c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
B_B <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)

B_H_count <- data.frame(count[,B_H])
B_B_count <- data.frame(count[,B_B])

B_H_design <- designtab[B_H,]
B_B_design <- designtab[B_B,]

write.table(B_H_count, file=file.path(datapath, 'B_H_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_H_design, file=file.path(datapath, 'B_H_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_B_count, file=file.path(datapath, 'B_B_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_B_design, file=file.path(datapath, 'B_B_design.txt'), quote=F, row.names=F, sep='\t')


M_H <- B_H + 32
M_B <- B_B + 32

M_H_count <- data.frame(count[,M_H])
M_B_count <- data.frame(count[,M_B])

M_H_design <- designtab[M_H,]
M_B_design <- designtab[M_B,]

write.table(M_H_count, file=file.path(datapath, 'M_H_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_H_design, file=file.path(datapath, 'M_H_design.txt'), quote=F, row.names=F, sep='\t')
write.table(M_B_count, file=file.path(datapath, 'M_B_count.txt'), quote=F, row.names=T, sep='\t')
write.table(M_B_design, file=file.path(datapath, 'M_B_design.txt'), quote=F, row.names=F, sep='\t')

E_H <- B_H + 64
E_B <- B_B + 64

E_H_count <- data.frame(count[,E_H])
E_B_count <- data.frame(count[,E_B])

E_H_design <- designtab[E_H,]
E_B_design <- designtab[E_B,]

write.table(E_H_count, file=file.path(datapath, 'E_H_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_H_design, file=file.path(datapath, 'E_H_design.txt'), quote=F, row.names=F, sep='\t')
write.table(E_B_count, file=file.path(datapath, 'E_B_count.txt'), quote=F, row.names=T, sep='\t')
write.table(E_B_design, file=file.path(datapath, 'E_B_design.txt'), quote=F, row.names=F, sep='\t')


# tissue analysis input design and counts generation for EMB, B, E, M contrasts - combined by courtship status
EMB_M <- c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,89,90,93,94)
EMB_F <- c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,87,88,91,92,95,96)

EMB_M_count <- data.frame(count[,EMB_M])
EMB_F_count <- data.frame(count[,EMB_F])

EMB_M_design <- designtab[EMB_M,]
EMB_F_design <- designtab[EMB_F,]

write.table(EMB_M_count, file=file.path(datapath, 'EMB_M_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_M_design, file=file.path(datapath, 'EMB_M_design.txt'), quote=F, row.names=F, sep='\t')
write.table(EMB_F_count, file=file.path(datapath, 'EMB_F_count.txt'), quote=F, row.names=T, sep='\t')
write.table(EMB_F_design, file=file.path(datapath, 'EMB_F_design.txt'), quote=F, row.names=F, sep='\t')


B_M <- c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30)
B_F <- c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32)

B_M_count <- data.frame(count[,B_M])
B_F_count <- data.frame(count[,B_F])

B_M_design <- designtab[B_M,]
B_F_design <- designtab[B_F,]

write.table(B_M_count, file=file.path(datapath, 'B_M_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_M_design, file=file.path(datapath, 'B_M_design.txt'), quote=F, row.names=F, sep='\t')
write.table(B_F_count, file=file.path(datapath, 'B_F_count.txt'), quote=F, row.names=T, sep='\t')
write.table(B_F_design, file=file.path(datapath, 'B_F_design.txt'), quote=F, row.names=F, sep='\t')

