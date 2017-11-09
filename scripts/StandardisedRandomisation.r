library(limma)
graphics.off()

# change to location where fileS3 is unzipped
datap="~/git/feminisation_direction/input/pseudo_counts"
setwd(datap)
dataf=list.files()[grep("EM[A-Z]",list.files())]
dataf=paste0(datap,"/",dataf)

outpath="~/git/feminisation_direction/output/results_standard_randomisation"
if(!file.exists(outpath)) dir.create(outpath)
itr=choose(4,2)*choose(4,2)

extrem=list("","FBgn0255061","","","","FBgn0255061","","")
extrem="FBgn0255061"

dataname <- substr(basename(dataf),1,5)

set.seed(1234)
p=matrix(1,length(dataf),6)
rownames(p)=dataname
colnames(p)=c("p.male.lower", "p.female.lower", "p.neutral.lower", "p.male.upper", "p.female.upper", "p.neutral.upper")

male=NULL
female=NULL
neutral=NULL

for(i in 1:length(dataf))
{
        x=read.delim(dataf[i],sep="\t",header=1)
        x=x[!(x[,1] %in% extrem),]
        # male and female bias gene index
        fidx=x[,"bias"] %in% "female"
        midx=x[,"bias"] %in% "male"
        nidx=x[,"bias"] %in% "neutral"
        # counrtship E and M column index 
        # counrtship E and M column index 
        ib.c=grep("B_[A-Z]_C",colnames(x))
        ie.c=grep("E_[A-Z]_C",colnames(x))
        im.c=grep("M_[A-Z]_C",colnames(x))
        # Virgin E and M column index 
        ib.v=grep("B_[A-Z]_V",colnames(x))
        ie.v=grep("E_[A-Z]_V",colnames(x))
        im.v=grep("M_[A-Z]_V",colnames(x))
###
	if(substr(dataname[i],4,4) == "C") x8=x[,c(im.c,ie.c)] else x8=x[,c(im.v,ie.v)] 	
	y=as.matrix(round(x8))
	rownames(y)=x[,1]
###
	comb1=combn(1:4,2)
	comb2=combn(5:8,2)
	comb=1:4
	for(j in 1:ncol(comb1)) for(k in 1:ncol(comb2)) comb=cbind(comb,c(comb1[,j],comb2[,k]))
        f=m=u=NULL 
        itr=0
	group=factor(rep("M",8),levels=c("M","E"))
	group[5:8]="E"
	sgrp=data.frame(grp=group)
	dmat=model.matrix(~grp,sgrp) 
	dmat[,2] = dmat[,2] - 0.5
	for(j in 1:(ncol(comb)-1)) {
	  if(j==1) {
#	  dmat=model.matrix(~grp,sgrp)
	  vcpm=voom(y,dmat,plot=FALSE)
	  fits=lmFit(vcpm,dmat)
	  efits=eBayes(fits)
	  pv=efits$p.value[,2,drop=FALSE]
	  lfc=efits$coefficients[,2,drop=FALSE]

	  stdlfc=sign(lfc)*abs(qnorm(pv/2+1E-16,0,1))
          f=c(f,mean(stdlfc[fidx]));
          m=c(m,mean(stdlfc[midx])); 
          u=c(u,mean(stdlfc[nidx]));
	  } else {
	  for(k in (j+1):ncol(comb)) {
	  itr=itr+1
	  z=y[,c(comb[,j],comb[,k])]
	  vcpm=voom(z,dmat,plot=FALSE)
	  fits=lmFit(vcpm,dmat)
	  efits=eBayes(fits)
	  pv=efits$p.value[,2,drop=FALSE]
	  lfc=efits$coefficients[,2,drop=FALSE]
	  
          colnames(pv)="PV_E.vs.C";colnames(lfc)="logFC_E.vs.C"
        
##logfc=log(E/C)       E/M is discussed
 
	stdlfc=sign(lfc)*abs(qnorm(pv/2+1E-16,0,1))
#        stdlfc[pv==1]=0
#        x11();
#        par(mfcol=c(2,2))
#        hist(stdlfc,30,main="ALL")
#        hist(stdlfc[midx],30,main="male")
#        hist(stdlfc[fidx],30,main="female")
#        hist(stdlfc[nidx],30,main="neutral")

        f=c(f,mean(stdlfc[fidx]));
        m=c(m,mean(stdlfc[midx])); 
        u=c(u,mean(stdlfc[nidx]));
print(c(i,j,k))
	}
	}
	}
######
#        n1=length(ff)
#        n2=length(mm)
#        n3=length(uu)
#        f=mean(ff);m=mean(mm);u=mean(uu)
#        for(j in 1:itr) {
#                s=sample(1:n1, n1 %/% 2)
#                tmp= (mean(ff[s]) - mean(ff[-s]))
#                f=c(f,tmp)
#
#                s=sample(1:n2, n2 %/% 2)
#                tmp= (mean(mm[s]) - mean(mm[-s]))
#                m=c(m,tmp)
#
#                s=sample(1:n3, n3 %/% 2)
#                tmp= (mean(uu[s]) - mean(uu[-s]))
#                u=c(u,tmp)
#        }
        pdf(paste0(outpath,"/",dataname[i],".pdf"), width=12, height=4)
			par(mfcol=c(1,3))
			par(mar=c(5,5,4,3))
                hist(m, 10, main="male-biased", xlab=""); abline(v=m[1], col=2, lwd=3)
                hist(f,10,main="female-biased",xlab=""); abline(v=f[1],col=2, lwd=3)
                hist(u,10,main="unbiased",xlab=""); abline(v=u[1],col=2, lwd=3)
                mtext(dataname[i],cex=1.4,outer=TRUE,line=-1.5)
        dev.off()

        png(paste0(outpath,"/",dataname[i],".png"))
                par(mfcol=c(2,2))
                hist(m,10,main="male biased gene set",xlab=""); abline(v=m[1],col=2)
                hist(f,10,main="female biased gene set",xlab=""); abline(v=f[1],col=2)
                hist(u,10,main="neutral gene set",xlab=""); abline(v=u[1],col=2)
                mtext(dataname[i],cex=1.4,outer=TRUE,line=-1.5)
        dev.off()
        print(paste0(dataname[i],": \n"))
        p[i,] = c(sum(m < m[1]),sum(f < f[1]),sum(u < u[1]),sum(m > m[1]),sum(f > f[1]),sum(u > u[1]))/itr
        male=cbind(male,m)
        female=cbind(female,f)
        neutral=cbind(neutral,u)
}

colnames(male) <- colnames(female) <- colnames(neutral) <- dataname
rownames(male) <- rownames(female) <- rownames(neutral) <- c("measured",paste0("iteration_",1:itr))
p = cbind(data=rownames(p),p)
Male = cbind(label=rownames(male),male)
Female = cbind(label=rownames(female),female)
Neutral = cbind(label=rownames(neutral),neutral)

write.table(p, paste0(outpath,"/p.values.txt"),sep="\t",row.names=F)
write.table(Male, paste0(outpath,"/male.txt"),sep="\t",row.names=F)
write.table(Female, paste0(outpath,"/female.txt"),sep="\t",row.names=F)
write.table(Neutral, paste0(outpath,"/neutral.txt"),sep="\t",row.names=F)

###  plot

png(paste0(outpath,"/hist_femaleSample.png"),width=960,height=720)
par(mfcol=c(3,4))
for(i in 1:4) {
if(i==1) {
hist(male[,i], 10, xlab="", ylab="Freq_maleSet", main=colnames(male)[i]); abline(v=male[1,i], col=2, lwd=2)
hist(female[,i], 10, xlab="", ylab="Freq_femaleSet", main=""); abline(v=female[1,i], col=2)
hist(neutral[,i], 10, xlab="Mean stdLogFC", ylab="Freq_neutralSet", main=""); abline(v=neutral[1,i], col=2)
} else {
hist(male[,i],10,xlab="",ylab="",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],10,xlab="",ylab="",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],10,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,i],col=2)
}
}
dev.off()

pdf(paste0(outpath,"/hist_femaleSample.pdf"),width=7.2,height=5)
par(mfcol=c(3,4))
for(i in 1:4) {
if(i==1) {
hist(male[,i],10,xlab="",ylab="Freq_maleSet",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],10,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],10,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,i],col=2)
} else {
hist(male[,i],10,xlab="",ylab="",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],10,xlab="",ylab="",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],10,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,i],col=2)
}
}
dev.off()

png(paste0(outpath,"/hist_maleSample.png"),width=960,height=720)
par(mfcol=c(3,4))
for(i in 1:4) {j=4+i
if(i==1) {
hist(male[,j],10,xlab="",ylab="Freq_maleSet",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],10,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],10,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,j],col=2)
} else {
hist(male[,j],10,xlab="",ylab="",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],10,xlab="",ylab="",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],10,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,j],col=2)
}
}
dev.off()

pdf(paste0(outpath,"/hist_maleSample.pdf"),width=960,height=720)
par(mfcol=c(3,4))
for(i in 1:4) {j=4+i
if(i==1) {
hist(male[,j],10,xlab="",ylab="Freq_maleSet",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],10,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],10,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,j],col=2)
} else {
hist(male[,j],10,xlab="",ylab="",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],10,xlab="",ylab="",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],10,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,j],col=2)
}
}
dev.off()

