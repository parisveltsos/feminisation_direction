graphics.off()

# change to location where fileS3 is unzipped
datap="~/git/feminisation_direction/input/pseudo_counts"
setwd(datap)
dataf=list.files()[grep("EM[A-Z]",list.files())]
dataf=paste0(datap,"/",dataf)

itr=10000
outpath="~/git/feminisation_direction/output/results_balanced_randomisation"
if(!file.exists(outpath)) dir.create(outpath)

extrem=list("","FBgn0255061","","","","FBgn0255061","","")
extrem="FBgn0255061"

dataname <- substr(basename(dataf),1,5)

set.seed(1234)
p=matrix(1,length(dataf),6)
rownames(p)=dataname
colnames(p)=c("p.male.lower", "p.female.lower", "p.neutral.lower", "p.male.upper", "p.female.upper", <D-≤>"p.neutral.upper")

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
	pv=x[,"PValue"]
	lfc=x[,"logFC"]
##logfc=log(M/E)       E/M is discussed
	stdlfc = -sign(lfc)*abs(qnorm(pv/2+1E-16,0,1))
	stdlfc[pv==1]=0
# 	x11();
	par(mfcol=c(2,2))	
	hist(stdlfc,30,main="ALL")
	hist(stdlfc[midx],30,main="male")
	hist(stdlfc[fidx],30,main="female")
	hist(stdlfc[nidx],30,main="neutral")

	ff=(stdlfc[fidx]); #ff=ff[!is.na(ff)]
	mm=(stdlfc[midx]); #mm=mm[!is.na(mm)]
	uu=(stdlfc[nidx]); #uu=uu[!is.na(uu)]

	n1=length(ff)
	n2=length(mm)
	n3=length(uu)
	f=mean(ff);m=mean(mm);u=mean(uu)
	for(j in 1:itr) {
		s=sample(1:n1, n1 %/% 2)
		tmp= (mean(ff[s]) - mean(ff[-s]))
		f=c(f,tmp)

		s=sample(1:n2, n2 %/% 2)
		tmp= (mean(mm[s]) - mean(mm[-s]))
		m=c(m,tmp)

		s=sample(1:n3, n3 %/% 2)
		tmp= (mean(uu[s]) - mean(uu[-s]))
		u=c(u,tmp)
	}
	pdf(paste0(outpath,"/",dataname[i],".pdf"), width=12, height=4)
		par(mfcol=c(1,3))
		par(mar=c(5,5,4,3))
		hist(m, 100, main="male-biased", xlab=""); abline(v=m[1], col=2, lwd=2.5)
		hist(f, 100, main="female-biased", xlab=""); abline(v=f[1], col=2, lwd=2.5)
		hist(u, 100, main="unbiased",xlab=""); abline(v=u[1], col=2, lwd=2.5)
		mtext(dataname[i], cex=1.4, outer=TRUE, line=-1.5)
	dev.off()

	png(paste0(outpath,"/",dataname[i],".png"))
		par(mfcol=c(2,2))
		hist(m,100,main="male biased gene set",xlab=""); abline(v=m[1],col=2)
		hist(f,100,main="female biased gene set",xlab=""); abline(v=f[1],col=2)
		hist(u,100,main="neutral gene set",xlab=""); abline(v=u[1],col=2)
		mtext(dataname[i],cex=1.4,outer=TRUE,line=-1.5)
	dev.off()
	print(paste0(dataname[i],": \n"))
	p[i,] = c(sum(m <= m[1])-1,sum(f <= f[1])-1,sum(u <= u[1])-1,sum(m > m[1]),sum(f > f[1]),sum(u > u[1]))/itr
	male=cbind(male,m)
	female=cbind(female,f)
	neutral=cbind(neutral,u)
}
colnames(male) <- colnames(female) <- colnames(neutral) <- dataname
rownames(male) <- rownames(female) <- rownames(neutral) <- c("measured",paste0("iteration_",1:itr))
p=cbind(data=rownames(p),p)
Male=cbind(label=rownames(male),male)
Female=cbind(label=rownames(female),female)
Neutral=cbind(label=rownames(neutral),neutral)

write.table(p, paste0(outpath,"/p.values.txt"),sep="\t",row.names=F)
write.table(Male, paste0(outpath,"/male.txt"),sep="\t",row.names=F)
write.table(Female, paste0(outpath,"/female.txt"),sep="\t",row.names=F)
write.table(Neutral, paste0(outpath,"/neutral.txt"),sep="\t",row.names=F)

###  plot

png(paste0(outpath,"/hist_femaleSample.png"),width=960,height=720)
par(mfcol=c(3,4))
for(i in 1:4) {
if(i==1) {
hist(male[,i],100,xlab="",ylab="Freq_maleSet",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],100,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],100,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,i],col=2)
} else {
hist(male[,i],100,xlab="",ylab="",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],100,xlab="",ylab="",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],100,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,i],col=2)
}
}
dev.off()

pdf(paste0(outpath,"/hist_femaleSample.pdf"),width=7.2,height=5)
par(mfcol=c(3,4))
for(i in 1:4) {
if(i==1) {
hist(male[,i],100,xlab="",ylab="Freq_maleSet",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],100,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],100,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,i],col=2)
} else {
hist(male[,i],100,xlab="",ylab="",main=colnames(male)[i]);abline(v=male[1,i],col=2)
hist(female[,i],100,xlab="",ylab="",main="");abline(v=female[1,i],col=2)
hist(neutral[,i],100,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,i],col=2)
}
}
dev.off()

png(paste0(outpath,"/hist_maleSample.png"),width=960,height=720)
par(mfcol=c(3,4))
for(i in 1:4) {j=4+i
if(i==1) { 
hist(male[,j],100,xlab="",ylab="Freq_maleSet",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],100,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],100,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,j],col=2)
} else {
hist(male[,j],100,xlab="",ylab="",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],100,xlab="",ylab="",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],100,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,j],col=2)
}
}
dev.off()

pdf(paste0(outpath,"/hist_maleSample.pdf"),width=7.2,height=5)
par(mfcol=c(3,4))
for(i in 1:4) {j=4+i
if(i==1) {
hist(male[,j],100,xlab="",ylab="Freq_maleSet",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],100,xlab="",ylab="Freq_femaleSet",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],100,xlab="Mean stdLogFC",ylab="Freq_neutralSet",main="");abline(v=neutral[1,j],col=2)
} else {
hist(male[,j],100,xlab="",ylab="",main=colnames(male)[j]);abline(v=male[1,j],col=2)
hist(female[,j],100,xlab="",ylab="",main="");abline(v=female[1,j],col=2)
hist(neutral[,j],100,xlab="Mean stdLogFC",ylab="",main="");abline(v=neutral[1,j],col=2)
}
}
dev.off()


## end

