#draw TMRCA plot 
rm(list=ls())
library(stats)
vcf_file="eg_vcf.vcf"
vcf_data=read.table(vcf_file)
hetro_base=vcf_data$V2[vcf_data$V7=="PASS"]
TMRCA_file="TMRCA" # this should be flexible
mydata=read.table(TMRCA_file)
n_row=dim(mydata)[1]
n_col=dim(mydata)[2]-1



base=mydata[,1]
base_min=min(base)
base_max=max(base)

TMRCA=as.matrix(mydata[,-1])
TMRCA_mean=rowMeans(mydata[,-1])
TMRCA_min=min(TMRCA_mean)
TMRCA_max=max(TMRCA_mean)
pdf("TMRCA.pdf")
#plot(c(base_min,base_max),c(TMRCA_min,TMRCA_max),type='n')
plot(base, TMRCA_mean, type="S",ylab="TMRCA")
#par(cin=10)
for (i in 1:(n_row-1)){
	ME = 1.96 * sqrt(var(TMRCA[i,])/n_col)
	lines(c(base[i],base[i]),c(TMRCA_mean[i]+ME,TMRCA_mean[i]-ME),type="p",pch=20,col="blue")
	#lines(c(base[i],base[i+1],base[i+1]),c(TMRCA_mean[i],TMRCA_mean[i],TMRCA_mean[i+1]), type="l")

}

lines(hetro_base,rep(TMRCA_min,length(hetro_base)),type="p",pch=20,col="red")
dev.off()

pdf("TMRCA-heat.pdf")

x=base
y=TMRCA[,1]
for (i in 2:(n_col)){
x=c(x,base)
y=c(y,TMRCA[,i])
}

df=data.frame(x,y)
k <- with(df,MASS:::kde2d(x,y))
filled.contour(k)

#for (i in 1:(n_row-1)){
#	ME = 1.96 * sqrt(var(TMRCA[i,])/n_col)
#	lines(c(base[i],base[i]),c(TMRCA_mean[i]+ME,TMRCA_mean[i]-ME),type="p",pch=20,col="blue")
#	lines(c(base[i],base[i+1],base[i+1]),c(TMRCA_mean[i],TMRCA_mean[i],TMRCA_mean[i+1]), type="l")

#}
dev.off()
