#this script will create the SNP principal components "snp.pca.dat.csv" for input into script "BenhamCheviron_MultiVar_analysis.R" 

vcf<-read.table("~/Data/MultitraitAnalyses/vcf_file_editing/SAVS_westcoast.filtered.012")
vcf<-as.data.frame(vcf)
vcf[vcf == -1] <- NA
for (i in 1:length(vcf[1,])){
	vcf[,i][is.na(vcf[,i])] = mean(vcf[,i], na.rm=TRUE)
}

snp.pca<-prcomp(vcf)
print(summary(snp.pca))
snps.prcomp<-snp.pca$x[,c(1:84)]
snp.pca.dat<-cbind(samples,snps.prcomp)

write.csv(snp.pca.dat, file="~/Data/MultitraitAnalyses/snp.pca.dat.csv")
par(mfrow=c(1,1))
snp.pca<-read.csv("~/Data/MultitraitAnalyses/snp.pca.dat.csv",header=TRUE)

#plot PCA of SNP data to assess patterns of population structure.
plot(snp.pca[,4], snp.pca[,5], type='n')
loc<-levels(snp.pca[,3])
colors<-c("red","blue","green","orange","purple","dark gray","dark red","yellow", "dark green", "black")
for (i in 1:length(loc)){
	loc.sub<-subset(snp.pca, V2==loc[i])
	points(loc.sub[,4],loc.sub[,5], pch=21, bg=colors[i])
}
print(length(snp.pca[1,]))
