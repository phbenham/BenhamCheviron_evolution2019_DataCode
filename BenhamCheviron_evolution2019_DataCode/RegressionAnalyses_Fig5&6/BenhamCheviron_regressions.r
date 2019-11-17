library(AICcmodavg)
library(ape)
library(ggplot2)
library(gridExtra)

phys.div<-read.csv("~/Dropbox/Dissertation/SAVS_physiology/TraitDivergence_tstat.csv", header=TRUE)

traits.z<-scale(phys.div[,c(3:14, 35,36,37)], center=TRUE, scale=TRUE)
#traits.z<-cbind(phys.div[,c(3:8)], traits.z)
traits.z<-as.data.frame(traits.z)

ylabels <- c("Kidney mass", "Medulla volume", "Plasma osmolality", "Urine osmolality", "Urine:Plasma ratio", "TEWL")
par(mfrow=c(3,2))
par(mar=c(4,4.5,2,2))

for (i in 1:6) {
	print(ylabels[i])
	lm1<-lm(traits.z[,i]~Migration, data=traits.z)
	print(summary(lm1))
	lm2<-lm(traits.z[,i] ~ poly(Migration, 2, raw=TRUE), data=traits.z)
	print(summary(lm2))
	
	trait.sum<-summary(lm1)
	
	r2<-trait.sum$adj.r.squared
	p<-trait.sum$coefficients[2,4]

	rp = vector('expression',2)
	rp[1] = substitute(expression(italic(R)^2 == value1),list(value1 = format(r2,dig=3)))[2]
	rp[2] = substitute(expression(italic(p) == value2),list(value2 = format(p,dig=3)))[2]
	
	mod.comp<-anova(lm1,lm2)
	print(mod.comp)
	
	newdat = data.frame(Migration <- seq(from=min(traits.z$Migration), to = max(traits.z$Migration), length.out=250))
	newdat$pred = predict(lm2, newdata=newdat)

	plot(traits.z[,i]~Migration, data=traits.z, xlab = "Migration", ylab=ylabels[i], pch=21, 
	bg="black", cex.axis=1.2, cex.lab=1.5, lwd=2)
	abline(lm1)
	if(i != 3){
		abline(lm1)
	}
	else {
		with(newdat,lines(x=Migration,y=pred))
	}
	
	if(i < 4){
		legend('topright', legend=rp, bty='n', cex=1.5)
	}	
	else {
		legend('bottomright', legend=rp, bty='n', cex=1.5)	
	}
}
