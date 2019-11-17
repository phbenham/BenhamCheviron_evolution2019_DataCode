library('compute.es')
library('pwr')

#script for calculating trait effect sizes, anova analyses and generating Figure 4 from Benham & Cheviron 2019
#Evolution.

phys.data2<-read.csv("~/Dropbox/Chap2_finalmaterials_forsubmission/CodeData/EffectSizeAnalyses_Fig4/BenhamCheviron_Evolution_CaliforniaPhysData.csv", header=TRUE)
attach(phys.data2) 

#subset all traits by habitat
TEWL.cal<-subset(phys.data2[,13], habitat=="tidal_cal")
TEWL.mex<-subset(phys.data2[,13], habitat=="tidal_mex")
TEWL.int<-subset(phys.data2[,13], habitat=="interior")

kidney.cal<-subset(phys.data2[,15], habitat=="tidal_cal")
kidney.mex<-subset(phys.data2[,15], habitat=="tidal_mex")
kidney.int<-subset(phys.data2[,15], habitat=="interior")

plasma.cal<-subset(phys.data2[,16], habitat=="tidal_cal")
plasma.mex<-subset(phys.data2[,16], habitat=="tidal_mex")
plasma.int<-subset(phys.data2[,16], habitat=="interior")

urine.cal<-subset(phys.data2[,17], habitat=="tidal_cal")
urine.mex<-subset(phys.data2[,17], habitat=="tidal_mex")
urine.int<-subset(phys.data2[,17], habitat=="interior")

UPratio.cal<-subset(phys.data2[,18], habitat=="tidal_cal")
UPratio.mex<-subset(phys.data2[,18], habitat=="tidal_mex")
UPratio.int<-subset(phys.data2[,18], habitat=="interior")

Medulla.cal<-subset(phys.data2[,19], habitat=="tidal_cal")
Medulla.mex<-subset(phys.data2[,19], habitat=="tidal_mex")
Medulla.int<-subset(phys.data2[,19], habitat=="interior")


TEWLcal.t<-t.test(TEWL.cal, TEWL.int)
TEWLmex.t<-t.test(TEWL.mex, TEWL.int)

kidneycal.t<-t.test(kidney.cal, kidney.int)
kidneymex.t<-t.test(kidney.mex, kidney.int)

plasmacal.t<-t.test(plasma.cal, plasma.int)
plasmamex.t<-t.test(plasma.mex, plasma.int)

urinecal.t<-t.test(urine.cal, urine.int)
urinemex.t<-t.test(urine.mex, urine.int)

UPratiocal.t<-t.test(UPratio.cal, UPratio.int)
UPratiomex.t<-t.test(UPratio.mex, UPratio.int)

medullacal.t<-t.test(Medulla.cal, Medulla.int)
medullamex.t<-t.test(Medulla.mex, Medulla.int)

tvalue<-c(kidneycal.t$statistic, kidneymex.t$statistic, medullacal.t$statistic, medullamex.t$statistic, plasmacal.t$statistic, plasmamex.t$statistic, urinecal.t$statistic, urinemex.t$statistic, UPratiocal.t$statistic, UPratiomex.t$statistic, TEWLcal.t$statistic, TEWLmex.t$statistic)
n.t<-c(14,38,14,22,24,25,20,25,13,25,26,29)
n.c<-c(17,17,14,14,16,16,15,15,14,14,16,16)

ttest.result<-data.frame(id=1:12, t=tvalue, n.t=n.t, n.c=n.c)
print(ttest.result)

d<-tes(t,n.1=n.t,n.2=n.c,level=95, dig=2, id=id, data=ttest.result)
print(d)

tewl.pwr<-pwr.t.test(d=-0.43, sig.level=0.05, power=0.8, type="two.sample")
kidmass.pwr<-pwr.t.test(d=3.12, sig.level=0.05, power=0.8, type="two.sample")
plasma.pwr<-pwr.t.test(d=2.50, sig.level=0.05, power=0.8, type="two.sample")
urine.pwr<-pwr.t.test(d=0.97, sig.level=0.05, power=0.8, type="two.sample")
UP.pwr<-pwr.t.test(d=0.68, sig.level=0.05, power=0.8, type="two.sample")
medulla.pwr<-pwr.t.test(d=2.67, sig.level=0.05, power=0.8, type="two.sample")

trait<-c("Kidney", "Medulla", "Plasma", "Urine", "U:P", "TEWL" )
sample.size<-c(kidmass.pwr$n, medulla.pwr$n, plasma.pwr$n, urine.pwr$n, UP.pwr$n, tewl.pwr$n)

trait.samplesizes<-cbind(trait, sample.size)
print(trait.samplesizes)


cal.d<-d[c(1,3,5,7,9,11),]
mex.d<-d[c(2,4,6,8,10,12),]

var.cal<-cal.d[,6]
sd.cal<-sqrt(var.cal)
CI.cal<-sd.cal*1.96

var.mex<-mex.d[,6]
sd.mex<-sqrt(var.mex)
CI.mex<-sd.mex*1.96

counta<-rbind(cal.d,mex.d)
count<-counta[,5]

location<-c("cal","cal","cal","cal","cal","cal","mex","mex","mex","mex","mex","mex")

trait<-c("A","B","C","D","E","F","A","B","C","D","E","F")


group.bar<-data.frame( location=location, trait=trait, count=count)
data=count
data=matrix(data,ncol=6,byrow=TRUE)
colnames(data)=levels(group.bar$trait)
rownames(data)=levels(group.bar$location)

ci.plot<-matrix(c(CI.cal,CI.mex),2,6,byrow=TRUE)
print(ci.plot)

error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) !=length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, lwd=0.5, length=length, ...)
 }


barz<-barplot(data, names.arg=c("Kidney mass","Medulla","Plasma","Urine","U:P ratio","TEWL"), ylim=c(-2,4), space=c(0,2), col=c("white","black"), beside=TRUE, ylab="Effect size (d)", lwd=2)
error.bar(barz,data,ci.plot)
abline(h=0, col="black", lwd=2)
legend(20,4, legend=c("California", "Mexico"), pch=c(21,21), pt.bg=c("white", "black"), cex=1.25)
detach(phys.data2)

#anova for each trait
medulla.aov<-aov(phys.data2$medulla~phys.data2$habitat)
print(TukeyHSD(medulla.aov))

kidney.aov<-aov(phys.data2$kidmass_corr~phys.data2$habitat)
print(TukeyHSD(kidney.aov))

plasma.aov<-aov(phys.data2$blood_osmo~phys.data2$habitat)
print(TukeyHSD(plasma.aov))

urine.aov<-aov(phys.data2$urine_osmo~phys.data2$habitat)
print(TukeyHSD(urine.aov))

UP.aov<-aov(phys.data2$UPratio~phys.data2$habitat)
print(TukeyHSD(UP.aov))

TEWL.aov<-aov(phys.data2$TEWL~phys.data2$habitat)
print(TEWL.aov)
print(TukeyHSD(TEWL.aov))


#power tests:
CalD<-cal.d$d
MexD<-mex.d$d

traits=c("Kidney mass","Medulla","Plasma","Urine","U:P ratio","TEWL")
for (i in 1:6) {
	cal.trait.pwr<-pwr.t.test(d=CalD[i], power=0.8, sig.level=0.05)	
	mex.trait.pwr<-pwr.t.test(d=MexD[i], power=0.8, sig.level=0.05)
	print(traits[i])
	print(cal.trait.pwr)
	print(mex.trait.pwr)
}
