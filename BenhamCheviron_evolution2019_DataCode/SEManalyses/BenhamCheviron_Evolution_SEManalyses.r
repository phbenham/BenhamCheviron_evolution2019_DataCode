library(AICcmodavg)
library(lavaan)
library(ape)
library(glue)

#to do correction for spatial autocorrelation requires downloading the r script:
#lavSpatialCorrect.R from: https://github.com/jebyrnes/spatial_correction_lavaan

source('~/spatial_correction_lavaan-master/lavSpatialCorrect.R', chdir = TRUE)

phys.div<-read.csv("~/Dropbox/Chap2_finalmaterials_forsubmission/CodeData/SEManalyses/BenhamCheviron_Evolution_TraitDivergence.csv",header=TRUE)

traits.z<-scale(phys.div[,c(3:36)], center=TRUE, scale=TRUE)
traits<-names(phys.div[,3:8])
for (i in 1:6){
	trait.name<-traits[i]
	print(trait.name)
	mod1<- 'eco.comp =~ mayjune_tmax + mayjune_tmin + mayjune_prec + min_H2O_osmo   
			{trait.name} ~ Migration + eco.comp
			Migration ~~ eco.comp'
	
	mod1<-glue(mod1)
		   
	fit <- sem(mod1, data=traits.z, missing = "ML")
	
	mod.sum<-summary(fit, standardized=TRUE, rsq=TRUE)			
	print(mod.sum)
	
	#Does not work with missing data. Script needs to be modified to exclude missing data for some values, e.g. kidney 		#mass
	#mod1.spatialcorrect<-lavSpatialCorrect(fit,phys.div$latitude,phys.div$longitude)
	#print(mod1.spatialcorrect)
}	