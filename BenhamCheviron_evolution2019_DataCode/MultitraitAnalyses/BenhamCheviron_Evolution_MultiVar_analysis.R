library(gdata)
library(ade4)
library(lavaan)
library(ncf)

#This script contributes to multivariat analyses presented in figure 7 of Benham & Cheviron Evolution 2019.
#Code and functions derived from Stuart et al. Nat Ecol & Evol r code and modified for the Savannah Sparrow data.

# === f.theta  ===
f.theta <- function(table.of.diff.in.means, unique.id, select.col){
  angle.matrix.radians <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      angle.matrix.radians[i,j] <- round(acos(cor(x = t(table.of.diff.in.means[i,select.col]), t(table.of.diff.in.means[j,select.col]), use = "na.or.complete")), 3) #
    }
  }
  print(angle.matrix.radians)
  rownames(angle.matrix.radians) <- unique.id
  colnames(angle.matrix.radians) <- unique.id
  angle.matrix.degrees <- round(angle.matrix.radians*(180/pi), 3)
  angle.output <- list(angle.matrix.radians, angle.matrix.degrees)
  names(angle.output) <- c("theta.radians", "theta.degrees")
  return(angle.output)
}
# === end) f.theta  ===

# === f.deltaL ====
f.deltaL <- function(table.of.diff.in.means, select.col, unique.id){
  get.vectlength <- function(vect){
    return(sqrt(sum(vect^2, na.rm = TRUE)))
  }
  table.of.diff.in.means.t <- table.of.diff.in.means[,select.col]
  length.of.vector.by.watershed <- apply(table.of.diff.in.means.t, 1 ,get.vectlength) #1 signifies by row. apply the function to every row in the input dataframe
  length.diff.matrix <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      length.diff.matrix[i,j] <-round(length.of.vector.by.watershed[i] - length.of.vector.by.watershed[j], 3)
    }
  }
  rownames(length.diff.matrix) <- unique.id
  colnames(length.diff.matrix) <- unique.id
  
  length.of.vector.by.watershed <- as.data.frame(length.of.vector.by.watershed)
  length.of.vector.by.watershed$wshd <- unique.id
  length.output <- list(length.diff.matrix, length.of.vector.by.watershed) 
  return(length.output)
}
# === end) f.deltaL ====
#calculate theta and deltaL for all datasets
SampleTable<-read.csv("~/Data/MultitraitAnalyses/Nonparallelism_physiologyData.csv", header=TRUE)
#calculate t.statistics, theta, deltaL, L for all traits
Locations<-levels(SampleTable$locality)
Cols<-colnames(SampleTable[,c(8:13)])
traits = length(Cols)
pops = length(Locations)-1 


trait.output <- matrix(ncol = traits, nrow = pops)

for (k in 1:(length(Locations)-1)) {
	print(Locations[k])
	TidMarsh1<-subset(SampleTable, locality == Locations[k] | locality == "Upland")
	
	for (m in 8:13){
		tryCatch({
		trait.t<-t.test(TidMarsh1[,m]~TidMarsh1$locality)
		trait.output[k,m-7]<-trait.t$statistic
		}, error = function(e){cat("ERROR :",conditionMessage(e),"\n")})
	}
}

trait.output<-data.frame(trait.output)
rownames(trait.output)<-Locations[1:9]
colnames(trait.output)<-Cols
print(trait.output)
#calculate theta, deltaL, and L for traits
#traits.theta
unique.id <- Locations[1:9]
select.col <- c(1:6)
traits.theta <- f.theta(trait.output, unique.id, select.col)
traits.theta.radians <- traits.theta[[1]]
traits.theta.degrees <- traits.theta[[2]]
all.trait.theta <- as.dist(traits.theta.radians, diag = FALSE, upper = FALSE)


(mean.theta.morpho.deg <- mean(all.trait.theta) * (180/pi))
(sd.theta.morpho.deg <- sd(all.trait.theta)*(180/pi))
(se.theta.morpho.deg <- (sd(all.trait.theta)/sqrt(length(all.trait.theta)))*(180/pi))
(range.theta.morpho.deg <- (c(min(all.trait.theta), max(all.trait.theta))) * (180/pi)) 

#traits deltaL and L
traits.delta.L <- f.deltaL(trait.output, select.col, unique.id)
traits.deltaL<-traits.delta.L[[1]]
traits.L<-traits.delta.L[[2]]
traits.deltaL.dist<-as.dist(traits.deltaL, diag=FALSE, upper = FALSE)

(mean.dL.morpho.deg <- mean(traits.deltaL.dist))
(sd.dL.morpho.deg <- sd(traits.deltaL.dist))
(se.dL.morpho.deg <- sd(traits.deltaL.dist)/sqrt(length(traits.deltaL.dist)))
(range.dL.morpho.deg <- c(min(traits.deltaL.dist), max(traits.deltaL.dist)))



par(mfrow=c(1,2))
hist(upperTriangle(traits.theta.degrees))
hist(upperTriangle(traits.deltaL))
########################################################################
#z-transform, calculate mean differences between upland and each tidal marsh, calculate theta, etc.
env.scale<-scale(SampleTable[,c(14:36)], center=TRUE, scale=TRUE)
SampleTable<-cbind(SampleTable[,c(1:13)],env.scale)
SampleTable<-as.data.frame(SampleTable)
Locations.env<-levels(SampleTable$locality)
Cols.env<-colnames(SampleTable[,c(14:36)])
env.var = length(Cols.env)
pops.env = length(Locations.env)-1
#14:36

env.output <- matrix(ncol = env.var, nrow = pops.env)

for (k in 1:(length(Locations.env)-1)) {
	print(Locations.env[k])
	TidMarsh1<-subset(SampleTable, locality == Locations.env[k] | locality == "Upland")

	for (m in 14:36){
		env.mean<-mean(TidMarsh1[,m][TidMarsh1$locality=="Upland"], na.rm=TRUE) - mean(TidMarsh1[,m][TidMarsh1$locality==Locations.env[k]], na.rm=TRUE)
		env.output[k,m-13]<-env.mean
		
	}
}	

env.output<-data.frame(env.output)
rownames(env.output)<-Locations.env[1:9]
colnames(env.output)<-Cols.env

#calculate theta, deltaL, and L for env variables
#traits.env
unique.id.env <- Locations.env[1:9]
select.col.env <- c(1:23)
env.theta <- f.theta(env.output, unique.id.env, select.col.env)
env.theta.radians <- env.theta[[1]]
env.theta.degrees <- env.theta[[2]]
all.env.theta <- as.dist(env.theta.degrees, diag = FALSE, upper = FALSE)

#traits deltaL and L
env.delta.L <- f.deltaL(env.output, select.col.env, unique.id.env)
env.deltaL<-env.delta.L[[1]]
env.L<-env.delta.L[[2]]
env.deltaL.dist<-as.dist(env.deltaL, diag=FALSE, upper = FALSE)

########################################################################
#calculate mean differences between upland and each tidal marsh, calculate theta, etc. for
#genetic data
snp.pca<-read.csv("~/Data/MultitraitAnalyses/snp.pca.dat.csv",header=TRUE)

Locations.gen<-levels(snp.pca$locality)
Cols.gen<-colnames(snp.pca[,c(4:63)])
gen.var = length(Cols.gen)
pops.gen = length(Locations.gen)-1

gen.output <- matrix(ncol = gen.var, nrow = pops.gen)

for (k in 1:(length(Locations.gen)-1)) {
	print(Locations.gen[k])
	TidMarsh1<-subset(snp.pca, locality == Locations.gen[k] | locality == "Upland")

	for (m in 4:63){
		gen.mean<-mean(TidMarsh1[,m][TidMarsh1$locality=="Upland"], na.rm=TRUE) - mean(TidMarsh1[,m][TidMarsh1$locality==Locations.gen[k]], na.rm=TRUE)
		gen.output[k,m-3]<-gen.mean	
	}
}	

gen.output<-data.frame(gen.output)
rownames(gen.output)<-Locations.gen[1:9]
colnames(gen.output)<-Cols.gen

#calculate theta, deltaL, and L for gen variables
#traits.gen
unique.id.gen <- Locations.gen[1:9]
select.col.gen <- c(1:60)
gen.theta <- f.theta(gen.output, unique.id.gen, select.col.gen)
gen.theta.radians <- gen.theta[[1]]
gen.theta.degrees <- gen.theta[[2]]
all.gen.theta <- as.dist(gen.theta.degrees, diag = FALSE, upper = FALSE)

#traits deltaL and L
gen.delta.L <- f.deltaL(gen.output, select.col.gen, unique.id.gen)
gen.deltaL<-gen.delta.L[[1]]
gen.L<-gen.delta.L[[2]]
gen.deltaL.dist<-as.dist(gen.deltaL, diag=FALSE, upper = FALSE)


########################################################################
Fst<-read.csv("~/Data/MultitraitAnalyses/FstMat.csv", header=TRUE)
Fst<-as.matrix(Fst[,-1])
Fst<-Fst/(1-Fst)
Fst.dist<-as.dist(Fst, diag=FALSE, upper = FALSE)
print(Fst.dist)

Fst_L<-c(0.113,0.12,0.04,0.114,0.059,0.078,0.067,0.101,0.061)

#par(mfrow=c(3,2))
#par(mar=c(4,4.5,2,2))
#plot(all.trait.theta~all.env.theta, pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
#abline(lm(all.trait.theta~all.env.theta))
#all.gene.theta
#plot(all.trait.theta~all.gen.theta, pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
#abline(lm(all.trait.theta~all.gen.theta))

#plot(traits.L[,1]~env.L[,1], pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
env.lm<-lm(traits.L[,1]~env.L[,1])
print(summary(env.lm))
#abline(env.lm)

#plot(traits.L[,1]~gen.L[,1], pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
gen.lm<-lm(traits.L[,1]~gen.L[,1])
print(summary(gen.lm))
#abline(gen.lm)

#plot(traits.deltaL.dist~env.deltaL.dist, pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
#abline(lm(traits.deltaL.dist~env.deltaL.dist))

#plot(traits.deltaL.dist~gen.deltaL.dist, pch = 21, bg="black", cex.axis = 1.5, cex.lab = 1.5)
#abline(lm(traits.deltaL.dist~gen.deltaL.dist))

##mantel tests
env.mantel.theta<- mantel.rtest(all.env.theta, all.trait.theta, nrepet = 9999)
gen.mantel.theta<- mantel.rtest(all.gen.theta, all.trait.theta, nrepet = 9999)
print("env theta mantel test")
print(env.mantel.theta)

print("gen theta mantel test")
print(gen.mantel.theta)


env.mantel.deltaL<- mantel.rtest(env.deltaL.dist, traits.deltaL.dist, nrepet = 9999)
gen.mantel.deltaL<- mantel.rtest(gen.deltaL.dist,  traits.deltaL.dist, nrepet = 9999)

print("env deltaL mantel test")
print(env.mantel.deltaL)

print("gen deltaL mantel test")
print(gen.mantel.deltaL) 

trait.theta<-as.matrix(all.trait.theta)
env.theta<-as.matrix(all.env.theta)
gen.theta<-as.matrix(all.gen.theta)

pmt<-partial.mantel.test(trait.theta,env.theta,gen.theta,resamp=10000)

env.deltaL<-as.matrix(env.deltaL.dist)
gen.deltaL<-as.matrix(gen.deltaL.dist)
trait.deltaL<-as.matrix(traits.deltaL.dist)

pmt2<-partial.mantel.test(trait.deltaL,env.deltaL,gen.deltaL, resamp=10000)

print(pmt)
print(pmt2)

print("theta")
print(mean.theta.morpho.deg)
print(sd.theta.morpho.deg)
print(se.theta.morpho.deg)
print(range.theta.morpho.deg)

print("deltaL")
print(mean.dL.morpho.deg)
print(sd.dL.morpho.deg)
print(se.dL.morpho.deg)
print(range.dL.morpho.deg) 

