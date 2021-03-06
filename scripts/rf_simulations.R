setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/scripts")
source("RF_functions.R")

#generic one chromosome simulations
set.seed(2341)
geno <- simgeno(124)
head(geno)
dim(geno)
set.seed(213756)
pheno <- simtrait(geno[, 200], 1, 0.4) + simtrait(geno[, 500], 0.75) - 
          simtrait(geno[,750], 1.5, 0.3)
pheno
dim(pheno)
length(pheno)
pdens(pheno, main = "distribution of quantitative trait")


library(qtl)
#some of this involves random draws
set.seed(123567)
# simulate brassica cross
# this is for version 1.2 NOT 1.5 of the map
n.ind <- 124
genlen <- c(87, 110, 137, 61, 85, 91, 102, 75, 144, 102)
brass_markers <- c(136, 100, 212, 89, 72, 134, 139, 71, 177, 143)

mymap <- sim.map(len = genlen, n.mar = brass_markers,
	eq.spacing = FALSE, include.x = FALSE, anchor.tel=TRUE)
plot(mymap)
summary(mymap)

?sim.cross
br_test_cross <- sim.cross(map = mymap, n.ind = 124, type = "riself")
summary(br_test_cross)

head(br_test_cross$geno)
genotypes <- pull.geno(br_test_cross)
head(genotypes)
str(genotypes)
genotypes[genotypes == 2] <- 0
head(genotypes)
colnames(genotypes)
geno.names <- dimnames(genotypes)[[2]]


# Demo values using brassica marker simulated data
# sample significant markers for our 4 simulated phenotypes
# or define them as you wish
m1 <- sample(geno.names, 3, replace = FALSE)
m1
m2 <- sample(geno.names, 2, replace = FALSE)
m3 <- sample(geno.names, 2, replace = FALSE)
m4 <- sample(geno.names, 1, replace = FALSE)

## get marker genotypes
g11 <- (genotypes[,m1[1]])
g11
g12 <- (genotypes[,m1[2]]) 
g13 <- (genotypes[,m1[3]])
g21 <- (genotypes[,m2[1]])
g22 <- (genotypes[,m2[2]])
g31 <- (genotypes[,m3[1]])
g32 <- (genotypes[,m3[2]])
g41 <- (genotypes[,m4[1]]) 

# function (orginal by Micheal Kuhn) to convert 0, 1 genotype assignments at markers into QTL with 
# tunable effect size and variance
br_traits <- function(geno, effect = 1, variance = 0.1){
 variance <- variance*abs(effect)
 geno <- as.logical(geno)
 pheno <- numeric(length(geno))
 pheno[!geno] <- rnorm(sum(!geno), 0, variance)
 pheno[geno] <- rnorm(sum(geno), effect, variance)
 return(pheno)
}

y1 <- br_traits(g11, effect = 2, variance = 0.1) + br_traits(g12, effect = 2, variance = 0.1)
plot(y1)

# replace trait columns what trait names
br_test_cross$pheno$trait1 <- br_traits(g11, effect = 0.3, variance = 0.5) + br_traits(g12, effect = 0.3, variance = 0.5)
br_test_cross$pheno

#epistatic model

m3[1] #g31
# [1] "D9M120"

m2[1] #g22
# [1] "D3M58"

m1[3] #for g13
# [1] "D2M33"

g31epi <- g31*0.20
g22epi <- g22*0.20
geno_epi <- g31epi + g22epi
geno_epi
br_test_cross$pheno$trait2 <- br_traits(g13, effect = 0.3, variance = 0.5) + geno_epi
br_test_cross$pheno





scanout <- scanone(br_test_cross, pheno.col = 2)
plot(scanout)

scanout2 <- scanone(br_test_cross, pheno.col = 3)
plot(scanout2)



scantwoout <- scantwo(br_test_cross, pheno.col = 3)
plot(scantwoout)

br_geno <- pull.geno(br_test_cross)
head(br_geno)
tail(br_geno)
br_geno[br_geno == 2] <- 0
dim(br_geno)
br_pheno <- as.matrix(pull.pheno(br_test_cross)[3])
br_pheno

br_geno_epi <- pull.geno(br_test_cross)
head(br_geno_epi)
tail(br_geno_epi)
br_geno_epi[br_geno_epi == 2] <- 0
dim(br_geno_epi)
br_pheno_epi <- as.matrix(pull.pheno(br_test_cross)[3])
br_pheno_epi


library(randomForest)
?randomForest



rf <- randomForest(y = br_pheno, x = br_geno, ntree = 4000)
sf <-  rfsf(rf)
plot(sf)

corr <- estBias(br_geno, 10000, verbose = F)
plot(corr, type = "h", ylab = "selection freq bias", main = "bias correction factor")

sf_corr <- sf - corr
sf_corr[sf_corr < 0 ] <- 0


par(mfrow = c(2, 1))
plot(sf, type = "h", ylab = "select freq", main = "RFSF, uncorrected")
points(c(200, 500, 750), sf[c(200, 500, 750)], col = "red", lwd = 1.5)



plot(sf_corr, type = "h", ylab = "adjusted select freq", main = "RFSF, corrected")
points(c(200, 500, 750), sf_corr[c(200, 500, 750)], col = "red", lwd = 1.5)
sf_corr
expr <- matrix(rnorm(8 * 100), 8, 100)
expr




################
rf_epi <- randomForest(y = br_pheno_epi, x = br_geno_epi, ntree = 4000)
sf_epi <-  rfsf(rf_epi)
plot(sf_epi)

corr_epi <- estBias(br_geno_epi, 10000, verbose = F)
plot(corr_epi, type = "h", ylab = "selection freq bias", main = "bias correction factor")

sf_corr_epi <- sf_epi - corr_epi
sf_corr_epi[sf_corr_epi < 0 ] <- 0


par(mfrow = c(2, 1))
plot(sf_epi, type = "h", ylab = "select freq", main = "RFSF, uncorrected")
points(c(200, 500, 750), sf[c(200, 500, 750)], col = "red", lwd = 1.5)



plot(sf_corr_epi, type = "h", ylab = "adjusted select freq", main = "RFSF, corrected")
points(c(200, 500, 750), sf_corr[c(200, 500, 750)], col = "red", lwd = 1.5)
sf_corr
expr <- matrix(rnorm(8 * 100), 8, 100)
expr







# compare to EBglmnet
setwd("/Users/Cody_2/git.repos/EBglmnet/data/")
library(EBglmnet)
data(BASIS)
head(BASIS)
data(y)
y
n <- 50
k <- 100
BASIS <- BASIS[1:n, 1:k]
y <- y[1:n]

CV <- EBlassoNEG.GaussianCV(BASIS, y, nFolds = 3, Epis = "no")
CV
output <- EBlassoNEG.Gaussian(BASIS, y, a_gamma = 0.1, b_gamma = 0.1, Epis = "yes")
output


CV2 <- EBlassoNEG.GaussianCV(geno, pheno, nFolds = 3, Epis = "no")
output <- EBlassoNEG.Gaussian(geno, pheno, a_gamma = 0.1, b_gamma = 0.1, Epis = "yes")
output
# gets close, but picks 1 or two bins over

CV3 <- EBlassoNEG.GaussianCV(br_geno, br_pheno, nFolds = 3, Epis = "no")
output3 <- EBlassoNEG.Gaussian(br_geno, br_pheno, a_gamma = 0.1, b_gamma = 0.1, Epis = "yes")
output3
plot

CV4 <- EBelasticNet.GaussianCV(br_geno, br_pheno, nFolds = 3, Epis = "no")
CV4
?EBelasticNet.Gaussian
output4 <- EBelasticNet.Gaussian(br_geno, br_pheno, lambda = 1, alpha = 0.1, Epis = "yes")
output4
# gets close, but picks 1 or two bins over
