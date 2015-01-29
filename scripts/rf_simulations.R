setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/scripts")
source("RF_functions.R")

#simulations

set.seed(2341)
geno <- simgeno(124)
head(geno)
dim(geno)
set.seed(213756)
pheno <- simtrait(geno[, 200], 1, 0.4) + simtrait(geno[, 500], 0.75) - 
          simtrait(geno[,750], 1.5, 0.3)
dim(pheno)
length(pheno)
pdens(pheno, main = "distribution of quantitative trait")

library(randomForest)
?randomForest
rf = randomForest(y = pheno, x = geno, ntree = 2000)
sf = rfsf(rf)
str(rf)
vu <- randomForest::varUsed(rf)
vu
sf2 <- vu/sum(vu)
rf$importance
plot(rf$importance)
sf2
plot(sf2)
corr = estBias(geno, 10000, verbose = F)
plot(corr, type = "h", ylab = "selection freq bias", main = "bias correction factor")

sf_corr <- sf - corr
sf_corr[sf_corr < 0 ] <- 0


par(mfrow = c(2, 1))
plot(sf, type = "h", ylab = "select freq", main = "RFSF, uncorrected")
points(c(200, 500, 750), sf[c(200, 500, 750)], col = "red", lwd = 1.5)

plot(sf_corr, type = "h", ylab = "adjusted select freq", main = "RFSF, corrected")
points(c(200, 500, 750), sf[c(200, 500, 750)], col = "red", lwd = 1.5)

expr <- matrix(rnorm(8 * 100), 8, 100)
expr

# compare to EBglmnet
setwd("/Users/Cody_2/git.repos/EBglmnet/data/")
library(EBglmnet)
data(BASIS)
head(BASIS)
data(y)
y
n = 50
k = 100
BASIS = BASIS[1:n, 1:k]
y = y[1:n]
CV = EBlassoNEG.GaussianCV(BASIS, y, nFolds = 3, Epis = "no")
CV
output = EBlassoNEG.Gaussian(BASIS, y, a_gamma = 0.1, b_gamma = 0.1, Epis = "yes")
output

CV2 <- EBlassoNEG.GaussianCV(geno, pheno, nFolds = 3, Epis = "no")
output <- EBlassoNEG.Gaussian(geno, pheno, a_gamma = 0.1, b_gamma = 0.1, Epis = "yes")
output
# gets close, but picks 1 or two bins over

