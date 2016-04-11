setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/scripts")
source("RF_functions.R")

library(randomForest)
library(qtl)
library(ggplot2)
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
geno.names
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

# function (original by Micheal Kuhn) to convert 0, 1 genotype assignments at markers into QTL with 
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


########################
scanout2 <- scanone(br_test_cross, pheno.col = 3)
plot(scanout2)
########################
rf_epi <- randomForest(y = br_pheno_epi, x = br_geno_epi, ntree = 4000)
sf_epi <-  rfsf(rf_epi)
plot(sf_epi)

corr_epi <- estBias(br_geno_epi, 10000, verbose = F)
plot(corr_epi, type = "h", ylab = "selection freq bias", main = "bias correction factor")

sf_corr_epi <- sf_epi - corr_epi
sf_corr_epi[sf_corr_epi < 0 ] <- 0


# comparison plots
par(mfrow = c(2, 1))
plot(scanout2)
# plot(sf_epi, type = "h", ylab = "select freq", main = "RFSF, uncorrected")

plot(sf_corr_epi, type = "h", ylab = "adjusted select freq", main = "RFSF, corrected")

#############################
# QTL on same chromosome for same trait
#############################
head(geno.names)
geno.names # get index for what markers you want

m6 <- geno.names[c(640, 650)]
m6
# [1] "D6M31"  "D6M121"
## get marker genotypes
g61 <- (genotypes[,m6[1]])
g61
g62 <- (genotypes[,m6[2]]) 
g62

m1[3] #for g13
# [1] "D2M33"

g61epi <- g61*0.80
g62epi <- g62*0.50
geno_epi <- g61epi + g62epi
geno_epi
br_test_cross$pheno$trait3 <- br_traits(g13, effect = 0.3, variance = 0.5) + geno_epi
br_test_cross$pheno

########################
scanout3 <- scanone(br_test_cross, pheno.col = 4)
plot(scanout3)

br_geno_epi <- pull.geno(br_test_cross)
head(br_geno_epi)
tail(br_geno_epi)
br_geno_epi[br_geno_epi == 2] <- 0
dim(br_geno_epi)
br_pheno_epi <- as.matrix(pull.pheno(br_test_cross)[4])
br_pheno_epi

rf_epi <- randomForest(y = br_pheno_epi, x = br_geno_epi, ntree = 4000)
sf_epi <-  rfsf(rf_epi)
plot(sf_epi)

corr_epi <- estBias(br_geno_epi, 10000, verbose = F)

sf_corr_epi <- sf_epi - corr_epi
sf_corr_epi[sf_corr_epi < 0 ] <- 0


# comparison plots
par(mfrow = c(1, 1))
plot(scanout3)
plot(scanout3, chr = 6)
# plot(sf_epi, type = "h", ylab = "select freq", main = "RFSF, uncorrected")

head(sf_corr_epi)
tail(sf_corr_epi)
plot(sf_corr_epi, type = "h", ylab = "adjusted select freq", main = "RFSF, corrected")
plot.map(br_test_cross, chr = 6)

scanout4 <- scanone(br_test_cross, pheno.col = 4)
set.seed(123454)
scanout4 <- cim(br_test_cross, pheno.col = 4)
plot(scanout4)

scanout4

flr_cim_A06 <- as.data.frame(subset(scanout4, chr = 6))
flr_cim_A06
plot(flr_cim_A06$pos, flr_cim_A06$lod)
peak2 <- 10
flr_A06 <- ggplot(flr_cim_A06)
flr_A06 <- flr_A06 +  theme_bw() + scale_y_continuous(limits=c(0, 40)) + 
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 2.87, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak2 * -0.02), yend = (peak2 * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") 
flr_A06
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/output")
ggsave("flr_A06_qtl_simulation.pdf", flr_A06, height = 10, width = 10)

head(sf_corr_epi)
qplot(sf_corr_epi)

