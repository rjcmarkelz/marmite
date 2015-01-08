###############
###############
# total eqtl no treatment
###############
###############
# br_blues_total_RQTL.csv was created as part of the data clean up in cr_un_eqtl.R
library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

brass_total <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="br_blues_total_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brass_total)

class(brass_total)[1] <- "riself"
brass_total <- jittermap(brass_total)
brass_total

brass_total <- est.rf(brass_total)
plot.rf(brass_total) 

#about a minute
brass_total <- calc.errorlod(brass_total, error.prob=0.001)

system.time(scanone_imp_tot <- scanone(brass_total, pheno.col = 1:35039, 
	         method = "imp", use="all.obs"))

save.image(file = "un_eqtl.RData", version = NULL,
 ascii = FALSE, safe = TRUE)

library(qtlhot)
set.seed(12345)

permtest <- scanone(brass_total, method = "imp", n.perm = 1000)
permtest

alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs <- summary(permtest, alphas)

# get 3% cutoff
lod.thrs
lod.thr <- lod.thrs[1]

#reduce object size by getting removing NS peaks
high1 <- highlod(scanone_imp_tot, lod.thr = lod.thr, drop.lod = 1.5)

max(high1, lod.thr = lod.thrs)
high1
hots1 <- hotsize(high1, lod.thr = lod.thr)
hots1

plot(hots1, cex.lab = 1.5, cex.axis = 1.5)





