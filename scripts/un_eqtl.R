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

plot(hots1, cex.lab = 1.5, cex.axis = 1.5, chr = "A03")

library(eqtl)
?define.peak
brassica_peaks <- define.peak(scanone_imp_tot, th = 4.0, si = 1.0,
                               lodcolumn= "all", chr=c("A03"))
attributes(brassica_peaks)$scanone


tail(brassica_peaks, 20)
dim(brassica_peaks)


# get significant peaks
brassica_peaks[2][1]
is.na(brassica_peaks[2][1])


test <- brassica_peaks[1:10]
class(test)
head(test)
?pull.pheno
testpheno <- colnames(pull.pheno(brass_total, pheno.col = 1:10))
testpheno

test2 <- as.data.frame(is.na(sapply(test,`[`,1)))
str(test2)
test2

test
test3 <- test[!is.na(sapply(test,`[`,1))]
str(test3)

# QUICK AND DIRTY, build back up peak object
attributes(test3)$class <- c("peak", "list")
attributes(test3)$features <- c("lod", "mname.peak", "peak.cM", 
        "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
attributes(test3)$scanone <- "scanone_imp_tot"
attributes(test3)$lod.th <- 4.0
attributes(test3)$si <- 1.0
attributes(test3)$window <- 20
attributes(test3)
# got it

peaks_red <- brassica_peaks[!is.na(sapply(brassica_peaks,`[`,1))]
attributes(peaks_red)$class <- c("peak", "list")
attributes(peaks_red)$features <- c("lod", "mname.peak", "peak.cM", 
        "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
attributes(peaks_red)$scanone <- "scanone_imp_tot"
attributes(peaks_red)$lod.th <- 4.0
attributes(peaks_red)$si <- 1.0
attributes(peaks_red)$window <- 20
attributes(peaks_red)
length(peaks_red)

