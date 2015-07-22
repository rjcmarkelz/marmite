
###########
# cody markelz 
# quick hacking to make some progress
# April 30, 2015
###########


#######
# 2010
#######
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
brass_2010 <- read.table("brassica_field_2010_blups.csv", sep = ",", header = TRUE)
head(brass_2010)
brass_2010$Line <- sub("(L_)(\\d+)", "RIL_\\2", brass_2010$Line)

#######
# 2012
#######
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
brass_2012 <- read.table("Brassica2012allmeans.csv", sep = ",", header = TRUE)
head(brass_2012)
brass_2012$id <- sub("(\\d+)", "RIL_\\1", brass_2012$id)

###############
# merge 2010 2012
###############
field_data <- merge(brass_2010, brass_2012, by.x = "Line", by.y = "id", all = TRUE)
dim(field_data)

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
write.table(field_data, "2010_2012_field_data.csv", sep = ",")

#make sure directory is set
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")

library(qtl)
library(qtlbim)
field_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="all_traits_RQTL.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(field_traits)
plot(field_traits)

brassica_traits_qb <- field_traits
brassica_traits_qb <- qb.genoprob(brassica_traits_qb, step=2, stepwidth = "variable")
summary(brassica_traits_qb)

brassica_traits_qb_leaflength <- qb.mcmc(brassica_traits_qb, pheno.col = 2, seed = 1616, epistasis = TRUE)
#3000 iterations
#3000 iterations
plot(brassica_traits_qb_leaflength)
?qb.hpdone
other_HDI <- qb.hpdone(brassica_traits_qb_leaflength, effects = "estimate")
str(other_HDI)
other_HDI
hist(brassica_traits_qb$pheno$abbiomass)
ph_names <- names(brassica_traits_qb$pheno)
ph_names

# fitness correlated across years
cor(field_traits$pheno[57], field_traits$pheno[28], use = "pairwise.complete.obs")

cor(field_traits$pheno[52], field_traits$pheno[28], use = "pairwise.complete.obs")
cor(field_traits$pheno[52], field_traits$pheno[57], use = "pairwise.complete.obs")

cor(field_traits$pheno[77], field_traits$pheno[57], use = "pairwise.complete.obs")

phenos <- field_traits$pheno
head(phenos)

library(corrgram)

corrgram(phenos, order=TRUE, lower.panel=panel.shade,
  upper.panel=panel.pie, text.panel=panel.txt,
  main="")


dim(phenos)

correlations <- as.data.frame(vapply(
  phenos[, -c(57,134)],
  function(x)
  {
    cor(phenos[, 57], x, use = "pairwise.complete.obs")
  },
  numeric(1)
))
head(correlations)
correlations[which.max(abs(correlations))]
str(correlations)
correlations


str(field_traits)
field_traits
set.seed(123)
pt <- scanone(field_traits, method = "hk", n.perm = 1000)
alphas <- seq(0.01, 0.10, by=0.01)
lod.thrs <- summary(pt, alphas)
lod.thrs
lod.thr <- lod.thrs[5]
scan1 <- scanone(field_traits, pheno.col = 1:133, method = "hk")
head(scan1)

library(qtlhot)
high1 <- highlod(scan1, lod.thr = min(lod.thrs), drop.lod = 1.5)
max(high1, lod.thr = lod.thrs)
hots1 <- hotsize(high1, lod.thr = lod.thr)
summary(hots1)
head(hots1)
hots1
str(hots1)
plot(hots1, cex.lab = 1.5, cex.axis = 1.5)
str(scan1)
head(scan1)

field_traits
?cim

cimout <- cim(field_traits, pheno.col = 1:133)
plot(cimout)

scanout <- t(as.data.frame(scan1["A10x12178600",]))
scanout

UN_a <- scanone(field_traits, pheno.col = "CRgerm_to_Flr", method = "hk")
plot(UN_a)

# found gene of interest on chromosome 10
# associated with many traits at that genomic position

# find genes in interval
# 









