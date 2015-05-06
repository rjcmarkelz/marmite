
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


############
# create RQTL dataset
############
head(field_data)
field_data$Line
# remove parental data in first two rows
field_data <- field_data[-c(1:2),]

field_data_t <- as.data.frame(t(field_data))
field_data_t
head(field_data_t)
colnames(field_data_t)
dim(field_data_t)
field_data_t[62,] <- field_data_t[1,]
rownames(field_data_t)[62] <- "id"
rownames(field_data_t)
field_data_t <- field_data_t[-1,]
head(field_data_t)
tail(field_data_t)

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
write.table(field_data_t, "2010_2012_field_data_RQTL.csv", sep = ",", col.names = FALSE)


library(qtl)
library(qtlnet)
field_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="2010_2012_field_data_RQTL.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(field_traits)

class(field_traits)[1] <- "riself"
field_traits <- jittermap(field_traits)
field_traits


field_traits <- sim.geno(field_traits, step = 1, n.draws = 64) 
#uses imputation
  #creates n.draw number of populations where the missing 
  #genotypes are filled in based on recombination frequencies

?calc.genoprob()
#just calculates the probabilities at different locations
field_traits <- calc.genoprob(field_traits, step = 1)
field_traits <- calc.genoprob(field_traits, error.prob=0.01)

summary(field_traits)
plot(field_traits)
names(field_data)

# grab marker names
# seed_so <- scanone(field_traits, method = "imp", pheno.col = 58, n.perm = 1000)
# plot(seed_so)
# summary(seed_so)

# seed_so_out <- scanone(field_traits, method = "imp", pheno.col = 58)
# plot(seed_so_out)

# seed_cim <- cim(field_traits, pheno.col = 58)
# plot(seed_cim)
# summary(seed_cim)


# frt_so <- scanone(field_traits, method = "imp", pheno.col = 57, n.perm = 1000)
# plot(frt_so)
# summary(frt_so)

# frt_so_out <- scanone(field_traits, method = "imp", pheno.col = 57)
# plot(frt_so_out)

frt_cim <- cim(field_traits, pheno.col = 57)
plot(frt_cim)
summary(frt_cim)
#              chr   pos   lod
# A03x5439663  A03  30.8 4.071
# A09x14496888 A09  68.9 5.590
# A10x11176662 A10  48.3 4.640

fit_so <- scanone(field_traits, method = "imp", pheno.col = 28)
plot(fit_so)
summary(fit_so)
fit_so_out <- scanone(field_traits, method = "imp", pheno.col = 28)
plot(fit_so_out)

fit_cim <- cim(field_traits, pheno.col = 28)
plot(fit_cim)
summary(fit_cim)[1,]
#              chr    pos   lod
# A01x26761091 A01  92.52 4.212
# A03x5068034  A03  28.95 4.686
# A09x18352859 A09  69.67 8.256



# fit_cr_so <- scanone(field_traits, method = "imp", pheno.col = 27)
# plot(fit_cr_so)
# summary(fit_cr_so)
# fit_cr_so_out <- scanone(field_traits, method = "imp", pheno.col = 27)
# plot(fit_cr_so_out)

# fit_cr_cim <- cim(field_traits, pheno.col = 27)
# plot(fit_cr_cim) # A10
# summary(fit)

field_QTLs <- list()
field_QTLs[[1]] <- makeqtl(field_traits, chr = "A03", pos = 30.8)
field_QTLs[[2]] <- makeqtl(field_traits, chr = "A09", pos = 69.67)
field_QTLs[[3]] <- makeqtl(field_traits, chr = "A10", pos = 48.3)
field_QTLs[[4]] <- makeqtl(field_traits, chr = "A01", pos = 92.52)
field_QTLs[[5]] <- makeqtl(field_traits, chr = "A03", pos = 30.8)
field_QTLs[[6]] <- makeqtl(field_traits, chr = "A09", pos = 69.67)

m1 <- c("A03x5439663", "A09x18352859", "A10x11176662")
m2 <- c("A01x26761091","A03x5439663", "A09x18352859")

# simulate phenotype data for expression for now
genotypes <- pull.geno(field_traits)
head(genotypes)
geno_names <- dimnames(genotypes)[[2]]

m3 <- sample(geno_names, 3, replace = FALSE)
m4 <- c("A09x18352859")

g31 <- genotypes[,m3[1]]; g32 <- genotypes[,m3[2]]; g33 <- genotypes[,m3[3]]; 
g41 <- genotypes[,m4[1]]

field_traits$pheno <- pull.pheno(field_traits)[-c("id")]
phenotypes <- pull.pheno(field_traits)
phenotypes
good_frt <- pull.pheno(field_traits)["good_frt_per_plant_UN2012"]
head(good_frt)
n.ind <- 124
expr <- runif(1,0.5,1)*good_frt + runif(2,0.5,1)[g31] + runif(2,0.5,1)[g41] + rnorm(n.ind)
expr
cor(expr, good_frt)


# there are some NAs in the phenotype matrix that I suspect is contributing to the problems.
# need to figure this out.




markers <- list(m1, m2)
names(markers) <- c("good_frt_per_plant_UN2012", "FitnessUN")
qdg_out <- qdg(cross= field_traits,
           phenotype.names     = c("good_frt_per_plant_UN2012","FitnessUN"),
           marker.names        = markers,
           QTL                 = field_QTLs,
           alpha               = 0.005,
           n.qdg.random.starts = 10,
           skel.method         = "pcskel")







