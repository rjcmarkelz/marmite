
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
field_data_t

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
str(field_traits)

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
str(field_traits)

# subset to remove RILS with many NA's
field_traits$pheno$ID <- field_traits$pheno$id
head(field_traits$pheno)

field_traits <- subset(field_traits, ind = c("-RIL_104", "-RIL_113", "-RIL_123", "-RIL_294", "-RIL_328", "-RIL_329" ))
head(field_traits$pheno)

summary(field_traits)
field_traits <- sim.geno(field_traits, step = 1, n.draws = 64) 

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

frt_so_out <- scanone(field_traits, method = "imp", pheno.col = 57)
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



"A08x19238852"

m1 <- c("A03x5439663", "A09x18352859", "A10x11176662")
m2 <- c("A01x26761091","A03x5439663", "A09x18352859")

# simulate phenotype data for expression for now
?fill.geno
#note caution with this so I do not have to throw out those indivduals or markers
field_traits <- fill.geno(field_traits, method = "argmax")
field_traits
summary(field_traits)
genotypes <- pull.geno(field_traits)
head(genotypes)
geno_names <- dimnames(genotypes)[[2]]

set.seed(12345)
m3 <- sample(geno_names, 3, replace = FALSE)
m4 <- c("A09x18352859")

g31 <- genotypes[,m3[1]]; g32 <- genotypes[,m3[2]]; g33 <- genotypes[,m3[3]]; 
g41 <- genotypes[,m4[1]]
g41


good_frt <- pull.pheno(field_traits)["good_frt_per_plant_UN2012"]
plot(good_frt$good_frt_per_plant_UN2012)

fit_un <- pull.pheno(field_traits)["FitnessUN"]

#make copy of cross to play with 
field_traits_red <- field_traits
summary(field_traits_red)
summary(field_traits)


n.ind <- 119


# 0.24 corr
set.seed(123)
rnorm(2,0.2,1)*good_frt
expr <- runif(2,.001,.1)*good_frt + runif(2,.001,.1)*fit_un + runif(2,50,100)[g41] 
expr
cor(expr, good_frt)
names(expr) <- paste("expr")
expr
cor(expr, good_frt)
cor(good_frt, fit_un)
cor(expr, fit_un)


# #.99 corr
# set.seed(12345)
# expr2 <- rnorm(2,0.2,1)*good_frt + rnorm(2,0.5,1)[g31] + rnorm(2,0.5,1)[g41]
# expr2
# cor(expr2, good_frt)

# # 0.6 corr
# set.seed(123456)
# expr3 <- rnorm(2,0.2,1)*good_frt + rnorm(2,0.5,1)[g31] + rnorm(2,0.5,1)[g41]
# expr3
# cor(expr3, good_frt)

# # -0.69 corr
# set.seed(123)
# expr3 <- rnorm(2,0.2,1)*good_frt + rnorm(2,0.5,1)[g31] + rnorm(2,0.5,1)[g41]
# expr3
# cor(expr3, good_frt)

# there are some NAs in the phenotype matrix that I suspect is contributing to the problems.
# need to figure this out.

markers <- list(m1, m2, m4)
markers

field_traits_red <- sim.geno(field_traits_red, step = 1, n.draws = 64) 

field_QTLs <- list()
m1.pos <- find.markerpos(field_traits_red, m1)
m1.pos
field_QTLs[[1]] <- makeqtl(field_traits_red, chr = m1.pos[,"chr"], pos = m1.pos[,"pos"])
m2.pos <- find.markerpos(field_traits_red, m2)
m2.pos
field_QTLs[[2]] <- makeqtl(field_traits_red, chr = m2.pos[,"chr"], pos = m2.pos[,"pos"])
m4.pos <- find.markerpos(field_traits_red, m4)
m4.pos
field_QTLs[[3]] <- makeqtl(field_traits_red, chr = m4.pos[,"chr"], pos = m4.pos[,"pos"])
names(field_QTLs) <- c("good_frt_per_plant_UN2012", "FitnessUN", "expr")
field_QTLs

field_traits_red$pheno <- data.frame(good_frt, fit_un, expr)
field_traits_red$pheno
plot(cim(field_traits_red, pheno.col = 3))
plot(cim(field_traits_red, pheno.col = 1))
plot(cim(field_traits_red, pheno.col = 2))

names(markers) <- c("good_frt_per_plant_UN2012", "FitnessUN", "expr")

field_traits_red <- sim.geno(field_traits_red, step = 1, n.draws = 64) 
summary(field_traits_red)

qdg_out <- qdg(cross= field_traits_red,
           phenotype.names     = c("good_frt_per_plant_UN2012","FitnessUN", "expr"),
           marker.names        = markers,
           QTL                 = field_QTLs,
           alpha               = 0.005,
           n.qdg.random.starts = 10,
           skel.method         = "pcskel")

summary(qdg_out)

graph2 <- graph.qdg(qdg_out)
graph2
plot(graph2)
tkplot(graph2)










#################
# aggregated traits for tiffany
#################
#################
setwd("~/git.repos/brassica_meta_analysis/Output/")                
aggregated <- read.table("brassica_blups.csv", sep = ",", header = TRUE) 
head(aggregated)
dim(aggregated)

aggregated$RILs

aggregated_t <- as.data.frame(t(aggregated))
aggregated_t
head(aggregated_t)
colnames(aggregated_t)
dim(aggregated_t)
aggregated_t[76,] <- aggregated_t[1,]
rownames(aggregated_t)[76] <- "id"
rownames(aggregated_t)
aggregated_t <- aggregated_t[-1,]
head(aggregated_t)
tail(aggregated_t)
rownames(aggregated_t)

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
write.table(aggregated_t, "aggregated_RQTL.csv", sep = ",", col.names = FALSE)


library(qtl)
library(qtlnet)
aggregated_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="aggregated_RQTL.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(aggregated_traits)
head(aggregated_traits)

class(aggregated_traits)[1] <- "riself"
aggregated_traits <- jittermap(aggregated_traits)
aggregated_traits


aggregated_traits <- sim.geno(aggregated_traits, step = 1, n.draws = 64) 
plot(aggregated_traits)
names(aggregated_traits$pheno)

leafarea_cim <- cim(aggregated_traits, pheno.col = 18)
plot(leafarea_cim)
setwd("~/git.repos/brassica_meta_analysis/Output/")
write.table(leafarea_cim, "specific_leaf_area.csv", sep = ",", col.names = TRUE, row.names = TRUE)

leaflength_cim <- cim(aggregated_traits, pheno.col = 37)
plot(leaflength_cim)
setwd("~/git.repos/brassica_meta_analysis/Output/")
write.table(leaflength_cim, "leaf_length.csv", sep = ",", col.names = TRUE, row.names = TRUE)

biomass_cim <- cim(aggregated_traits, pheno.col = 22)
plot(biomass_cim)
setwd("~/git.repos/brassica_meta_analysis/Output/")
write.table(biomass_cim, "biomass.csv", sep = ",", col.names = TRUE, row.names = TRUE)

period18C_cim <- cim(aggregated_traits, pheno.col = 1)
plot(period18C_cim)
setwd("~/git.repos/brassica_meta_analysis/Output/")
write.table(period18C_cim, "period18C.csv", sep = ",", col.names = TRUE, row.names = TRUE)

height_cim <- cim(aggregated_traits, pheno.col = 13)
plot(height_cim)
setwd("~/git.repos/brassica_meta_analysis/Output/")
write.table(height_cim, "height.csv", sep = ",", col.names = TRUE, row.names = TRUE)


biomass_so <- scanone(aggregated_traits, method = "imp", pheno.col = 22, n.perm = 1000)
plot(biomass_so)






