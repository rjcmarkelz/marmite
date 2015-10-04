########
# NEXT
########

#load eQTL data
summary(brassica_genes)
# get intervals for yield QTL
# subset genes in interval
# Build large model with all genes and yield QTLs


#####
# Data integration examples
#####

#################
#################
# aggregated traits
#################
#################
setwd("~/git.repos/brassica_meta_analysis/Output/")                
aggregated <- read.table("brassica_blups.csv", sep = ",", header = TRUE) 
head(aggregated)
dim(aggregated)

#################
#################
# branching
#################
#################
setwd("/Users/Cody_2/git.repos/brassica_branching/raw_data")

brass_branch <- read.cross("csv", file = 
	                  "Brassica_F8_v2.1_Brassica2010Mapping_Rubin_for_Cody.csv")

class(brass_branch)[1]  <-  "riself"
brass_branch <- jittermap(brass_branch)
head(brass_branch$pheno)
head(brass_branch$geno)
brass_branch <- calc.genoprob(brass_branch, step = 1, error.prob = 0.001)
branch <- brass_branch$pheno
head(branch)
dim(branch)
brass_geno <- pull.geno(brass_branch)
head(brass_geno)[, 1:10]
str(brass_branch)

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
head(field_data)

############
# create RQTL dataset
############
head(field_data)
field_data$Line
# remove parental data in first two rows
# field_data <- field_data[-c(1:2),]

# merge two large dataframes
head(aggregated)
head(field_data)
dim(aggregated)
dim(field_data)

field_data2 <- merge(field_data, aggregated, by.x = "Line", by.y = "RILs", all.x = TRUE)
dim(field_data2)
# [1] 125 135
head(field_data2)

field_data_t <- as.data.frame(t(field_data2))
field_data_t
head(field_data_t)
colnames(field_data_t)
dim(field_data_t)
field_data_t[135,] <- field_data_t[1,]
rownames(field_data_t)[135] <- "id"
rownames(field_data_t)
field_data_t <- field_data_t[-1,]
head(field_data_t)
tail(field_data_t)
field_data_t

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
write.table(field_data_t, "all_traits_RQTL.csv", sep = ",", col.names = FALSE)
write.table(field_data2, "all_traits.csv", sep = ",", col.names = TRUE)


#make sure directory is set
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")

library(qtl)
library(qtlbim)
field_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="all_traits_RQTL.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(field_traits)
plot(field_traits)

names(field_traits$pheno)
plot(cim(field_traits, pheno.col = 73))
plot((), ())
cor(field_traits$pheno[2], field_traits$pheno[118], use = "pairwise.complete.obs")
cor(field_traits$pheno[73], field_traits$pheno[74], use = "pairwise.complete.obs")
sec <- pull.pheno(field_traits, pheno.col = 74)
height <- pull.pheno(field_traits, pheno.col = 73)
hist(height)
hist(sec)
plot(height, sec)

cor(field_traits$pheno[57], field_traits$pheno[28], use = "pairwise.complete.obs")
# 0.8257112
cor(field_traits$pheno[82], field_traits$pheno[28], use = "pairwise.complete.obs")
cor(field_traits$pheno[82], field_traits$pheno[57], use = "pairwise.complete.obs")

cor(field_traits$pheno[86], field_traits$pheno[57], use = "pairwise.complete.obs")
cor(field_traits$pheno[86], field_traits$pheno[28], use = "pairwise.complete.obs")

cor(field_traits$pheno[73], field_traits$pheno[86], use = "pairwise.complete.obs")
cor(field_traits$pheno[2], field_traits$pheno[28], use = "pairwise.complete.obs")
cor(field_traits$pheno[4], field_traits$pheno[28], use = "pairwise.complete.obs")

#leaf traits
cor(field_traits$pheno[79], field_traits$pheno[77], use = "pairwise.complete.obs")
cor(field_traits$pheno, field_traits$pheno)









brassica_traits_qb <- field_traits
brassica_traits_qb <- qb.genoprob(brassica_traits_qb, step=2, stepwidth = "variable")
summary(brassica_traits_qb)

brassica_traits_qb_leaflength <- qb.mcmc(brassica_traits_qb, pheno.col = 2, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(brassica_traits_qb_leaflength)
str(brassica_traits_qb_leaflength)
plot(brassica_traits_qb_leaflength)
plot(qb.coda(brassica_traits_qb_leaflength))
best <- qb.best(brassica_traits_qb_leaflength)
best

?qb.scanone
temp <- qb.scanone(brassica_traits_qb_leaflength, type.scan = "posterior",  epistasis = FALSE,)
plot(temp)

brassica_traits_qb_leaflength2 <- qb.mcmc(brassica_traits_qb, pheno.col = 4, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(brassica_traits_qb_leaflength2)
str(brassica_traits_qb_leaflength2)
plot(brassica_traits_qb_leaflength2)
plot(qb.coda(brassica_traits_qb_leaflength2))
best <- qb.best(brassica_traits_qb_leaflength2)
best

?qb.scanone
temp2 <- qb.scanone(brassica_traits_qb_leaflength2, type.scan = "posterior",  epistasis = FALSE)
plot(temp2)



brassica_fitness_qb <- qb.mcmc(brassica_traits_qb, pheno.col = 28, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(brassica_fitness_qb)
str(brassica_fitness_qb)
plot(brassica_fitness_qb)
plot(qb.coda(brassica_fitness_qb))
best <- qb.best(brassica_fitness_qb)
best

?qb.scanone
temp3 <- qb.scanone(brassica_fitness_qb, type.scan = "posterior",  epistasis = FALSE)
plot(temp3)

brassica_fitness2_qb <- qb.mcmc(brassica_traits_qb, pheno.col = 57, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(brassica_fitness2_qb)
str(brassica_fitness2_qb)
plot(brassica_fitness2_qb)
plot(qb.coda(brassica_fitness2_qb))
best <- qb.best(brassica_fitness2_qb)
best

?qb.scanone
temp4 <- qb.scanone(brassica_fitness2_qb, type.scan = "posterior", epistasis = FALSE)
plot(temp4)

fit <- pull.pheno(brassica_traits_qb, pheno.col = 28)
fit2 <- pull.pheno(brassica_traits_qb, pheno.col = 57)
leafl1 <- pull.pheno(brassica_traits_qb, pheno.col = 2)
leafl2 <- pull.pheno(brassica_traits_qb, pheno.col = 4)
?cor()
cor(fit, leafl1, use = "pairwise.complete.obs")
cor(fit, leafl2, use = "pairwise.complete.obs")
cor(fit2, leafl1, use = "pairwise.complete.obs")
cor(fit2, leafl2, use = "pairwise.complete.obs")
plot(fit, fit2)
# library(qtlnet)

########
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

library(ggplot2)
library(qtl)


# this is a quick script to make some plots for the PAG poster.
#old new compare
brassica_traits <- read.cross("csvsr", genfile ="old_map_rqtl_missing_RILS_removed.csv", 
                         phefile="br_blups_RQTL.csv", genotypes=c("0","2"))

class(brassica_traits)[1] <- "riself"
brassica_traits <- jittermap(brassica_traits)
brassica_traits
traits_list <- pull.pheno(brassica_traits)
colnames(traits_list)

so.perm <- scanone(brassica_traits, method = "imp", n.perm = 1000) 
summary(so.perm)
so.perm95 <- summary(so.perm)[1] #keep 95%
# LOD thresholds (1000 permutations)
#      lod
# 5%  2.44
# 10% 2.13

X2011_STP
scanone_2011_STP <- scanone(brassica_traits, pheno.col = "X2011_STP", method = "imp", use="all.obs", chr = 3)
scanone_2011_STP

plot(scanone_2011_STP, chr = 3)


peak <- max(scanone_2011_STP$lod)
oldmapplot <- ggplot(scanone_2011_STP)
oldmapplot <- oldmapplot +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
oldmapplot
# new version of the map
br_phys <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="br_blups_RQTL.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
class(br_phys)[1] <- "riself"
br_phys
so.perm2 <- scanone(br_phys, method = "imp", n.perm = 1000) 
summary(so.perm2)
so.perm295 <- summary(so.perm2)[1] #keep 95%
# LOD thresholds (1000 permutations)
#      lod
# 5%  3.23
# 10% 2.81

?scanone
scanone_2011_STP <- scanone(br_phys, pheno.col = "X2011_STP", method = "imp", use="all.obs", chr = "A03")
plot(scanone_2011_STP, chr = "A03")

peak2 <- max(scanone_2011_STP$lod)
newplot_map <- ggplot(scanone_2011_STP)
newplot_map <- newplot_map +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 3.08, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
newplot_map

##############
##############
##############

library(qtlbim)
br_phys_qb2 <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="br_blups_RQTL.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
br_phys_qb2 <- jittermap(br_phys_qb2)
br_phys_qb <- qb.genoprob(br_phys_qb2, step=2, stepwidth = "variable")
summary(br_phys_qb)

?qb.mcmc
br_phy_height <- qb.mcmc(br_phys_qb, epistasis = TRUE, pheno.col = "X2011_STP", seed = 1616,  n.iter = 3000)
plot(qb.loci(br_phy_height))
?qb.scanone
height_post <- qb.scanone(br_phy_height, type.scan = "posterior",  epistasis = FALSE, smooth = 2)
height_post <- as.data.frame(height_post)
height_post <- subset(height_post, chr == "A03")
height_post
sum(height_post$main)
dim(height_post)

peak3 <- max(height_post$main)
heightbayes <- ggplot(height_post)
heightbayes <- heightbayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak3 * -0.02), yend = (peak3 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
heightbayes
###############
###############
###############

spec_traits <- read.cross("csvsr", genfile ="spec_gen.csv", 
                         phefile="spec_phe.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
summary(spec_traits)
spec_traits



spec_traits_qb <- qb.genoprob(spec_traits, step=2, stepwidth = "variable")
summary(spec_traits_qb)
str(spec_traits_qb)

spec_pri <- qb.mcmc(spec_traits_qb, epistasis = TRUE, pheno.col = 6, seed = 1616,  n.iter = 10000)
plot(spec_pri)
summary(spec_pri)
?qb.scanone

spec_pri_post <- qb.scanone(spec_pri, type.scan = "posterior",  epistasis = FALSE, smooth = 2)
spec_pri_post <- as.data.frame(spec_pri_post)
spec_pri_post <- subset(spec_pri_post, chr == "A10")
spec_pri_post
sum(spec_pri_post$main)
dim(spec_pri_post)

peak4 <- max(spec_pri_post$main)
spec_pribayes <- ggplot(spec_pri_post)
spec_pribayes <- spec_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
spec_pribayes

data(qbExample)
str(qbExample)
spec_pri_HDI <- qb.hpdone(spec_pri, effects = "estimate")
spec_pri_HDI
test
head(test)






# note ran into some map merge issues
# need to replace some values here and rerun for the spec_traits using the combined object
combined$pheno
spec_pri <- qb.mcmc(combined, epistasis = TRUE, pheno.col = 4, seed = 1616,  n.iter = 3000)
plot(spec_pri)
summary(spec_pri)
?qb.scanone

spec_pri_post <- qb.scanone(spec_pri, type.scan = "posterior",  epistasis = FALSE, smooth = 2)
spec_pri_post <- as.data.frame(spec_pri_post)
spec_pri_post <- subset(spec_pri_post, chr == "A10")
spec_pri_post
sum(spec_pri_post$main)
dim(spec_pri_post)

peak4 <- max(spec_pri_post$main)
spec_pribayes <- ggplot(spec_pri_post)
spec_pribayes <- spec_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
spec_pribayes

data(qbExample)
str(qbExample)
spec_pri_HDI <- qb.hpdone(spec_pri, effects = "estimate")
spec_pri_HDI
test
###############
###############
###############

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
#load data
load("eqtl_PAG.RData")
ls()

scanone_imp_p450 <- scanone(brass_total, pheno.col = "Bra009312", method = "imp", use="all.obs", chr = "A10")
plot(scanone_imp_p450)
scanone_imp_p450

peak <- max(scanone_imp_p450$lod)
peak 

p450_plot <- ggplot(scanone_imp_p450)
p450_plot <- p450_plot +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) + 
                        geom_hline(yintercept = 3.23, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20))
p450_plot
head(br_phys)



###bayes####
p450_traits_qb <- brass_total
class(p450_traits_qb)[1] <- "bc"

p450_traits_qb <- qb.genoprob(p450_traits_qb, step=2, stepwidth = "variable")
summary(p450_traits_qb)

summary()

p450_pri <- qb.mcmc(p450_traits_qb, epistasis = TRUE, pheno.col = "Bra009312", seed = 1616,  n.iter = 3000)
plot(p450_pri)
p450_pri$pairs
?qb.scanone

p450_pri_post <- qb.scanone(p450_pri, type.scan = "posterior", epistasis = FALSE, smooth = 2)
p450_pri_post <- as.data.frame(p450_pri_post)
p450_pri_post <- subset(p450_pri_post, chr == "A10")
p450_pri_post
sum(p450_pri_post$main)
dim(p450_pri_post)

peak4 <- max(p450_pri_post$main)
p450_pribayes <- ggplot(p450_pri_post)
p450_pribayes <- p450_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
p450_pribayes

str(brass_total_qb)

p450_HDI <- qb.hpdone(p450_pri, effects = "estimate")
p450_HDI

#################
#################
# reflectance gene candidates
#################
#################
# gene candidates from 95% bayes confidence interval ----------

setwd("/Users/Cody_2/git.repos/brassica_reflectance/")
specgenes <- read.csv("spectral_gene_candidates.csv", header = TRUE)
head(specgenes)
str(specgenes)
dim(specgenes)
# [1] 45  4
specgenes$Gene <- as.character(specgenes$Gene)


# QTL mapping on the genes in the reflectance QTL interval -----------
scanone_imp_p450_all <- scanone(brass_total, pheno.col = specgenes$Gene,
                         method = "imp", use = "all.obs")
#####Error in scanone(brass_total, pheno.col = specgenes$Gene, method = "imp", : 
# there are genes that are not expressed in our dataset.

# drop these non-expressed genes -----------
drops <- as.character(c("Bra009297", "Bra009304", "Bra1002808", "Bra1002809",
                        "Bra1002810", "Bra1002811", "Bra1002812",
                         "Bra1002813", "Bra1002814", "Bra009314", "Bra009324"))
drops
str(drops)
specgenes2 <- subset(specgenes, !(specgenes$Gene %in% drops))
specgenes2

# QTL mapping on the cleaned data -----------
scanone_imp_p450_all <- scanone(brass_total, pheno.col = specgenes2$Gene,
                           method = "imp", use="all.obs", chr = "A10")
head(scanone_imp_p450_all)
tail(scanone_imp_p450_all)

# plot a few LOD traces to get a feel for things -----------
plot(scanone_imp_p450_all$Bra009331)
plot(scanone_imp_p450_all$Bra009326)
plot(scanone_imp_p450_all$Bra009312)


#########
# other cis-eqtl
other_pri <- qb.mcmc(p450_traits_qb, epistasis = TRUE, pheno.col = "Bra009331", seed = 1616,  n.iter = 3000)
plot(other_pri)
other_pri$pairs
?qb.scanone

other_pri_post <- qb.scanone(other_pri, type.scan = "posterior", epistasis = FALSE, smooth = 2)
other_pri_post <- as.data.frame(other_pri_post)
other_pri_post <- subset(other_pri_post, chr == "A10")
other_pri_post
sum(other_pri_post$main)
dim(other_pri_post)

peak4 <- max(other_pri_post$main)
other_pribayes <- ggplot(other_pri_post)
other_pribayes <- other_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
other_pribayes

str(brass_total_qb)

other_HDI <- qb.hpdone(other_pri, effects = "estimate")
str(other_HDI)
other_HDI

###########
# QDG model with 95 HDI estimates
###########
# get closest markers to QTLs
# get phenotypes into one dataframe
# make the QTL 

spec_summary  <- as.data.frame(summary(spec_pri_HDI))
p450_summary  <- as.data.frame(summary(p450_HDI))
other_summary <- as.data.frame(summary(other_HDI))

p450_summary
?find.marker

mspec <- find.marker(combined, chr = spec_summary$chr, pos = spec_summary$pos)
m450 <-  find.marker(p450_traits_qb, chr = p450_summary$chr, pos = p450_summary$pos)
m_other <- find.marker(p450_traits_qb, chr = other_summary$chr, pos = other_summary$pos)
mspec
m_other
m450
markers <- list(mspec, m_other, m450)
markers
names(markers) <- c("npqi_UN", "Bra009331", "Bra009312")




?pull.pheno
head(pull.pheno(spec_traits_qb))

pull.pheno(p450_traits_qb)[,1:5]
?write.cross
write.cross(spec_traits_qb)
write.cross(p450_traits_qb, chr = "A10", format = "csvsr" )
summary(p450_traits_qb)
summary(spec_traits_qb)


spec_traits_qb$pheno$id
p450_traits_qb$pheno$id2 <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
p450_traits_qb$pheno$id2
`
combined <- p450_traits_qb
p450_red <- pull.pheno(p450_traits_qb, pheno.col = c("id", "id2", "Bra009331", "Bra009312"))
p450_red

spec_traits_qb$pheno
spec_red <- pull.pheno(spec_traits_qb, pheno.col = c(1, 6))
spec_red
merged_pheno <- merge(p450_red, spec_red, by.x = "id2", by.y = "id")
merged_pheno
merged_pheno <- merged_pheno[order(merged_pheno$id),]
merged_pheno$doublcheck <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
merged_pheno
combined$pheno <- merged_pheno[,c(2,5,3,4)]
combined$pheno


combined <- fill.geno(combined, method = "argmax")
combined

names(markers) <- c("npqi_UN", "Bra009331", "Bra009312")
mspec
m_other
m450


combined

combined <- sim.geno(combined, step = 1, n.draws = 64) 
m1.pos <- find.markerpos(combined, mspec)
m1.pos
field_QTLs <- list()
field_QTLs[[1]] <- makeqtl(combined, chr = m1.pos[,"chr"], pos = m1.pos[,"pos"])
m2.pos <- find.markerpos(combined, m_other)
m2.pos
field_QTLs[[2]] <- makeqtl(combined, chr = m2.pos[,"chr"], pos = m2.pos[,"pos"])
m4.pos <- find.markerpos(combined, m450)
m4.pos
field_QTLs[[3]] <- makeqtl(combined, chr = m4.pos[,"chr"], pos = m4.pos[,"pos"])
names(field_QTLs) <- c("npqi_UN", "Bra009331", "Bra009312")
field_QTLs
markers
plot(cim(combined, pheno.col = 3))
plot(cim(combined, pheno.col = 4))
plot(cim(combined, pheno.col = 2))


combined$pheno
subset
combined <- subset(combined, ind = c("-RIL_329","-RIL_311", "-RIL_294", "-RIL_234", "-RIL_104", "-RIL_113", "-RIL_123"))


##########
##########
# so close, but need to finish talk
##########
cor(combined$pheno$npqi_UN, combined$pheno$Bra009312)
cor(combined$pheno$npqi_UN, combined$pheno$Bra009331)
cor(combined$pheno$Bra009331, combined$pheno$Bra009312)


library(qtlnet)

qdg_out <- qdg(cross= combined,
           phenotype.names     = c( "Bra009331", "Bra009312", "npqi_UN"),
           marker.names        = markers,
           QTL                 = field_QTLs,
           alpha               = 0.0005,
           n.qdg.random.starts = 30,
           skel.method         = "pcskel")




?qdg
summary(qdg_out)
?graph.qdg
graph2 <- graph.qdg(qdg_out)
graph2
plot(graph2)
plot.igraph(graph2, edge.color="black")
tkplot(graph2, edge.color = "black")

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
save.image()
ls()




















