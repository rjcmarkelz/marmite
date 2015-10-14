##########
# Cody Markelz
# markelz@gmail.com
# voom/limma analysis pipeline for brassica RNAseq data
##########

# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

# large file (112 MB) takes about 1 minute to get into memory
mapped_counts <- read.delim("RIL_v1.5_mapping.tsv", header = TRUE, sep = "\t")
dim(mapped_counts)
# [1] 43151   843

colnames(mapped_counts)

#replace all NA values with 0 
mapped_counts[is.na(mapped_counts)] <- 0
head(mapped_counts)
tail(mapped_counts)

#remove first row
mapped_counts <- mapped_counts[-1,]
head(mapped_counts)[,1:10]

# use this cleaned table for analysis
write.table(mapped_counts, file="RIL_mapped_counts_cleaned.csv", sep=",") 

# move the first column to row names
mapped_counts[,1]
rownames(mapped_counts) <- mapped_counts[,1]
mapped_counts <- mapped_counts[,-1]
head(mapped_counts)
dim(mapped_counts)

#remove these columns from full dataset because either CR or UN are missing
# Br_trtCR:Br_RIL234 columns 340:343
# Br_trtCR:Br_RIL311 columns 556:558
# causing failures to converge below in lmFit()
# need to remove now and redue the design matrix
colnames(mapped_counts)
dim(mapped_counts)
mapped_counts <- mapped_counts[,-c(556:558)]
dim(mapped_counts)
mapped_counts <- mapped_counts[,-c(340:343)]
dim(mapped_counts)
colnames(mapped_counts)


#######
# design matrix
#######
samples <- names(mapped_counts)
samples

Br_group <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4", colnames(mapped_counts)))
Br_group2 <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4\\5\\7", colnames(mapped_counts)))
Br_RIL   <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2", colnames(mapped_counts)))
Br_trt <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\4", colnames(mapped_counts)))
Br_group 
Br_group2
Br_RIL  # 122 levels
length(Br_RIL)
Br_trt <- relevel(Br_trt, ref = "UN")
Br_trt

# examine some potential contrasts
test <- contrasts(Br_trt)
test <- contr.sum(Br_trt)
head(test)
?contr.sum
contr.sum(Br_trt)
contrasts(Br_trt) <- contr.sum(2)
contrasts(Br_RIL) <- contr.sum(122)


# full model
design <- model.matrix(~ Br_RIL*Br_trt)
colnames(design)
head(design)

# group model for eQTL
group_design <- model.matrix(~ 0 + Br_group)
head(group_design)

######
# edgeR and Limma
######

# one time and will take a while to install all dependencies
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

# load libs
library(edgeR)
library(limma)

# whitney
# made directory structure the same on whitney
load('~/git.repos/brassica_eqtl_v1.5/data/brassica_voom_whitney.RData')
Br_group <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4", colnames(mapped_counts)))
brassica_DE <- DGEList(counts = mapped_counts, group = Br_group)
brassica_DE <- DGEList(counts = mapped_counts, group = Br_group)
brassica_DE$samples
dim(brassica_DE)
# [1] 43150   835

# keep genes with at least 1 count per million in at least 20 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
#[1] 35039  835

brassica_DE <- calcNormFactors(brassica_DE)
system.time(brass_voom <- voom(brassica_DE, design, plot = FALSE))
system.time(fit1 <-lmFit(brass_voom, design)) # currently running full model


brassica_DE <- DGEList(counts = mapped_counts, group = Br_group)
brassica_DE$samples
dim(brassica_DE)
# [1] 43150   835

# keep genes with at least 1 count per million in at least 20 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
#[1] 35039  835

brassica_DE <- calcNormFactors(brassica_DE)

# output plots, they are much to large to fit into memory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/output")
png(file="Brassica_MDS.png",width=1000,height=1000,res=100)
plotMDS(brassica_DE)
dev.off()

system.time(brass_voom <- voom(brassica_DE, group_design, plot = TRUE))
# Coefficients not estimable: ril234:trtCR ril311:trtCR 
# removed these above form analysis because causes issues in convergence

system.time(fit1 <-lmFit(brass_voom, group_design))
 #    user   system  elapsed 
 # 991.949  106.047 1092.513
fit1 <- eBayes(fit1)
toptable(fit1)
head(fit1)

?topTable

toptable(fit1, coef = 3)
topTable(fit1, coef = 1:20)

head(brass_voom)


#check out coeffs
brassica_BLUEs <- fit1$coeff
head(brassica_BLUEs)
str(brassica_BLUEs)
brassica_BLUE_df <- as.data.frame(fit1$coeff)
head(brassica_BLUE_df)


plot(brassica_BLUE_df)[1,]
hist(brassica_BLUE_df[,3])

brassica_BLUE_df_t <- as.data.frame(t(brassica_BLUE_df))
head(brassica_BLUE_df_t)[,1:10]

# examine a few gene expression values in linear to see how they are not normally
# distributed, column 26 as an example
hist(brassica_BLUE_df_t[,26])
test <- 2^(brassica_BLUE_df_t[,26])
hist(test)

# 2015_1_6
# I am convinced that I should leave these values in log2 for eQTL mapping

#export BLUES for eQTL mapping
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
?write.table
write.table(brassica_BLUE_df, file="brassica_blues.csv", sep=",",
           row.names = TRUE, col.names = TRUE) 

# an alternative way to do the design but I wanted to use a classical interpretation of
# the design matrix where the intercept is the mean value for each gene
# sample_info <- data.frame(sample = samples,
#                           ril = sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+", "\\2", samples),
#                           trt = sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+", "\\4", samples),
#                           rep = sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+", "\\7", samples),
#                           group = sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+", "\\2\\3\\4", samples)
# )
# head(sample_info)
# sample_info$trt <- relevel(sample_info$trt, ref = "UN")
#write.csv(sample_info, "sample_info.csv")

