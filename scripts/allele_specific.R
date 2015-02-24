##########
# Cody Markelz
# markelz@gmail.com
# calculate allele specific expression
# 2015_02_19
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
samples <- as.data.frame(names(mapped_counts))
samples

Br_group  <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4", colnames(mapped_counts)))
Br_group2 <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4\\5\\7", colnames(mapped_counts)))
Br_RIL    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2", colnames(mapped_counts)))
Br_trt    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\4", colnames(mapped_counts)))
Br_rep    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\7", colnames(mapped_counts)))
Br_rep
Br_group 
Br_group2
Br_RIL  # 122 levels
length(Br_RIL)
Br_trt <- relevel(Br_trt, ref = "UN")
Br_trt

# full model
design <- model.matrix(~ 0 + Br_RIL*Br_trt)
colnames(design)
head(design)

# group model for eQTL
group_design <- model.matrix(~ 0 + Br_group)
head(group_design)
colnames(group_design)


# edgeR and Limma
library(edgeR)
library(limma)

brassica_DE <- DGEList(counts = mapped_counts, group = Br_group)
brassica_DE$samples

dim(brassica_DE)
# [1] 43150   835

# keep genes with at least 1 count per million in at least 20 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
#[1] 35039   835

brassica_DE <- calcNormFactors(brassica_DE)
brassica_DE
str(brassica_DE)
# output plots, they are much to large to fit into memory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/output")
png(file = "Brassica_MDS_allele_specific.png", width = 1000, height = 1000, res = 100)
plotMDS(brassica_DE)
dev.off()

system.time(brass_voom <- voom(brassica_DE, design, plot = TRUE))
# Coefficients not estimable: ril234:trtCR ril311:trtCR 
# removed these above from analysis because causes issues in convergence
str(brass_voom)

?lmFit
?duplicateCorrelation
# look for correlation between replicates
# must account for in linear model
system.time(brass_dup <- duplicateCorrelation(brass_voom, design = design, block = Br_rep))
#     user   system  elapsed 
# 1509.105  114.936 1614.638

str(brass_dup)
brass_dup$consensus.correlation

# this should take a while
system.time(brass_fit <- lmFit(brass_voom, block = Br_rep,
             design = design, correlation=brass_dup$consensus.correlation))












