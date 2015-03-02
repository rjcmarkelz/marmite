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
             design = design, correlation = brass_dup$consensus.correlation))

summary(brass_fit)
brass_fit <- eBayes(brass_fit)
topTable(brass_fit)
TT_brass_fit <- topTable(brass_fit)
head(TT_brass_fit)[,1:5]


# group model for eQTL
system.time(brass_voom_group <- voom(brassica_DE, group_design, plot = TRUE))

system.time(brass_group_dup <- duplicateCorrelation(brass_voom, design = group_design, block = Br_rep))
    # user   system  elapsed 
# 1406.919  120.018 1518.139 

str(brass_dup)
brass_dup$consensus.correlation

# this should take a while
system.time(brass_group_fit <- lmFit(brass_voom_group, block = Br_rep,
             design = group_design, correlation = brass_group_dup$consensus.correlation))

#     user   system  elapsed 
# 8855.623  544.291 9359.772 

summary(brass_group_fit)
brass_group_fit <- eBayes(brass_group_fit)
?topTable
TT_brass_group <- topTable(brass_group_fit)
head(TT_brass_group)[,1:5]
summary(brass_group_fit)
dim(brass_group_fit$coefficients)
# [1] 35039   244
# write the coefficients to file

brass_group_coef <- brass_group_fit$coefficients
dim(brass_group_coef)
head(brass_group_coef)

write.table()

# output model fit that takes into account correlation between replicates of the same 
# RIL 3-5 reps per genotype.
write.table(brass_group_coef, "brassica_mixed_group_fit_coef.csv", row.names = TRUE, sep = ",")
#save output object 
save(brass_group_fit, file = "brass_group_fit_object.RData")
# also save the tstats
brass_group_tstats <- brass_group_fit$t
head(brass_group_tstats)
write.table(brass_group_tstats, "brassica_mixed_group_fit_tstats.csv", row.names = TRUE, sep = ",")



colnames(brass_group_fit$coefficients)
br_un_cr    <- factor(sub("(Br_)(group)(\\d+)(_)(\\w+)",
                       "\\5", colnames(brass_group_fit$coefficients)))
br_un_cr
br_un_cr_num <- as.numeric(br_un_cr)
str(br_un_cr_num)
br_un_cr_num <- replace(br_un_cr_num, br_un_cr_num == -1, 0)
?makeContrasts


group_fit_trt <- contrasts.fit(brass_group_fit, br_un_cr_num)
group_fit_trt <- eBayes(group_fit_trt)
topTable(group_fit_trt, number = 20)
str(group_fit_trt)
summary(group_fit_trt)
topTable(brass_group_fit, levels = group_design, )
brass_group_fit$coefficients["Bra011398",1:10]

# sanity check
colnames(mapped_counts)
brass_subset <- mapped_counts[,1:51]
head(brass_subset)

#FINALLY SOME ASE
# pull in contrast matrix that has the 1 or 0 based on genotype of each gene in each RIL
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Input")
gene_contrasts <- read.table("gene_marker_contrast_matrix_long.csv", sep = ",", header = TRUE)
head(gene_contrasts)

#order by gene name to make small test
gene_contrasts <- gene_contrasts[order(gene_contrasts$tx_name),]
head(gene_contrasts)

gene_cont_sub <- head(gene_contrasts, 10)
gene_cont_sub
dim(gene_cont_sub)

# order by gene name to get small test working
head(brass_group_coef)
brass_group_coef <- brass_group_coef[order(rownames(brass_group_coef)),]
brass_group_sub <- head(brass_group_coef, 10)

head(brass_group_sub)
head(gene_cont_sub)

rownames(gene_cont_sub) <- gene_cont_sub$tx_name

# subset for each treatment
brass_group_un <- as.data.frame(brass_group_sub[, grep("*_UN", colnames(brass_group_sub))])
brass_group_un 
str(brass_group_un)
brass_group_cr <- brass_group_sub[, grep("*_CR", colnames(brass_group_sub))]
brass_group_cr

colnames(brass_group_un) <- sub("(Br_group)(\\d+)(_)(UN)", "RIL_\\2", colnames(brass_group_un))
head(brass_group_un)
dim(brass_group_un)
dim(gene_cont_sub)

#remove tx_name column
gene_cont_sub <- gene_cont_sub[, -1]
dim(gene_cont_sub)

gene_cont_sub <- subset(gene_cont_sub, select = colnames(brass_group_un))
dim(gene_cont_sub)
head(gene_cont_sub)
head(brass_group_un)
library(genefilter)
# ?rowttests
# rowttests
# rowttests(brass_group_un, gene_cont_sub, tstatOnly = FALSE)
?tapply
str(brass_group_un)
brass_group_un <- as.data.frame(t(brass_group_un))
str(brass_group_un)
gene_cont_sub <- as.data.frame(t(gene_cont_sub))
str(gene_cont_sub)
str(brass_group_un[1,])
gene_cont_sub[2,]

gene_cont_sub[2,]
brass_group_un[1,]
gene_cont_sub$Bra000002
brass_group_un$Bra000002
brass_names <- names(brass_group_un)
str(brass_names)
for()

brass_names <- paste("Bra000002")
brass_names
brass_group_un
gene_cont_sub
tapply(brass_group_un$Bra000010, INDEX = gene_cont_sub$Bra000010,  FUN = mean)

length(brass_names)
brass_names
brass_names <- brass_names[1:3]
brass_names[2]

for (j in brass_names){
    print(j)
    test <- tapply(brass_group_un[,j], INDEX = gene_cont_sub[,j],  FUN = mean)
    print(test)
}


	
    test
test
test2

tapply(brass_group_un[1,], INDEX = gene_cont_sub[2,],  FUN = mean)
?mean


x  <- matrix(runif(40), nrow=4, ncol=10)
x
f2 <- factor(floor(runif(ncol(x))*2))
f2
f4 <- factor(floor(runif(ncol(x))*4))
f4
r1 <- rowttests(x)
r2 <- rowttests(x, f2)
r4 <- rowFtests(x, f4)
r1
r2
r4
#issues
# genes not represented in each row of the matrix and the coef df


t.test_results <- mapply(t.test, x= df1_t, y = df2_t, SIMPLIFY = F)





