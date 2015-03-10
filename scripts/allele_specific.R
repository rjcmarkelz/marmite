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

#################
#################
#FINALLY SOME ASE
# pull in contrast matrix that has the 1 or 0 based on genotype of each gene in each RIL
#################
#################
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Input")
gene_contrasts <- read.table("gene_marker_contrast_matrix_long.csv", sep = ",", header = TRUE)
head(gene_contrasts)

#order by gene name to make small test
gene_contrasts <- gene_contrasts[order(gene_contrasts$tx_name),]
head(gene_contrasts)



# order by gene name to get small test working
head(brass_group_coef)
dim(brass_group_coef)

# small subset for testing
#brass_group_coef <- brass_group_coef[order(rownames(brass_group_coef)),]
#brass_group_sub <- head(brass_group_coef, 10)
# gene_cont_sub <- head(gene_contrasts, 10)
# gene_cont_sub
# dim(gene_cont_sub)

head(brass_group_coef)
head(gene_contrasts)

# subset for each treatment
brass_group_un <- as.data.frame(brass_group_coef[, grep("*_UN", colnames(brass_group_coef))])
head(brass_group_un)
str(brass_group_un)
brass_group_cr <- brass_group_coef[, grep("*_CR", colnames(brass_group_coef))]
head(brass_group_cr)

#########
# write files
#########
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(brass_group_un, "brass_group_uncrowded_coef.csv", sep = ",", col.names = TRUE, row.names = TRUE)
write.table(brass_group_cr, "brass_group_crowded_coef.csv", sep = ",", col.names = TRUE, row.names = TRUE)


#rename columns
colnames(brass_group_un) <- sub("(Br_group)(\\d+)(_)(UN)", "RIL_\\2", colnames(brass_group_un))
head(brass_group_un)
dim(brass_group_un)
dim(gene_contrasts)

# make first column the rownames
rownames(gene_contrasts) <- gene_contrasts$tx_name

#remove tx_name column
gene_contrasts <- gene_contrasts[, -1]
dim(gene_contrasts)

# make sure both have the same RILs
ril_names <- intersect(colnames(brass_group_un), colnames(gene_contrasts))
gene_contrasts <- subset(gene_contrasts, select = ril_names)
brass_group_un <- subset(brass_group_un, select = ril_names)


dim(gene_contrasts)
head(gene_contrasts)
head(brass_group_un)
dim(brass_group_un)

#transpose to apply column wide functions
brass_group_un <- as.data.frame(t(brass_group_un))
str(brass_group_un)
gene_contrasts <- as.data.frame(t(gene_contrasts))
str(gene_contrasts)

# only compared genes in both sets
colnames(brass_group_un)
colnames(gene_contrasts)
brass_names <- intersect(colnames(brass_group_un), colnames(gene_contrasts))

# vector of common gene names
length(brass_names)
# 20851

# make some empty lists
# might take a while to run as code is not optimized for speed for this one time task
mylist <- list()
brass_pvalues <- list()
brass_tvalues <- list()

for (i in brass_names){
    #calculate means
    exp_means <- as.data.frame(tapply(brass_group_un[,i], INDEX = list(gene_contrasts[,i]),  FUN = mean))
    names(exp_means) <- i
    mylist <- c(mylist, exp_means)
    #t.test for pvalues
    brassfac <- as.factor(gene_contrasts[,i])
    testfacout <- t.test(brass_group_un[,i] ~ brassfac)
    out_p <- testfacout$p.value
    names(out_p) <- i
    brass_pvalues <- c(brass_pvalues,out_p)
    out_t <- testfacout$statistic
    names(out_t) <- i
    brass_tvalues <- c(brass_tvalues, out_t)
    
}

#make dataframes from list outputs
head(mylist)
mlistdf <- data.frame(t(sapply(mylist, c)))
brass_p_df <- data.frame(sapply(brass_pvalues, c))
brass_t_df <- data.frame(sapply(brass_tvalues, c))
head(mlistdf)
head(brass_t_df)
head(brass_p_df)
#double check before merge
rownames(mlistdf)
rownames(brass_p_df)

#merge dataframes based on row name
brass_merge <- merge(mlistdf, brass_t_df, by = "row.names")
brass_merge
brass_merge <- merge(brass_merge, brass_p_df, by.x = "Row.names", by.y = "row.names")
brass_merge
colnames(brass_merge) <- c("gene_name", "R500", "IMB211", "t_stat", "p_value")

# adjust for multiple testing
brass_merge$p_adjusted <- p.adjust(brass_merge$p_value, method = "BY", n = length(brass_merge$p_value))
head(brass_merge)

#print dataframe to .csv file
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(brass_merge, "allele_specific_test_p_adjusted.csv", sep = ",", col.names = TRUE, row.names = FALSE)


# quick data subset for brassica UV QTL
setwd("/Users/Cody_2/git.repos/brassica_UV/data/")
#redue for marc to get wide format
brass_UV_7 <- read.table("gr7_3_wide_candidates.txt", sep = " ")
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
brass_merge <- read.table("allele_specific_test_p_adjusted.csv", sep = ",", header = TRUE)
head(brass_UV_7)
dim(brass_UV_7)
# [1] 131   1
str(brass_UV_7)
brass_UV_7 <- as.character(brass_UV_7$x)

brass_UV_8 <- read.table("gr8_2_candidates.txt", sep = " ")
head(brass_UV_8)
dim(brass_UV_8)
# [1] 171   1
brass_UV_8 <- as.character(brass_UV_8$x)

brass_UV_7_exp <- brass_merge[brass_merge$gene_name %in% brass_UV_7,]
head(brass_UV_7_exp)
dim(brass_UV_7_exp)
# [1] 106   6
write.table(brass_UV_7_exp, "brass_UV_A07_peak_genes.csv", sep = ",", col.names = TRUE, row.names = FALSE)


brass_UV_8_exp <- brass_merge[brass_merge$gene_name %in% brass_UV_8,]
head(brass_UV_8_exp)
dim(brass_UV_8_exp)
# [1] 102   6
write.table(brass_UV_8_exp, "brass_UV_A08_peak_genes.csv", sep = ",", col.names = TRUE, row.names = FALSE)



