##########
# Cody Markelz
# markelz@gmail.com
# double check contrast statements 
# this is the same pipeline as allele specific
# 2015_02_25
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

mapped_counts_sub <- mapped_counts[,1:51]

#######
# design matrix
#######
samples <- as.data.frame(names(mapped_counts_sub))
samples

Br_group  <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4", colnames(mapped_counts)))
Br_group2 <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4\\5\\7", colnames(mapped_counts_sub)))
Br_RIL    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2", colnames(mapped_counts_sub)))
Br_trt    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\4", colnames(mapped_counts)))
Br_rep    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\7", colnames(mapped_counts_sub)))
Br_rep
Br_group 
Br_group2
Br_RIL  
length(Br_RIL)
Br_trt <- relevel(Br_trt, ref = "UN")
Br_trt

# full model
design <- model.matrix(~ 0 + Br_RIL*Br_trt)
colnames(design)
head(design)

trt_design <- model.matrix(~ 0 + Br_trt)
colnames(trt_design)
head(trt_design)
dim(trt_design)

# group model for eQTL
group_design <- model.matrix(~ 0 + Br_group)


# edgeR and Limma
library(edgeR)
library(limma)


##############
# group model for eQTL
##############
system.time(brass_voom_group <- voom(brassica_DE, group_design, plot = TRUE))

system.time(brass_group_fit <- lmFit(brass_voom_group, design = group_design))


summary(brass_group_fit)
brass_group_fit <- eBayes(brass_group_fit)
topTable(brass_group_fit)
TT_brass_group <- topTable(brass_group_fit)
head(TT_brass_group)[,1:5]
summary(brass_group_fit)
dim(brass_group_fit$coefficients)
# [1] 35039   244

colnames(brass_group_fit$coefficients)
br_un_cr    <- factor(sub("(Br_)(group)(\\d+)(_)(\\w+)",
                       "\\5", colnames(brass_group_fit$coefficients)))
br_un_cr_num <- as.numeric(br_un_cr)
str(br_un_cr_num)
br_un_cr_num <- replace(br_un_cr_num, br_un_cr_num == 2, 0)
?makeContrasts

group_fit_trt <- contrasts.fit(brass_group_fit, br_un_cr_num)
group_fit_trt <- eBayes(group_fit_trt)
topTable(group_fit_trt, number = 20)

gft_treat <- treat(group_fit_trt, lfc = 2)
topTreat(gft_treat)

#########
# treatment model
#########
system.time(brass_v_trt <- voom(brassica_DE, trt_design, plot = TRUE))
# Coefficients not estimable: ril234:trtCR ril311:trtCR 
# removed these above from analysis because causes issues in convergence
str(brass_v_trt)

system.time(brass_trt_fit <- lmFit(brass_v_trt, design = trt_design))
brass_trt_fit <- eBayes(brass_trt_fit)
topTable(brass_trt_fit)
brass_trt_fit$coefficients[1:2,1:2]

trt_design
contrast.matrix <- makeContrasts(Br_trtUN - Br_trtCR, levels = trt_design)
contrast.matrix

trt_cont <- contrasts.fit(brass_trt_fit, contrast.matrix)
trt_cont <- eBayes(trt_cont)
topTable(trt_cont, coef=1, adjust="BH")


##########
# whole dataset treatment model
##########
Br_group  <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2\\3\\4", colnames(mapped_counts)))
Br_RIL    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2", colnames(mapped_counts_sub)))
Br_trt    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\4", colnames(mapped_counts)))
Br_rep    <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\7", colnames(mapped_counts_sub)))
Br_rep
Br_group 
Br_group2
Br_RIL  
length(Br_RIL)
Br_trt <- relevel(Br_trt, ref = "UN")
Br_trt

# full model
design <- model.matrix(~ 0 + Br_RIL*Br_trt)
colnames(design)
head(design)

trt_design <- model.matrix(~ 0 + Br_trt)
colnames(trt_design)
head(trt_design)
dim(trt_design)

brassica_DE <- DGEList(counts = mapped_counts, group = Br_group)
brassica_DE$samples

dim(brassica_DE)
# [1] 43150    835

# keep genes with at least 1 count per million in at least 20 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
# [1] 35039   835

brassica_DE <- calcNormFactors(brassica_DE)
brassica_DE
str(brassica_DE)

system.time(brass_voom <- voom(brassica_DE, trt_design, plot = TRUE))
str(brass_voom)

?lmFit
# ?duplicateCorrelation
# look for correlation between replicates
# must account for in linear model
# system.time(brass_dup <- duplicateCorrelation(brass_voom, design = design, block = Br_rep))
# str(brass_dup)
# brass_dup$consensus.correlation

system.time(brass_fit <- lmFit(brass_voom, design = trt_design))

brass_fit <- eBayes(brass_fit)
topTable(brass_fit)
brass_fit$coefficients[1:2,1:2]

trt_design
contrast.matrix <- makeContrasts(Br_trtUN - Br_trtCR, levels = trt_design)
contrast.matrix

trt_cont <- contrasts.fit(brass_fit, contrast.matrix)
trt_cont <- eBayes(trt_cont)
topTable(trt_cont, coef = 1, adjust = "BH", number = 100)

#########
#### go enrichment cr vs un
#########
library(Biostrings)
library(GO.db)
library(goseq)

trt_cont_go <- topTable(trt_cont, coef = 1, adjust = "BH", number = Inf)
head(trt_cont_go)
dim(trt_cont_go)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
?write.csv
write.table(trt_cont_go, "un_vs_cr_DE_nodupcorr_model.csv", sep = ",", col.names = TRUE, row.names = TRUE)


# infile and manipulate go annotation file
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
brass_go <- read.table("Brassica_rapa_v1.5_final_annot.wego", header = FALSE)
head(brass_go)
colnames(brass_go) <- c("Gene", "GO")
dim(brass_go)
# make list object to store data
brass_go_list <- strsplit(as.character(brass_go[,2]),split=",",fixed=T)
head(brass_go_list)
names(brass_go_list) <- as.character(brass_go[,1])
head(brass_go_list)
length(brass_go_list)

# also need to infile the gene length data as goseq uses this to estimate any biases based on gene
# lengths in the dataset. 
brass_gene_lengths <- read.table("Brassica_rapa_v1.5_final_gene_lengths", header = FALSE)
head(brass_gene_lengths)
colnames(brass_gene_lengths) <- c("Gene", "length")
dim(brass_gene_lengths)
str(brass_gene_lengths)
# [1] 43463     2
brass_gene_lengths$Gene <- as.character(brass_gene_lengths$Gene)

# subset brass_gene_lengths for signifcant genes that were actually tested in the linear model
trt_cont_go$Gene <- rownames(trt_cont_go)
keeps <- as.character(trt_cont_go$Gene)
head(keeps)
head(trt_cont_go)
str(trt_cont_go)
brass_gene_lengths_red <- brass_gene_lengths[brass_gene_lengths$Gene %in% keeps,]
dim(brass_gene_lengths_red)
# [1] 35039     2

dim(trt_cont_go)
# [1] 35039     7
# conservative subset
brass_genes_sig <- subset(trt_cont_go, adj.P.Val < 0.0001)
dim(brass_genes_sig)
# [1] 5257    7

?nullp




