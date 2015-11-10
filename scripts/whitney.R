#######
###seems strange no gxe here
#######
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
fit1 <- eBayes(fit1)
toptable(fit1)
head(fit1)
head(fit1$design,20)

check <- as.data.frame(fit1$design)
check[1]

glmLRT(fit1)

summary(fit1)
toptable(fit1, coef = 1)
topTable(fit1, coef = 1:22)
topTable(fit1, coef = "Br_trtCR")

summary(decideTests(fit1)) #take a look at this output again
head(brass_voom)
?decideTests
decideTests(fit1)

load('~/git.repos/brassica_eqtl_v1.5/data/brass_group_fit_object.RData')
summary(brass_group_fit)
head(brass_group_fit$design)
group_design <- brass_group_fit$design
head(group_design)

?diffSplice
cont.matrix <- cbind(SvsUinWT=c(0,0,-2,-2),SvsUinMu=c(0,0,-2,2),Diff=c(0,0,0,4))
cont.matrix

treatment <- factor(sub("(Br_group)(\\d+)(_)(\\w+)",
                       "\\4", colnames(group_design)))
treatment

geno <- factor(sub("(Br_group)(\\d+)(_)(\\w+)",
                       "RIL_\\2", colnames(group_design)))
geno

trt <- model.matrix(~treatment)
colnames(trt)
trt2 <- trt[,2]
trt2


f
fit2 <- contrasts.fit(brass_group_fit, trt)


treatment
contrasts(treatment)

######sunday nov 8##########
load('~/git.repos/brassica_eqtl_v1.5/data/eQTL_full_model.RData')
ls()
design[,123]
colnames(design)
Br_RIL
design2 <- model.matrix(~0 + Br_trt:Br_RIL)
system.time(brass_voom2 <- voom(brassica_DE, design2, plot = FALSE))
system.time(fit2 <-lmFit(brass_voom2, design2)) # currently running full model
fit2 <- eBayes(fit2)
?topTable
topTable(fit2, number = 20000)

Br_trt
geno <- model.matrix(~0 +Br_RIL)
head(geno)
dim(geno)
geno[,112]
design3 <- model.matrix(~Br_trt:Br_RIL)
head(design3)
design3[,2]

system.time(brass_voom3 <- voom(brassica_DE, design3, plot = FALSE))
head(brass_voom3)
system.time(fit3 <-lmFit(brass_voom3, design3)) # currently running full model

head(fit1)
?exprs
?ns
design4 <- model.matrix(~Br_trt)
system.time(brass_voom4 <- voom(brassica_DE, design4, plot = FALSE))
head(brass_voom4)
system.time(fit4 <-lmFit(brass_voom4, design4)) # currently running full model

head(fit4)
head(Lambda)
lam_rows <- as.data.frame(rownames(Lambda))
lam_rows[1] <- as.character(lam_rows[1])
str(lam_rows)
colnames(lam_rows) <- paste("test")
lam_rows$test <- as.character(lam_rows$test)

fit4_coefs <- as.data.frame(fit4$coefficients)
head(rownames(fit4_coefs))
fit4_coefs$gene <- rownames(fit4_coefs)
head(fit4_coefs)

str(fit4_coefs)
fit4_red <- fit4_coefs[fit4_coefs$gene %in% lam_rows$test,]
dim(fit4_red)
fit4_red


write.table(fit4_red, "shade_genes_trt.csv", sep = ",")




colnames(fit1$coef )[123:244]
?topTable
gxe <- topTable(fit1, coef = 123:244, p.value = 0.05, number = 2000)
dim(gxe)
gxe
gxe_genes <- rownames(gxe)
gxe_genes
setwd('~/git.repos/brassica_eqtl_v1.5/data/')
?write.table
write.table(gxe_genes, "gxe_genes.csv", row.names = FALSE, sep = ",")




library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
load('~/git.repos/brassica_eqtl_v1.5/data/cr_un_eqtl.RData')
ls()

#takes a minute 
#GxE analysis
brassica_genes <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="br_blues_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brassica_genes)

class(brassica_genes)[1] <- "riself"
brassica_genes <- jittermap(brassica_genes)
str(brassica_genes)

#replace names to be compatable with eQTL package
names(brassica_genes$geno) <- paste(1:10)
names(brassica_genes$geno)

?read.table

gxegenes <- read.table("gxe_genes.csv", sep = ",", header = TRUE )
gxegenes
ls()
str(scanone.imp.1)
scanone_2 <- scanone.imp.1
head(scanone_2)[1:10]

scanone_2$marker <- rownames(scanone_2)
dim(scanone_2)
scanone_2 <- scanone_2[c(35042,1:35041)]

# a few genes have multiple arabidopsis hits
gxe_qtl <- scanone_2[c(1,2,3, gxegenes$x)]
str(gxe_qtl)
dim(gxe_qtl)
head(gxe_qtl)
plot(gxe_qtl)

# gxe_qtl <- gxe_qtl[-c(4:7)]

library(reshape2)
library(ggplot2)
melt?
?melt
head(gxe_qtl)[1:10]
dim(gxe_qtl)

# do partial melts of the dataframes
# follows where the cis-eqtl are
library(reshape2)
library(ggplot2)
gxe_melt <- melt(gxe_qtl[c(1:100)], id = c("marker", "chr", "pos"))
head(gxe_melt)
colnames(gxe_qtl)
dim(gxe_melt)

peak <- max(gxe_melt$value)
gxe_plot <- ggplot(gxe_melt)
gxe_plot <- gxe_plot +  theme_bw() + geom_point(aes(x = pos, y = value), size = 2) +
                        # geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        # geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(chr ~ .) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
gxe_plot



#compare shade genes with gxe genes
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_shade <- read.delim("br_shade_genes.csv", header = TRUE, sep = ",")
head(br_shade)
dim(br_shade)

head(scanone_imp_tot)[1:10]
br_shade
gxegenes

merged <- merge(br_shade, gxegenes, by.x = "V1", by.y = "x")
merged


head(scanone_2)[1:10]
merged_qtl <- scanone_2[c(1,2,3, merged$V1)]
str(merged_qtl)
dim(merged_qtl)
head(merged_qtl)
plot(merged_qtl)


merged_melt <- melt(merged_qtl, id = c("marker", "chr", "pos"))
head(merged_melt)
colnames(merged_qtl)
dim(merged_melt)

peak <- max(merged_melt$value)
merged_plot <- ggplot(merged_melt)
merged_plot <- merged_plot +  theme_bw() + geom_point(aes(x = pos, y = value), size = 2) +
                        # geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        # geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(chr ~ .) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
merged_plot
