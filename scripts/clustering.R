# a majority of this code is from various tutorials 
# on the WGCNA site. 
library(WGCNA)

# install dependencies
source("http://bioconductor.org/biocLite.R") 
biocLite("impute") 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

brassica_blues <- read.delim("brassica_blues.csv", 
	                    header = TRUE, row.names = 1, sep = ",")
head(brassica_blues)

# uncrowded only
brassica_blues <- brassica_blues[,grepl("*\\UN$",names(brassica_blues))]
dim(brassica_blues)

# need to transpose for WGCNA and change to RIL names
brexp <- as.data.frame(t(brassica_blues))
head(brexp)[,1:10]

# add RIL numbers
row.names(brexp) <- sub("(Br_group)(\\d+)(_)(UN)", "RIL_\\2", row.names(brexp))
head(brexp)[,1:10]
dim(brexp)

# random sample of values for testing
set.seed(1257)
?sample
sampvec <- sample(1:35039, 5000, replace = FALSE) 
sampvec
brexp_rand <- brexp[,sampvec]
dim(brexp_rand)
head(brexp_rand)[,1:10]


# WGCNA analysis pipeline begin
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
powers
# determine soft threshold
soft <- pickSoftThreshold(brexp_rand, powerVector = powers, verbose = 5)
str(soft)

# plotting to determine where the values are scale free
# in this case it is 6
plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
	xlab = "Soft Threshold Power", ylab = "Scale Free Topology", type = "n",
	main = paste("scale independence"))
text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
	labels = powers , cex = 0.9, col = "red")
abline(h = 0.90, col = "red")

adj_mtx <- adjacency(brexp_rand, power = 6)
head(adj_mtx)[,1:10]
system.time(TOM <- TOMsimilarity(adj_mtx))
dissTOM <- 1 - TOM

geneTree <- flashClust(as.dist(dissTOM), method = "average")

plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
	labels = FALSE, hang = 0.04)

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
	                         deepSplit = 2, pamRespectsDendro = FALSE, 
	                         minClusterSize = 30)

table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
	                dendroLabels = FALSE, hang = 0.03,
	                addGuide = TRUE, guideHang = 0.05,
	                main = "Gene dendrogram and module Colors")

MEList <- moduleEigengenes(brexp_rand, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss = 1 - cor(MEs)


METree <- flashClust(as.dist(MEDiss), method = "average")
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module Eigengenes")
#add threshold of 0.85 cutoff for merging modules
abline(h=0.15, col = "red")

merged_clusters <- mergeCloseModules(brexp_rand, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merged_clusters$colors
mergedMEs <- merged_clusters$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                     c("Dynamic Tree Cut", "merged dynamic"),
	                dendroLabels = FALSE, hang = 0.03,
	                addGuide = TRUE, guideHang = 0.05)

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs


# now for all the traits
system.time( bwnet <- blockwiseModules(brexp, maxBlockSize = 5000, power = 6, 
	                      TOMType = "signed", networkType = "signed",
	                      minModuleSize = , reassignThreshold = 0, 
	                      mergeCutHeight = 0.25, numericLabels = TRUE,
	                      saveTOMs = TRUE, saveTOMFileBase ="brass_network_TOM",
	                      nThreads = 3, verbose = 3))

dynamicColors <- labels2colors(bwnet$colors)
table(dynamicColors)

MEsO <- bwnet$MEs 

MEDiss = 1 - cor(MEs)
METree <- flashClust(as.dist(MEDiss), method = "average")
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module Eigengenes")
#add threshold of 0.85 cutoff for merging modules
abline(h=0.15, col = "red")


setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/") 
brassica_blups <- read.table("Brapa2012allmeans_reduced_ASPB_poster.csv", header=TRUE, sep = ",")
head(brassica_blups)
tail(brassica_blups)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_blups2 <- read.delim("2011_phenology_heights_blups.csv", 
	                               header = TRUE, row.names = 1, sep = ",")
head(br_blups2)

#for Brapa2012allmeans-1.csv
colnames(brassica_blups)[1] <- paste("RILs")
brassica_blups$RILs <- gsub("(\\d+)","RIL_\\1", brassica_blups$RILs) 
head(brassica_blups)
row.names(brassica_blups) <- brassica_blups$RILs
brassica_blups <- brassica_blups[,-1]

keeps <- rownames(brexp)
keeps
blup_keeps <- brassica_blups[match(keeps, row.names(brexp)), ]
blup_keeps
dim(brassica_blups)
dim(brexp)
dim(blup_keeps)


br_blups2_keep <- br_blups2[match(keeps, row.names(brexp)), ]
br_blups2_keep
dim(brassica_blups)
dim(brexp)
dim(br_blups2_keep)

nGenes <- ncol(brexp)
nSamples <- nrow(brexp)

#correlate triats and eigen values
moduleTraitCor <- cor(MEs, blup_keeps, use = "p")
?corPvalueStudent
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCor <- cor(MEs, br_blups2_keep, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#plot correlations
textMatrix <- paste(signif(moduleTraitCor,2), "\n(",
	                signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c())

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
png(file="expression_pheno_clusters_2011.png",
	width=2500,height=1500, res=100)
labeledHeatmap(Matrix = moduleTraitCor,
	           xLabels = names(br_blups2_keep),
	           yLabels = names(MEs),
	           ySymbols = names(MEs),
	           colorLabels = FALSE,
	           colors = blueWhiteRed(50),
	           textMatrix = textMatrix,
	           setStdMargins = FALSE,
	           cex.text = 0.5,
	           zlim = c(-1,1),
	           main = paste("Module-trait Relationships")
	           )
dev.off()
?labeledHeatmap

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
save.image(file = "clustering.RData", version = NULL,
 ascii = FALSE, safe = TRUE)

