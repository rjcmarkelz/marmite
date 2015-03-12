# Quick script to take a look at some shallow sequencing data for DK and JC
setwd("/Users/Cody_2/git.repos/brassica_field_2014/raw_data")

counts <- read.table("Gene_counts_corwin.csv", sep = ",", header = TRUE)
head(counts)
str(counts)
dim(counts)

counts_names <- names(counts)
counts_names

head(counts)[97:100]

lane1 <- counts[1:97]
head(lane1)

row.names(lane1) <- lane1$geneID
lane1 <- lane1[,-1]
?rowSums
lane1$total <- rowSums(lane1, na.rm = TRUE)
lane1$means <- rowMeans(lane1, na.rm = TRUE)
lane1$means
hist(lane1$total)
hist(lane1$means)

row.names(lane1)
lane1$dup_names <- as.character(row.names(lane1))
lane1$dup_names
lane1$dup_names <- sub("(.+)(\\.\\d)","\\1", lane1$dup_names)
head(lane1$dup_names)
str(lane1$dup_names)

dim(lane1)
# [1] 200  99


lane1 <- lane1[!duplicated(lane1$dup_names),]
dim(lane1)
# [1] 90 99

lane1 <- lane1[lane1$gene_unique,]
head(lane1)

hist(lane1$total)
hist(lane1$means)
lane1


sum(is.na(lane1))
dim(lane1)

lane1$row_NA <- apply(lane1, 1, function(x) length(which(is.na(x))))
colnames(lane1)
plot(row_NA ~ means, data = lane1)
head(lane1)

gene_func <- read.table("C_N_metabolism_arabidopsis_corwin.csv", sep = ",", header = TRUE)
head(gene_func)
dim(gene_func)
# [1] 213   2
gene_func <- gene_func[!duplicated(gene_func$GeneID),]
dim(gene_func)
# [1] 96  2


merged <- merge(lane1, gene_func, by.x = "dup_names", by.y = "GeneID", all.x = TRUE)
dim(merged)
head(merged)

setwd("/Users/Cody_2/git.repos/brassica_field_2014/raw_data")
write.table(merged, "C_N_metabolism_subset_shallow.csv", sep = ",", row.names = FALSE)

