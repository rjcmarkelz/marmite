###############
###############
# total eqtl no treatment
###############
###############
# br_blues_total_RQTL.csv was created as part of the data clean up in cr_un_eqtl.R
library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

brass_total <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="br_blues_total_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brass_total)

class(brass_total)[1] <- "riself"
brass_total <- jittermap(brass_total)
brass_total

brass_total <- est.rf(brass_total)
plot.rf(brass_total) 

#about a minute
brass_total <- calc.errorlod(brass_total, error.prob=0.001)

system.time(scanone_imp_tot <- scanone(brass_total, pheno.col = 1:35039, 
	         method = "imp", use="all.obs"))

system.time(scanone_imp_tot <- scanone(brass_total, pheno.col = 1:35039, 
	         method = "imp", use="all.obs"))

save.image(file = "un_eqtl.RData", version = NULL,
 ascii = FALSE, safe = TRUE)

library(qtlhot)
set.seed(12345)

permtest <- scanone(brass_total, method = "imp", n.perm = 1000)
permtest

alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs <- summary(permtest, alphas)

# get 3% cutoff
lod.thrs
lod.thr <- lod.thrs[1]


#reduce object size by getting removing NS peaks
?highlod
high1 <- highlod(scanone_imp_tot, lod.thr = lod.thr, drop.lod = 1.5)
dim(scanone_imp_tot)
head(scanone_imp_tot)[1:10]

length(high1)
class(high1)
attributes(high1)
max(high1, lod.thr = lod.thrs)
head(high1)
hots1 <- hotsize(high1, lod.thr = lod.thr)
hots1

head(high1$chr.pos)

head(high1$highlod, 20)
head(high1$names)

high1$chr.pos

setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")

transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)
# [1] 43463     6 
# includes scaffolds
transcripts$tx_name <- as.character(transcripts$tx_name)
transcripts$tx_chrom <- as.character(transcripts$tx_chrom)





cistrans_calc <- function(scanone_out, highlod_out, coord ){
	so_rownames <- rownames(scanone)
	hl_genenames <- highlod$names
	# match gene names with genomic position
	gene_pos <- coord[coord$tx_name %in% hl_genenames,]


}
head(transcripts)
cistrans(high1, transcripts)
hl_genenames

hl_genenames <- high1$names
str(hl_genenames)
gene_pos <- transcripts[transcripts$tx_name %in% hl_genenames,]
head(gene_pos)
dim(gene_pos)
rownames(scanone_imp_tot)[10]
high1$names[2]



lods <- as.data.frame(high1$highlod)
head(lods)

markers <- as.data.frame(rownames(scanone_imp_tot))
head(markers)
markers[1]
names(markers) <- "marker_name"
markers$index <- rownames(markers)

gene_names <- as.data.frame(high1$names)
head(gene_names)
names(gene_names) <- "tx_name"
gene_names$old <- rownames(gene_names)

gene_lods <- merge(gene_names, lods, by.x = "old", by.y = "phenos", all.y = TRUE)
head(gene_lods)
dim(gene_lods)

gene_lod_mark <- merge(gene_lods, markers, by.x = "row", by.y = "index", all.x = TRUE)
head(gene_lod_mark)
dim(gene_lod_mark)

cistrans_df <- merge(gene_lod_mark, transcripts, by.x = "tx_name", by.y = "tx_name", all.x = TRUE)
head(cistrans_df)
dim(cistrans_df)
tail(cistrans_df)

library(plyr)
chr_pos <- ldply(strsplit(as.character(cistrans_df$marker_name), split = "x"))
head(chr_pos)
cistrans_df$qtl_chrom <- chr_pos$V1
cistrans_df$qtl_pos <- chr_pos$V2

head(cistrans_df)
str(cistrans_df)

# ignoring the scaffolds for a moment
cistrans_df$cis_trans <- ifelse(cistrans_df$tx_chrom == cistrans_df$qtl_chrom, paste("cis"), paste("trans"))
cistrans_df[1500:1600,]
tail(cistrans_df)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(cistrans_df, "cis_trans_scanone_un.csv", sep = ",")












#################
# eQTL package
#################

library(eqtl)
?define.peak
brassica_peaks <- define.peak(scanone_imp_tot, th = 4.0, si = 2.0,
                               lodcolumn= "all", chr=c("A03"))
brassica_peaks_total <- define.peak(scanone_imp_tot, th = 4.0, si = 1.0,
                               lodcolumn= "all")
dim(brassica_peaks_total)
head(brassica_peaks_total)

class(brassica_peaks)
length(brassica_peaks)

peaks_red <- brassica_peaks[!is.na(sapply(brassica_peaks,`[`,1))]
attributes(peaks_red)$class <- c("peak", "list")
attributes(peaks_red)$features <- c("lod", "mname.peak", "peak.cM", 
        "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
attributes(peaks_red)$scanone <- "scanone_imp_tot"
attributes(peaks_red)$lod.th <- 4.0
attributes(peaks_red)$si <- 1.0
attributes(peaks_red)$window <- 20
attributes(peaks_red)
length(peaks_red)

# 3253 eQTL meeting sig threshold on chromosome 3
head(peaks_red)
names(peaks_red)
tail(peaks_red)




tester <- as.data.frame(peaks_red$Bra000007)
tester2 <- rbind(tester, tester)
tester2

tester2$gene <- paste("Bra000007")
tester2




# attributes(brassica_peaks)$scanone


# tail(brassica_peaks, 20)
# dim(brassica_peaks)


# # get significant peaks
# brassica_peaks[2][1]
# is.na(brassica_peaks[2][1])


# test <- brassica_peaks[1:10]
# class(test)
# head(test)
# ?pull.pheno
# testpheno <- colnames(pull.pheno(brass_total, pheno.col = 1:10))
# testpheno

# test2 <- as.data.frame(is.na(sapply(test,`[`,1)))
# str(test2)
# test2

# test
# test3 <- test[!is.na(sapply(test,`[`,1))]
# str(test3)

# # QUICK AND DIRTY, build back up peak object
# attributes(test3)$class <- c("peak", "list")
# attributes(test3)$features <- c("lod", "mname.peak", "peak.cM", 
#         "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
# attributes(test3)$scanone <- "scanone_imp_tot"
# attributes(test3)$lod.th <- 4.0
# attributes(test3)$si <- 1.0
# attributes(test3)$window <- 20
# attributes(test3)
# # got it

# peaks_red <- brassica_peaks[!is.na(sapply(brassica_peaks,`[`,1))]
# attributes(peaks_red)$class <- c("peak", "list")
# attributes(peaks_red)$features <- c("lod", "mname.peak", "peak.cM", 
#         "mname.inf", "inf.cM", "mname.sup", "sup.cM", "si.quality")
# attributes(peaks_red)$scanone <- "scanone_imp_tot"
# attributes(peaks_red)$lod.th <- 4.0
# attributes(peaks_red)$si <- 1.0
# attributes(peaks_red)$window <- 20
# attributes(peaks_red)
# length(peaks_red)
# head(peaks_red)
# names(peaks_red)


# peaks_A01 <- brassica_peaks_total[is.na(sapply(brassica_peaks_total,`[`,1))]
# peaks_red_A01 <- brassica_peaks_total[!is.na(sapply(brassica_peaks_total,`[`,1))]

# length(peaks_red_total)

# summary(scanone_imp_tot)

