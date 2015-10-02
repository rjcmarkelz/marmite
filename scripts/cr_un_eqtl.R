
###########
###########
#RQTL MAPPING
###########
###########
library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

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

brassica_genes <- est.rf(brassica_genes)
plot.rf(brassica_genes) 

#about a minute
brassica_genes <- calc.errorlod(brassica_genes, error.prob=0.001)

# test
system.time(scanone.imp.1 <- scanone(brassica_genes, pheno.col = 1,
 method = "imp", use="all.obs"))

#rerun whole dataset after chr reassignment
system.time(scanone.imp.1 <- scanone(brassica_genes, pheno.col = 1:35039,
 method = "imp", use="all.obs"))

#    user  system elapsed 
# 845.002  11.957 854.041 

head(scanone.imp.1)


# save output in data directory
save.image(file = "cr_un_eqtl.RData", version = NULL,
 ascii = FALSE, safe = TRUE)


# eQTL hotspot analysis

library(qtlhot)
set.seed(12345)

permtest <- scanone(brassica_genes, method = "imp", n.perm = 1000)
permtest
alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs <- summary(permtest, alphas)

# get 5% cutoff
lod.thrs
lod.thr <- lod.thrs[5]

#reduce object size by getting removing NS peaks
high1 <- highlod(scanone.imp.1, lod.thr = lod.thr, drop.lod = 1.5)
max(high1, lod.thr = lod.thrs)
high1
hots1 <- hotsize(high1, lod.thr = lod.thr)
high1
summary(hots1)
# LOD threshold: 3.125673 
#              chr  pos max.N
# A01x9141005  A01 34.2    15
# A02x20560232 A02 88.0    15
# A03x4610765  A03 27.3    53
# A04x13790497 A04 42.1    10
# A05x25053034 A05 87.8     4
# A06x16894473 A06 42.7     6
# A07x11085819 A07 24.4    39
# A08x17544268 A08 43.1   231
# A09x14278053 A09 68.0   244
# A10x9187380  A10 40.9    30

set.seed(12345)
plot(hots1, cex.lab = 1.5, cex.axis = 1.5)

str(hots1)

#remove this phenotype
find.pheno(brassica_genes, "id")
brassica_genes$pheno <- brassica_genes$pheno[1:35039]

?hotperm
hotperm1 <- hotperm(cross = brassica_genes,
	         n.quant = 300, n.perm = 100, lod.thrs = lod.thrs,
	         alpha.levels = alphas, drop.lod = 1.5, verbose = FALSE)
summary(hotperm1)

quant1 <- quantile(hotperm1, 0.05, lod.thr = lod.thrs)
quant1
plot(high1, quant.level = quant1, sliding = TRUE)
hotsq1 <- hotsize(high1, lod = lod.thr, window = 5)
plot(hotsq1)
hotsq1


save.image(file = "cr_un_eqtl.RData", version = NULL,
 ascii = FALSE, safe = TRUE)


# cis trans analysis
AO1_eqtl <- highlod(scanone.imp.1, chr = "A01" , lod.thr = lod.thr, drop.lod = 1.5)
str(AO1_eqtl)
highlodA01 <- as.data.frame(AO1_eqtl$highlod)
head(highlodA01)
dim(highlodA01)

summary(scanone.imp.1)
system.time(brassica_genes_sig <- as.data.frame(summary(scanone.imp.1, 
	                          threshold = 3.15, format = "allpeaks" )))

system.time(scanone.imp.4 <- scanone(brassica_genes, pheno.col = 1:100,
 method = "imp", use="all.obs"))
system.time(brassica_genes_sig <- as.data.frame(summary(scanone.imp.4, 
	                          threshold = 3.15, format = "allpeaks" )))
head(brassica_genes_sig)


library(eqtl)
#define the peaks from the scanone object
brassica_peaks <- define.peak(scanone.imp.1, lodcolumn= 'all', th = 3.5, si = 1.5)
class(brassica_peaks)
attributes(brassica_peaks)
summary(brassica_peaks)
str(brassica_peaks)

define.peak(scanone.imp.1, th = 3.5, lodcolumn= 1:35000, graph = TRUE, chr=c("A01"));

# make a copy to modify
str(scanone.imp.1)
scanone.imp.2 <- scanone.imp.1
str(scanone.imp.2)

scanone.imp.2$chr <- as.numeric(sub("(A)(.+)","\\2", scanone.imp.2$chr))
scanone.imp.2$chr

head(scanone.imp.1)


head(ATH.coord)
str(ATH.coord)
#infile genome data
#used the BEDtoGFF.pl script to convert the BED file to GTF
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
# update with new genome mapping
br_bed <- read.table("Brassica_rapa_v1.5_final.bed", header = FALSE)
head(br_bed)
str(br_bed)

#transcripts data frame
transcripts <- br_bed[,c("V4","V4","V1","V6","V2","V3")]
head(transcripts, 10)
tail(transcripts, 20)

#make unique id for each gene based on row number for now 1-44239
transcripts$V4 <- rownames(transcripts)
head(transcripts)

#replace names of new 
colnames(transcripts) <- c("tx_id", "tx_name", "tx_chrom", "tx_strand",
	                         "tx_start", "tx_end")
head(transcripts)
str(transcripts)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
# write table for eQTL cis/trans
write.table(transcripts, "transcripts_eqtl_start_stop_eqtl.csv", sep = ",", col.names = TRUE, row.names = FALSE)

#format for eQTL
transcripts <- transcripts[c(2,3,5,6)]
head(transcripts)
tail(transcripts)

transcripts$tx_chrom  <- sub("(A0)(\\d)","\\2", transcripts$tx_chrom)
transcripts$tx_chrom  <- sub("(A10)","10", transcripts$tx_chrom)
transcripts <- transcripts[!grepl("Scaff", transcripts$tx_chrom),]
colnames(transcripts) <- c("etrait.name", "chr", "start", "stop")
str(transcripts)
transcripts$etrait.name <- as.character(transcripts$etrait.name)
transcripts$chr <- as.numeric(transcripts$chr)

######
#
######
brass_local <- localize.qtl(cross = scanone.imp.1, peak = brassica_peaks, data.gmap = transcripts);


# save output in data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
save.image(file = "cr_un_eqtl.RData", version = NULL,
 ascii = FALSE, safe = TRUE)









