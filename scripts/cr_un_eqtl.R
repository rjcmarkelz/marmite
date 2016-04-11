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
str(scanone.imp.1)

brassica_genes <- calc.genoprob(brassica_genes)

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
?highlod
high1 <- highlod(scanone.imp.1, lod.thr = lod.thr, drop.lod = 1.5)
max(high1, lod.thr = lod.thrs)
high1
hots1 <- hotsize(high1, lod.thr = lod.thr)
high1

high_9 <- highlod(scanone.imp.1, chr = "A09",lod.thr = lod.thr, drop.lod = 1.5 )
high_9


high_9 <- pull.highlod(high1, chr = A09, pos =68.0)
high_9
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
?highlod
str(scanone.imp.1)
scanone.imp.1[2:3,]
head(rownames(scanone.imp.1))
tail(rownames(scanone.imp.1))
summary(scanone.imp.1)

scanone_rows <- as.data.frame(rownames(scanone.imp.1))
head(scanone_rows)
rownames(scanone_rows)
scanone_rows$index <- as.numeric(rownames(scanone_rows))
head(scanone_rows)
#rownames of scanone objects are marker names

A01_eqtl <- highlod(scanone.imp.1, chr = "A01" , lod.thr = lod.thr, drop.lod = 0.1, restrict.lod = TRUE)
str(A01_eqtl)
max(A01_eqtl$highlod$row)

# pull eqtl
highlod_out <- as.data.frame(A01_eqtl$highlod)
max(highlod_out$lod)
#[1] 6.632117
head(highlod_out)
dim(highlod_out)

# pull marker indexs 
head(scanone_rows)

# pull gene indexs
gene_indx <- as.data.frame(A01_eqtl$names)
head(gene_indx)
gene_indx$gene_index <- as.numeric(rownames(gene_indx))

highlod_marker <- merge(highlod_out, scanone_rows, by.x = "row", by.y = "index", all.x = TRUE)
dim(highlod_marker)
head(highlod_marker)

highlod_marker_gene <- merge(highlod_marker, gene_indx, by.x = "phenos", by.y = "gene_index", all.x = TRUE)
dim(highlod_marker)
head(highlod_marker_gene)
colnames(highlod_marker_gene)[4:5] <- c("marker_names", "gene_names")

# might as well include the marker genetic map data
# pull gene indexs
marker_indx <- as.data.frame(A01_eqtl$chr.pos)
head(marker_indx)
marker_indx$marker_index <- as.character(rownames(marker_indx))


gxe_genes_markers <- merge(highlod_marker_gene, marker_indx, by.x = "marker_names", by.y = "marker_index", all.x = TRUE)
head(gxe_genes_markers)

# add physical position to compare with other plots
gxe_genes_markers$Mbp <- as.numeric(sub("(A\\w+)(x)(\\d+)", "\\3", gxe_genes_markers$marker_names))
head(gxe_genes_markers)
gxe_genes_markers$Mbp <- gxe_genes_markers$Mbp/1000000


setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
# write table for eQTL cis/trans
write.table(gxe_genes_markers, "gxe_genes_markers.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# save output in data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
save.image(file = "brassica_genes_sig.RData", version = NULL,
 ascii = FALSE, safe = TRUE)

# now for some plotting
dim(gxe_genes_markers)

# load data


library(ggplot2)


gxe_plot <- ggplot(gxe_genes_markers)
gxe_plot <- gxe_plot +  theme_bw() + geom_point(aes(x = Mbp, y = lod,  alpha = 0.1), size = 2) +
                        facet_grid(chr ~ .) + 
                        ggtitle("G by E eQTL distribution") +
                        xlab("Genomic Position (Mbp)") +
                        ylab("LOD Score") +
                        scale_y_continuous(limits= c(3, 7)) +
                        theme(axis.title=element_text(face="bold",size="20"),
                              axis.text=element_text(face="bold", size="12"),
                              legend.position = "none")  
gxe_plot
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/output")
ggsave(gxe_plot, file = "gxe_hotspots.pdf", width = 10, height = 20 )

# infile subgenomes
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

trans_cent <- read.table("trans_cent.csv", sep = ",", header = TRUE)
head(trans_cent)

gxe_trans_cent <- merge(gxe_genes_markers, trans_cent, by.x = "gene_names", by.y = "geneID", all.x = TRUE)
head(gxe_trans_cent)
dim(gxe_trans_cent)


######start here########
######start here########
######start here########
######start here########
gxe_trans_cent_plot <- ggplot(data = trans_cent[!(trans_cent$sub_genome == "NA"),])
gxe_trans_cent_plot <- gxe_trans_cent_plot +  theme_bw() + 
                    geom_rect(data = trans_cent[!(trans_cent$sub_genome == "NA"),],
                       aes(xmin = tx_start, xmax = tx_end, ymin = -1, ymax = 3, color = sub_genome), size = 4) +
                    geom_point(data = gxe_genes_markers,
                    	aes(x = Mbp, y = lod), size = 3, alpha = 0.3, color = "black") +
                    facet_grid(Chr ~ . ) +
                    # geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("trans-eQTL number") 

gxe_trans_cent_plot
######start here########




setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/output")
ggsave(gxe_trans_cent_plot, file = "gxe_hotspots.pdf", width = 10, height = 20 )


# Do the GxE hotspots have an overlap with known shade genes?
# load cr_un_eqtl
# eqtl dataset
brassica_genes

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_shade <- read.delim("br_shade_genes.csv", header = TRUE, sep = ",")
head(br_shade)
head(scanone.imp.1)[1:10]
br_shade
dim(br_shade)

scanone.imp.1[c(1,2)]
scanone.imp.2 <- scanone.imp.1
scanone.imp.2$marker <- rownames(scanone.imp.2)
dim(scanone.imp.2)
scanone.imp.2 <- scanone.imp.2[c(35042,1:35041)]


# a few genes have multiple arabidopsis hits
shade_qtl <- scanone.imp.2[c(1,2,3, br_shade$V1)]
str(shade_qtl)
dim(shade_qtl)

head(shade_qtl)

shade_qtl <- shade_qtl[-c(4:6)]

library(reshape2)
library(ggplot2)


shade_melt <- melt(shade_qtl , id = c("marker", "chr", "pos"))
head(shade_melt)


peak <- max(shade_melt$value)
shade_plot <- ggplot(shade_melt)
shade_plot <- shade_plot +  theme_bw() + geom_line(aes(x = pos, y = value, color = variable), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(chr ~ .) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome (cM)") +
                        ylab("LOD Score") 
shade_plot



















#####
# did not end up using this
#####
# summary(scanone.imp.1)
# system.time(brassica_genes_sig <- as.data.frame(summary(scanone.imp.1, 
# 	                          threshold = 3.15, format = "allpeaks" )))

# system.time(scanone.imp.4 <- scanone(brassica_genes, pheno.col = 1:100,
#  method = "imp", use="all.obs"))
# system.time(brassica_genes_sig <- as.data.frame(summary(scanone.imp.4, 
# 	                          threshold = 3.15, format = "allpeaks" )))
# head(brassica_genes_sig)

# library(eqtl)
# #define the peaks from the scanone object
# brassica_peaks <- define.peak(scanone.imp.1, lodcolumn= 'all', th = 3.5, si = 1.5)
# class(brassica_peaks)
# attributes(brassica_peaks)
# summary(brassica_peaks)
# str(brassica_peaks)

# define.peak(scanone.imp.1, th = 3.5, lodcolumn= 1:35000, graph = TRUE, chr=c("A01"));

# # make a copy to modify
# str(scanone.imp.1)
# scanone.imp.2 <- scanone.imp.1
# str(scanone.imp.2)

# scanone.imp.2$chr <- as.numeric(sub("(A)(.+)","\\2", scanone.imp.2$chr))
# scanone.imp.2$chr

# head(scanone.imp.1)


# head(ATH.coord)
# str(ATH.coord)

# #infile genome data
# #used the BEDtoGFF.pl script to convert the BED file to GTF
# setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
# # update with new genome mapping
# br_bed <- read.table("Brassica_rapa_v1.5_final.bed", header = FALSE)
# head(br_bed)
# str(br_bed)

# #transcripts data frame
# transcripts <- br_bed[,c("V4","V4","V1","V6","V2","V3")]
# head(transcripts, 10)
# tail(transcripts, 20)

# #make unique id for each gene based on row number for now 1-44239
# transcripts$V4 <- rownames(transcripts)
# head(transcripts)

# #replace names of new 
# colnames(transcripts) <- c("tx_id", "tx_name", "tx_chrom", "tx_strand",
# 	                         "tx_start", "tx_end")
# head(transcripts)
# str(transcripts)

# setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
# # write table for eQTL cis/trans
# write.table(transcripts, "transcripts_eqtl_start_stop_eqtl.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# #format for eQTL
# transcripts <- transcripts[c(2,3,5,6)]
# head(transcripts)
# tail(transcripts)

# transcripts$tx_chrom  <- sub("(A0)(\\d)","\\2", transcripts$tx_chrom)
# transcripts$tx_chrom  <- sub("(A10)","10", transcripts$tx_chrom)
# transcripts <- transcripts[!grepl("Scaff", transcripts$tx_chrom),]
# colnames(transcripts) <- c("etrait.name", "chr", "start", "stop")
# str(transcripts)
# transcripts$etrait.name <- as.character(transcripts$etrait.name)
# transcripts$chr <- as.numeric(transcripts$chr)

# ######
# #
# ######
# ?localize.qtl
# brass_local <- localize.qtl(cross = brassica_genes, peak = brassica_peaks, data.gmap = transcripts)
# brass_local

# # save output in data directory
# setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
# save.image(file = "cr_un_eqtl.RData", version = NULL,
#  ascii = FALSE, safe = TRUE)









