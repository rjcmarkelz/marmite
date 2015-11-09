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


#use hotspot library to extract qtl for each gene
library(qtlhot)
set.seed(12345)

permtest <- scanone(brass_total, method = "imp", n.perm = 1000)
permtest

alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs <- summary(permtest, alphas)

# get 3% cutoff
lod.thrs
lod.thr <- lod.thrs[1]
# decided to use more stringent cutoff, will use permuations once stopped running

#reduce object size by getting removing NS peaks
ls()
?highlod
high1 <- highlod(scanone_imp_tot, lod.thr = 5, drop.lod = 0)
head(high1)
dim(scanone_imp_tot)
head(scanone_imp_tot)[1:10]

#exploring high1 object to figure out how to subset it for what I want
length(high1)
class(high1)
attributes(high1)
max(high1, lod.thr = lod.thr)
head(high1)
head(high1$chr.pos)
head(high1$highlod, 20)
head(high1$names)


#take a look at hotspots
hots1 <- hotsize(high1, lod.thr = lod.thr)
hots1

#infile genomic coordinates of genes
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")

transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)
# [1] 43463     6 
# includes scaffolds, need to figure this out
transcripts$tx_name <- as.character(transcripts$tx_name)
transcripts$tx_chrom <- as.character(transcripts$tx_chrom)

# it is just easier to use merge to get the type of dataframe that is useful for plotting
# and other applications for now
markers <- as.data.frame(rownames(scanone_imp_tot))
head(markers)
markers[1]
names(markers) <- "marker_name"
markers$index <- rownames(markers)

lods <- as.data.frame(high1$highlod)
head(lods)
dim(lods)

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

#load plyr library to do a quick string split
library(plyr)
chr_pos <- ldply(strsplit(as.character(cistrans_df$marker_name), split = "x"))
head(chr_pos)
cistrans_df$qtl_chrom <- chr_pos$V1
cistrans_df$qtl_pos <- chr_pos$V2

head(cistrans_df)
str(cistrans_df)

# ignoring the scaffolds for a moment
cistrans_df$cis_trans <- ifelse(cistrans_df$tx_chrom == cistrans_df$qtl_chrom, paste("cis"), paste("trans"))
#take a quick look
head(cistrans_df)
cistrans_df[10000:10100,]
tail(cistrans_df)

# with LOD threshold of 5 and 0 lod support interval
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(cistrans_df, "cis_trans_scanone_un.csv", sep = ",")
save.image(file = "un_eqtl.RData", version = NULL, ascii = FALSE, safe = TRUE)


###########
# 2014_4_11
###########
dim(cistrans_df)

# [1] 40034    13
unique(cistrans_df$tx_name)
# 35039

cis_df <- subset(cistrans_df, cis_trans == "cis")
head(cis_df)
dim(cis_df)
# [1] 26964    13

trans_df <- subset(cistrans_df, cis_trans == "trans")
head(trans_df)
dim(trans_df)
# [1] 13070    13

str(cis_df)
cis_df$qtl_pos <- as.numeric(cis_df$qtl_pos)
cis_df$distMbp <- abs(cis_df$tx_start - cis_df$qtl_pos)/1000000
head(cis_df, 25)
plot(cis_df$distMbp)
plot(cis_df$lod)
#################
# Use mid-gene coordinate instead of startsite. Some are plus and some are minus strand, so loss of information
# if just using the start site from the "top" of the chromosome.
#################



# in the future this might be useful, but would have to change some the S3 object attributes and 
# methods to deal with scaffolds and character names for the chromosomes
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

##### 
#explore shade genes
#####
# use un_eqtl.R



scanone_imp_tot[c(1,2)]
scanone_imp_tot2 <- scanone_imp_tot
scanone_imp_tot2$marker <- rownames(scanone_imp_tot2)
dim(scanone_imp_tot2)
scanone_imp_tot2 <- scanone_imp_tot2[c(35042,1:35041)]


# a few genes have multiple arabidopsis hits
shade_qtl <- scanone_imp_tot2[c(1,2,3, br_shade$V1)]
str(shade_qtl)
dim(shade_qtl)
plot(shade_qtl)

shade_qtl <- shade_qtl[-c(4:6)]

library(reshape2)
library(ggplot2)

head(shade_qtl)[1:10]
shade_melt <- melt(shade_qtl , id = c("marker", "chr", "pos"))
head(shade_melt)
colnames(shade_qtl)

shade_melt_sub <- subset(shade_melt, chr == "A03")
head(shade_melt_sub)


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
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
shade_plot

peak <- max(shade_melt_sub$value)
shade_sub_plot <- ggplot(shade_melt_sub)
shade_sub_plot <- shade_sub_plot +  theme_bw() + geom_line(aes(x = pos, y = value, color = variable), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # facet_grid(chr ~ .) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
shade_sub_plot

#check out flowering genes
ase_cent[ase_cent$gene_name =="Bra009055",] #FLC chr A10
ase_cent[ase_cent$gene_name =="Bra028599",] #FLC
ase_cent[ase_cent$gene_name =="Bra006051",] #FLC

setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
gene_coord <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(gene_coord)
head(br_shade)

shade_coord <- merge(br_shade, gene_coord, by.x = "V1", by.y = "tx_name", all.x = TRUE)
dim(shade_coord)
head(shade_coord)
shade_coord

trans_shade <- subset(shade_coord, tx_chrom != "A03")
head(trans_shade)
dim(trans_shade)
trans_shade

# a few genes have multiple arabidopsis hits
shade_qtl <- scanone_imp_tot2[c(1,2,3, trans_shade$V1)]
str(shade_qtl)
dim(shade_qtl)
plot(shade_qtl)
head(shade_qtl)

head(shade_qtl)[1:10]
shade_melt <- melt(shade_qtl , id = c("marker", "chr", "pos"))
head(shade_melt)
colnames(shade_qtl)


shade_melt_sub2 <- subset(shade_melt, chr == "A03")
head(shade_melt_sub2)
shade_melt_sub2

shade_melt_sub2$physical <- as.numeric(sub("(A03)(x)(\\d+)", "\\3", shade_melt_sub2$marker))/1000000
head(shade_melt_sub2)
str(shade_melt_sub2)

peak <- max(shade_melt_sub2$value)
trans_shade_plot <- ggplot(shade_melt_sub2)
trans_shade_plot <- trans_shade_plot +  theme_bw() + geom_line(aes(x = physical, y = value, color = variable), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = physical, xend = physical), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(chr ~ .) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        theme(legend.position="none") +
                        xlab("Genomic Position (Mbp)") +
                        ylab("LOD Score") 
trans_shade_plot




#metabolism genes
# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

br_met <- read.delim("br_met_genes.csv", header = TRUE, sep = ",")
head(br_met)
head(scanone_imp_tot)[1:10]
br_shade

scanone_imp_tot[c(1,2)]
scanone_imp_tot2 <- scanone_imp_tot
scanone_imp_tot2$marker <- rownames(scanone_imp_tot2)
dim(scanone_imp_tot2)
scanone_imp_tot2 <- scanone_imp_tot2[c(35042,1:35041)]


# a few genes have multiple arabidopsis hits
met_qtl <- scanone_imp_tot2[c(1,2,3, br_met$V1)]
str(met_qtl)
dim(met_qtl)
plot(met_qtl)

met_qtl <- scanone_imp_tot2[c(1,2,3, br_met$V1)]
str(met_qtl)
dim(met_qtl)
plot(met_qtl)

met_qtl <- met_qtl[-c(4:7)]

library(reshape2)
library(ggplot2)
melt?
?melt
head(met_qtl)[1:10]
dim(met_qtl)

# do partial melts of the dataframes
# follows where the cis-eqtl are
met_melt <- melt(met_qtl[c(1:300)] , id = c("marker", "chr", "pos"))
head(met_melt)
colnames(met_qtl)
dim(met_melt)

peak <- max(met_melt$value)
met_plot <- ggplot(met_melt)
met_plot <- met_plot +  theme_bw() + geom_point(aes(x = pos, y = value), size = 2) +
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
?ggsave
met_plot
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/output")
ggsave("met_plot.pdf", met_plot, width = 20, height = 30)
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

