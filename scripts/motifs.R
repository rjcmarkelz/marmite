##### motif search #####
# shade gene list
# promotors of these genes
# eQTL of these genes
# trans eqtl hotspot
# large cis eQTL genes in hotspot?
# TFs in that region?
library(reshape2)
library(ggplot2)
library(Biostrings)
library(MotifDb)
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)


load('~/git.repos/brassica_eqtl_v1.5/data/un_eqtl.RData')
head(cistrans_df)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_shade <- read.delim("br_shade_genes.csv", header = TRUE, sep = ",")

# combine pieces of information
scanone_imp_tot[c(1,2)]
scanone_imp_tot2 <- scanone_imp_tot
scanone_imp_tot2$marker <- rownames(scanone_imp_tot2)
dim(scanone_imp_tot2)
scanone_imp_tot2 <- scanone_imp_tot2[c(35042,1:35041)]

# a few genes have multiple arabidopsis hits
shade_qtl <- scanone_imp_tot2[c(1,2,3, br_shade$V1)]
str(shade_qtl)
dim(shade_qtl)
head(shade_qtl)

# remove duplicated data
shade_qtl <- shade_qtl[-c(4:6)]

# melt for plotting
head(shade_qtl)[1:10]
shade_melt <- melt(shade_qtl , id = c("marker", "chr", "pos"))
head(shade_melt)
colnames(shade_qtl)

# plot all chromosomes
peak <- max(shade_melt$value)
shade_plot <- ggplot(shade_melt)
shade_plot <- shade_plot +  theme_bw() + geom_line(aes(x = pos, y = value, color = variable), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(chr ~ .) +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
shade_plot


# take a look at a specific chromosome
shade_melt_sub <- subset(shade_melt, chr == "A03")
head(shade_melt_sub)

# subset plot of chromosomes for closer look
peak <- max(shade_melt_sub$value)
shade_sub_plot <- ggplot(shade_melt_sub)
shade_sub_plot <- shade_sub_plot +  theme_bw() + geom_line(aes(x = pos, y = value, color = variable), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        theme(legend.position="none") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
shade_sub_plot

# promotors
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")
promoters <- readDNAStringSet("Brapa_1000bp_upstream_3.fa")
# promoters <- readFASTA("Brapa_1000bp_upstream_3.fa")

# subset based on shade genes
# subset based on treatment effects
# subset based on metabolism
# trans hotspots

# for some reason does not accept last entry
# check specifically last entry
br_shade
sh_prom <- promoters[head(br_shade$V1, 188)]
sh_prom <- unique(sh_prom)
sh_prom
# 170
head(MotifDb)
str(MotifDb)
?motifEnrichment
?query
data(MotifDb.Dmel.PFM)
MotifDb.Dmel.PFM
arab_motif <- grep('Athaliana', values (MotifDb)$organism)
str(arab_motif)
arab_motif_2 <- MotifDb[arab_motif]
head(str(arab_motif_2))
test <- matrix(as.list(arab_motif_2[1]))
test <- as.list(arab_motif_2@listData)
str(test)
?toPWM
pwm <- toPWM(test)
str(pwm)
pwm

# Needs list of PWMs using toPWM
?motifEnrichment
res <- motifEnrichment(sh_prom_l, pwm)
res

report <- groupReport(res)
report
plot(report[1:10], fontsize=7, id.fontsize=5)
warnings()

load('~/git.repos/brassica_genome_db/raw_data/brassica_gene_db.RData')
qtl_range <- GRanges(seqnames = "A03", ranges = IRanges(start=9500000, end=10000000))
qtl_range 
?transcriptsByOverlaps
str(brassica_db)
qtl_range_candidates <- transcriptsByOverlaps(brassica_db, qtl_range)
flr_range <- as.data.frame(qtl_range_candidates$tx_name)
qtl_range_candidates
flr_range
br_flr <- read.delim("br_flowering_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_flr)
colnames(br_flr)
str(br_flr[9])
colnames(flr_range) <- "y"
colnames(br_flr[9]) <- "x"

br_flr[9] <- as.character(br_flr[9])
flr_range <- as.character(flr_range)
colnames(flr_range) <- "y"
flr_range

merge(br_flr, flr_range, by.x = "Gene.name.in.B..rapa", by.y = "y")
