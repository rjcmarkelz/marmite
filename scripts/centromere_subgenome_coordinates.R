library(ggplot2)
install.packages("lsr")
library(lsr)
#infile genomic coordinates of genes
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")


transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)

centromeres <- read.table("centromere_subgenome_coordinates.csv", sep = ",", header = TRUE)
head(centromeres)
str(centromeres)
dim(centromeres)

trans_cent <- merge(centromeres, transcripts, by.x = "geneID", by.y = "tx_name", all.x = TRUE)
head(trans_cent)
tail(trans_cent)
dim(trans_cent)

plot(trans_cent$Start, trans_cent$tx_start)
# one to one (nearly)

trans_cent$Mbp <- (trans_cent$tx_start/1000000)
trans_cent$y_hold <- as.numeric(paste(1))
head(trans_cent)
str(trans_cent)
trans_cent$sub_genome <- "NA"
trans_cent$placeholder <- 1
trans_cent$placeholder 

#figure 2 found it in another paper
# trans_cent$sub_genome[trans_cent$Block %in% c("U", "Tb", "Qb", "Xa", "D", "Ta", "S", "J", "Fa")] <- "LF"
# #notice Ma block is both LF and MF2 on chrm A06
# trans_cent$sub_genome[trans_cent$Block %in% c("Cb", "Ab", "Ba", "Ma", "Xb", "Qa", "Wa", "L", "K", "V")] <- "LF"
# trans_cent$sub_genome[trans_cent$Block %in% c("H", "G", "Fb", "Xc", "Ib", "E",)] <- "LF"
# trans_cent$sub_genome

# not finished yet
trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = placeholder, group = Block, color = Block, y = Mbp)) +
                        facet_grid(~ tx_chrom) 
trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = Mbp, y = placeholder, color = Block)) +
                        facet_grid(~ tx_chrom) 
trans_cent_plot

#genome coordinates
?read.table
subgenomes <- read.table("Cheng_et_al_2012_supple_subgenomes.csv", sep = ",", header = TRUE, 
	                   stringsAsFactors = FALSE   )
head(subgenomes)
str(subgenomes)
LF1 <- subgenomes$LF
LF1

MF1 <- subgenomes$MF1
MF1

MF2 <- subgenomes$MF2
MF2

head(trans_cent)
str(trans_cent)
trans_cent$geneID <- as.character(trans_cent$geneID)
trans_cent$sub_genome[trans_cent$geneID %in% LF1] <- "LF"
head(trans_cent)
trans_cent$sub_genome

trans_cent$sub_genome[trans_cent$geneID %in% MF1] <- "MF1"
trans_cent$sub_genome[trans_cent$geneID %in% MF2] <- "MF2"
trans_cent$sub_genome[trans_cent$sub_genome == "1"] <- "NA"



trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = Mbp, y = placeholder, color = sub_genome)) +
                        facet_grid(~ tx_chrom) 

trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_tile(data = trans_cent[(trans_cent$tx_chrom == "A05") &
                     !(trans_cent$sub_genome == "NA"),], 
	                 aes(x = Mbp, y = placeholder, color = sub_genome)) 

trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_tile(data = trans_cent[!(trans_cent$Block == "NULL"),],
                    aes(x = Mbp, y = placeholder, color = Block)) +
                    facet_grid(tx_chrom ~ .) 


trans_cent$Centromere[1:150]
trans_cent$centplot <- trans_cent$Centromere
trans_cent$centplot 

#change name so there is a yes or no
trans_cent$centplot <- sub("Centromere A01", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A02", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A03", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A04", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A05", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A06", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A07", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A08", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A09", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("centromere A10", "centromere", trans_cent$centplot)
trans_cent$centplot <- sub("Trace", "ancestral centromere", trans_cent$centplot)

trans_cent$tx_chrom
trans_cent$chromplot <- as.character(trans_cent$tx_chrom)

trans_cent$chromplot <- sub("(Scaffold)(.)+", "ignore", trans_cent$chromplot)
trans_cent$chromplot
dim(trans_cent)
str(trans_cent)
# [1] 36932    23
trans_cent <- trans_cent[!is.na(trans_cent$chromplot),]
dim(trans_cent)
# [1] 36633    23 
trans_cent$chromplot <- as.factor(trans_cent$chromplot)
str(trans_cent)

levels(trans_cent$Centromere)
levels(trans_cent$centplot)

# there were both sub W blocks in the this dataset, merge to one
trans_cent$Block <- sub("Wa", "W", trans_cent$Block)
trans_cent$Block <- sub("Wb", "W", trans_cent$Block)

######
#NEAR FINAL PLOT
######
head(trans_cent)
trans_cent$tx_start <- trans_cent$tx_start/1000000
trans_cent$tx_end <- trans_cent$tx_end/1000000
head(trans_cent)

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_rect(data = trans_cent[!(trans_cent$Block == "NULL"),],
                       aes(xmin = tx_start, xmax = tx_end, ymin = -1, ymax = 1, color = Block), size = 3) +
                    geom_point(aes(x = tx_start, y = 0, shape = centplot), size = 5, alpha = 0.1, color = "black") +
                    scale_shape_manual(values=c(23, 15)) +
                    facet_grid(chromplot ~ . ) + 
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("") + 
                    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major.y = element_blank(), 
                    	panel.grid.minor.y = element_blank())
trans_cent_plot

# the same plot as above, but just the LF1, MF1, MF2
head(trans_cent)
str(trans_cent)
subgenome_blocks <- ggplot(trans_cent)
subgenome_blocks <- subgenome_blocks +  theme_bw() + 
                    geom_rect(data = trans_cent[!(trans_cent$sub_genome == "NA"),],
                       aes(xmin = tx_start, xmax = tx_end, ymin = -1, ymax = 1, color = sub_genome), size = 2) +
                    geom_point(aes(x = tx_start, y = 0, shape = centplot), size = 5, alpha = 0.1, color = "black") +
                    scale_shape_manual(values=c(23, 15)) +
                    facet_grid(chromplot ~ . ) + 
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("") + 
                    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major.y = element_blank(), 
                        panel.grid.minor.y = element_blank())
subgenome_blocks



# infile allele specific data
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.table(trans_cent, "trans_cent.csv", sep = ",", col.names = TRUE, row.names = FALSE)

ase <- read.table("allele_specific_test_p_adjusted.csv", header = TRUE, sep = ",")
trans_cent <- read.table("trans_cent.csv", header = TRUE, sep = ",")
dim(trans_cent)
head(ase)
dim(ase)

head(trans_cent)
tail(trans_cent)
head(ase)
tail(ase)
# merge available data
ase_cent <- merge(ase, trans_cent, by.x = "gene_name", by.y = "geneID", all.x = TRUE)
dim(ase_cent)
head(ase_cent)
write.table(ase_cent, "ase-cent.csv")

ase_cent$abs_t <- abs(ase_cent$t_stat)
ase_cent <- ase_cent[!is.na(ase_cent$Chr),]
plot(ase_cent$abs_t)

?ifelse
ase_cent$large <- ifelse(ase_cent$abs_t >= 10, 1, 0)
head(ase_cent)
sum(ase_cent$large)

ase_dist1 <- ggplot(ase_cent)
ase_dist1 <- ase_dist1 +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$sub_genome == "NA"),],
                        aes(xmin = tx_start, xmax = tx_end, ymin = -10, ymax = 0, color = sub_genome), size = 4) +
                    geom_point(
                        aes(x = tx_start, y = abs_t), size = 1, alpha = 0.3, color = "black") +
                    facet_grid(Chr ~ . ) +
                    # geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("t-statistic of cis-eQTL") 
ase_dist1

ase_dist2 <- ggplot(ase_cent)
ase_dist2 <- ase_dist2 +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$sub_genome == "NA"),],
                        aes(xmin = tx_start, xmax = tx_end, ymin = -3, ymax = 3, color = sub_genome), size = 3) +
                    geom_point(
                        aes(x = tx_start, y = t_stat), size = 1, alpha = 0.3, color = "black") +
                    facet_grid(Chr ~ . ) +
                    # geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("t-statistic of cis-eQTL") 
ase_dist2

ase_dist3 <- ggplot(ase_cent)
ase_dist3 <- ase_dist3 +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$sub_genome == "NA"),],
                        aes(xmin = tx_start, xmax = tx_end, ymin = -3, ymax = 3, color = sub_genome), size = 3) +
                    geom_point(aes(x = tx_start, y = abs_t), size = 1, alpha = 0.3, color = "black") +
                    geom_smooth(aes(x = tx_start, y = abs_t), colour="darkgoldenrod1", size=1.0, method="loess", degree=0, span=0.001, se=FALSE) +
                    facet_grid(Chr ~ . ) +
                    # geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("t-statistic of cis-eQTL") 
ase_dist3



############
#load un_eqtl.RData
############
peak <- max(hots1$max.N)
hottest <- head(hots1, 50)
hottest
hotspots <- ggplot(hots1)
hotspots <- hotspots + geom_line(aes(x = pos, y = max.N), size = 0.5) +
                        facet_grid(~ chr) +
                        geom_hline(yintercept = 100, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme_bw()  +
                        theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
hotspots

head(hots1)
dim(hots1)

#convert marker names to physical positions
hots1$position <- row.names(hots1)
hots1$position <- as.numeric(sub("(A)(\\d+)(x)(\\d+)", "\\4", hots1$position))
head(hots1)
tail(hots1)
hots1$position <- hots1$position/1000000
hots1$Chr      <- hots1$chr



######
#NEAR FINAL PLOT
######
# final 3 sub genomes
ase_cent_plot <- ggplot(data = ase_cent[!(ase_cent$sub_genome == "NA"),])
ase_cent_plot <- ase_cent_plot +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$sub_genome == "NA"),],
                       aes(xmin = tx_start, xmax = tx_end, ymin = -1, ymax = 3, color = sub_genome), size = 4) +
                    geom_point(data = hots1,
                    	aes(x = position, y = max.N), size = 3, alpha = 0.3, color = "black") +
                    facet_grid(Chr ~ . ) +
                    geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("trans-eQTL number") 

ase_cent_plot


#final all ranges
ase_cent_plot <- ggplot(data = ase_cent[!(ase_cent$Block == "NA"),])
ase_cent_plot <- ase_cent_plot +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$Block == "NULL"),],
                       aes(xmin = tx_start, xmax = tx_end, ymin = -1, ymax = 3, color = Block), size = 4) +
                    geom_point(data = hots1,
                    	aes(x = position, y = max.N), size = 3, alpha = 0.3, color = "black") +
                    facet_grid(Chr ~ . ) +
                    geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("trans-eQTL number") 

ase_cent_plot

#check out flowering genes
ase_cent[ase_cent$gene_name =="Bra009055",]
ase_cent[ase_cent$gene_name =="Bra028599",]
ase_cent[ase_cent$gene_name =="Bra006051",]


# 2015_06_30
# rolling averages
head(ase_cent)
library(zoo)


ase_dist3 <- ggplot(ase_cent)
ase_dist3 <- ase_dist3 +  theme_bw() + 
                    geom_rect(data = ase_cent[!(ase_cent$sub_genome == "NA"),],
                        aes(xmin = tx_start, xmax = tx_end, ymin = -3, ymax = 3, color = sub_genome), size = 3) +
                    geom_point(aes(x = tx_start, y = abs_t), size = 1, alpha = 0.3, color = "black") +
                    geom_smooth(aes(x = tx_start, y = abs_t), colour="darkgoldenrod1", size=1.0, method="loess", degree=0, span=0.001, se=FALSE) +
                    facet_grid(Chr ~ . ) +
                    # geom_hline(yintercept = 150, color = "black", size = 1) +
                    xlab("Genomic Position of Gene Start Site (Mbp)") +
                    ylab("t-statistic of cis-eQTL") 
ase_dist3

transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)

A03 <- transcripts[transcripts$tx_chrom == "A03",]
dim(A03)

ase_cis <- merge(ase, transcripts, by.x = "gene_name", by.y = "tx_name", all.x = TRUE)
head(ase_cis)
tail(ase_cis)
dim(ase_cis)
A03 <- ase_cis[ase_cis$tx_chrom == "A03",]
dim(A03)


?rollmean
?subset
tempA03 <- subset(ase_cent, Chr == "A03")
head(tempA03)
str(tempA03)
dim(tempA03)
tempA03$pos <- as.numeric(rownames(tempA03))
length(unique(ase_cent$gene_name))


?zoo
tempseries <- zoo(tempA03$t_stat, tempA03$pos)
!unique(tempA03$tx_start)

head(tempseries)
length(tempseries)
?rollmean
plot(tempseries)
temp_roll <- rollmean(tempseries, 3, fill = NA)
head(temp_roll)
length(temp_roll)
plot(temp_roll)
x.Date <- as.Date("2003-02-01") + c(1, 3, 7, 9, 14) - 1
x <- zoo(rnorm(5), x.Date)
head(x)
plot(x)
time(x)
plot(tempA03$)

