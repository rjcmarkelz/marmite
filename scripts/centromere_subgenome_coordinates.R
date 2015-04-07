library(ggplot2)
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

#figure 2 found it in another paper
# trans_cent$sub_genome[trans_cent$Block %in% c("U", "Tb", "Qb", "Xa", "D", "Ta", "S", "J", "Fa")] <- "LF"
# #notice Ma block is both LF and MF2 on chrm A06
# trans_cent$sub_genome[trans_cent$Block %in% c("Cb", "Ab", "Ba", "Ma", "Xb", "Qa", "Wa", "L", "K", "V")] <- "LF"
# trans_cent$sub_genome[trans_cent$Block %in% c("H", "G", "Fb", "Xc", "Ib", "E",)] <- "LF"
# trans_cent$sub_genome

# not finished yet
trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = y_hold, group = Block, color = Block, y = Mbp)) +
                        facet_grid(~ tx_chrom) 

trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = Mbp, y = y_hold, color = Block)) +
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
trans_cent_plot <- trans_cent_plot +  theme_bw() + geom_tile(aes(x = Mbp, y = y_hold, color = sub_genome)) +
                        facet_grid(~ tx_chrom) 

trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_tile(data = trans_cent[(trans_cent$tx_chrom == "A09") &
                     !(trans_cent$sub_genome == "NA"),], 
	                 aes(x = Mbp, y = y_hold, color = sub_genome)) 

trans_cent_plot

na.omit(data)


, tx_chrom = "A01"










