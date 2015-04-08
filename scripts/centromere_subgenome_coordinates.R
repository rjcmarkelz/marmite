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



?replace
trans_cent$Centromere[1:150]
trans_cent$centplot <- trans_cent$Centromere
trans_cent$centplot 

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

trans_cent$Block <- sub("Wa", "W", trans_cent$Block)
trans_cent$Block <- sub("Wb", "W", trans_cent$Block)



trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                   geom_tile(data = trans_cent[!(trans_cent$Block == "NULL"),],
                    aes(x = Mbp, y = placeholder, color = Block)) +
                   geom_point(data = trans_cent[!(trans_cent$centplot == "NA"),],
                   	aes(x = Mbp, y = placeholder, fill = centplot)) +
                    facet_grid(tx_chrom ~ .) 

trans_cent_plot


trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_tile(data = trans_cent[!(trans_cent$Block == "NULL"),],
                    aes(x = Mbp, y = placeholder, color = Block)) +
                   geom_point(data = trans_cent[!(trans_cent$centplot == "NA"),],
                   	aes(x = Mbp, y = placeholder, shape = centplot)) +
                    scale_shape_manual(values=c(21,24)) + 
                    facet_grid( ~ chromplot ) 

trans_cent_plot

trans_cent_plot <- ggplot(trans_cent)
trans_cent_plot <- trans_cent_plot +  theme_bw() + 
                    geom_tile(data = trans_cent[!(trans_cent$Block == "NULL"),],
                       aes(x = Mbp, y = placeholder, color = Block)) +
                    geom_point(aes(x = Mbp, y = placeholder, shape = centplot), size = 4, alpha = 0.1, color = "black") +
                    scale_shape_manual(values=c(20, 15)) +
                    facet_grid(chromplot ~ . )

trans_cent_plot


