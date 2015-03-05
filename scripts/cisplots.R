library(ggplot2)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

ase <- read.table("allele_specific_test_p_adjusted.csv", header = TRUE, sep = ",")
head(ase)
dim(ase)

#infile genomic coordinates of genes
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")

transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)


ase_cis <- merge(ase, transcripts, by.x = "gene_name", by.y = "tx_name", all.x = TRUE)
head(ase_cis)
tail(ase_cis)
dim(ase_cis)

ase_cis$abs_t <- abs(ase_cis$t_stat)
head(ase_cis,100)
qplot(ase_cis$abs_t)

#make megabase
ase_cis$Mbp <- ase_cis$tx_start/1000000

#alleles
ase_cis_plot <- ggplot(ase_cis)
ase_cis_plot <- ase_cis_plot +  theme_bw() + geom_point(aes(x = Mbp, y = t_stat,  alpha = 0.1), size = 2) +
                        facet_grid(~ tx_chrom) +
                        geom_hline(yintercept = 10, color = "red", size = 1) +
                        geom_hline(yintercept = -10, color = "red", size = 1) +
                        geom_hline(yintercept = 0, color = "black", size = 1) +
                        ggtitle("t-statistic Chromosome Distributions") +
                        xlab("Genomic Position of Gene Start Site (Mbp)") +
                        ylab("t-statistic") +
                        theme(axis.title=element_text(face="bold",size="20"),
                              axis.text=element_text(face="bold", size="12"),
                              legend.position = "none")  
ase_cis_plot

#grouped
ase_cis_abs <- ggplot(ase_cis)
ase_cis_abs <- ase_cis_abs +  theme_bw() + geom_point(aes(x = Mbp, y = abs_t, alpha = 0.1), size = 2) +
                        facet_grid(~ tx_chrom) +
                        geom_hline(yintercept = 10, color = "red", size = 1) +
                        ggtitle("t-statistic Chromosome Distributions") +
                        xlab("Genomic Position of Gene Start Site (Mbp)") +
                        ylab("t-statistic") +
                        theme(axis.title=element_text(face="bold",size="20"),
                              axis.text=element_text(face="bold", size="12"),
                              legend.position = "none")                      
ase_cis_abs

