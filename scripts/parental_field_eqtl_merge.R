###########
# Cody Markelz
# markelz@gmail.com
# Modified March 17, 2016
###########

library(dplyr)
library(data.table)
library(ggplot2)

# load dataset
load('~/git.repos/brassica_eqtl_paper/input/scanone-eqtl.RData')
ls()

# choose between cis-peaks
# trans peaks go to scaffolds
# this is why I need the scaffold data from Mike
dim(cistrans_df)
head(cistrans_df)
tail(cistrans_df)

scaffolds <- cistrans_df[grepl("^Sc", cistrans_df$tx_chrom),]
dim(scaffolds)
head(scaffolds)

cis_df <- subset(cistrans_df, cis_trans == "cis")
dim(cis_df)
head(cis_df)
length(unique(cis_df$tx_name))
str(cis_df)

cis_ag <- aggregate(lod ~ tx_name, data = cis_df, max)
dim(cis_ag)

cis_df <- merge(cis_ag, cis_df, by = c("tx_name", "lod"))
dim(cis_df)
head(cis_df)
length(unique(cis_df$tx_name))
cis_df <- cis_df[!duplicated(cis_df$tx_name),]

dim(cis_df)
head(cis_df)
tail(cis_df)
cis_df$tx_chrom
str(cis_df)
cis_df$qtl_pos <- as.numeric(cis_df$qtl_pos)

cis_plot <- ggplot(cis_df)
cis_plot <- cis_plot +  theme_bw() + geom_point(aes(x = tx_start, y = lod), size = 1.5) +
                        facet_grid(qtl_chrom ~ . ) +
                        theme(text = element_text(size = 20))
cis_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_eqtl_plot.png", width = 10, height = 15)

#get large cis effect genes
# arbitrary cutoff
large_cis <- cis_df[cis_df$lod > 100,]
dim(large_cis)
large_cis
#73
setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(large_cis, "large_effect_cis.csv", sep = ",", col.names = TRUE, row.names = TRUE)


# trans eQTL
trans_df <- subset(cistrans_df, cis_trans == "trans")
dim(trans_df)
trans_df <- trans_df[!grepl("^Sc", trans_df$tx_chrom),]
dim(trans_df)
# [1] 11024    13

head(trans_df)
tail(trans_df)

large_trans <- trans_df[trans_df$lod > 100,]
dim(large_trans)
large_trans
#17
setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(large_trans, "large_effect_trans.csv", sep = ",", col.names = TRUE, row.names = TRUE)

#########
# TODO
# need to make sure that there are not 2 eQTL per chromosome
# 1 MB cutoff?
#########

# we really do not have the resolution to care about how EXACT the qtl is next
# to the start site for the cis so we should just choose largest lod
# that is the statistical information that we have for the population

trans_ag <- aggregate(lod ~ tx_name, data = trans_df, max)
dim(trans_ag)
head(trans_ag)
trans_df <- merge(trans_ag, trans_df, by = c("tx_name", "lod"))
dim(trans_df)
head(trans_df)
trans_df$qtl_pos <- as.numeric(trans_df$qtl_pos)

trans_plot <- ggplot(trans_df)
trans_plot <- trans_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
trans_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("trans_eqtl_plot.png", width = 10, height = 10)

head(trans_df)
head(cis_df)
str(cis_df)
str(trans_df)

# put data together
ct_merge <- rbind(cis_df, trans_df)
head(ct_merge)
dim(ct_merge)
str(ct_merge)

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod, color = cis_trans), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
merge_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_trans_eqtl_plot.png", width = 10, height = 10)

# cis trans plot
merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot + geom_point(aes(x = qtl_pos, y = tx_start, color = lod ), size = 1.5) +
                        scale_y_reverse() +
                        facet_grid(tx_chrom ~ qtl_chrom) + theme_bw() + 
                        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_blank())
merge_plot
?ggsave
ggsave("cis_diagonal.png", width = 10, height = 10, dpi = 300)

############
# flowering gene distribution
############

setwd('~/git.repos/brassica_eqtl_v1.5/data/')
br_flr <- read.delim("br_flowering_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_flr)
dim(br_flr)
colnames(br_flr)

str(flr_genes)
flr_genes <- br_flr[9]
colnames(flr_genes)[1] <- "gene_name"
head(flr_genes)

trans_df$tx_name <- as.character(trans_df$tx_name)
cis_df$tx_name <- as.character(cis_df$tx_name)

trans_flr_df <- trans_df[trans_df$tx_name %in% flr_genes$gene_name,]
head(trans_flr_df)

trans_flr_plot <- ggplot(trans_df)
trans_flr_plot <- trans_flr_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_flr_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_flr_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("trans_flr_plot.png", width = 10, height = 10)

cis_flr_df <- cis_df[cis_df$tx_name %in% flr_genes$gene_name,]
head(cis_flr_df)
head(trans_df)
dim(cis_flr_df)

cis_flr_plot <- ggplot(cis_df)
cis_flr_plot <- cis_flr_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_flr_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_flr_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_flr_plot.png", width = 10, height = 10)

############
# parental analysis data
############
setwd('~/git.repos/brassica_parents/data/')

br_fruits <- read.delim("parental_fruit_field_DE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_fruits)
tail(br_fruits)
dim(br_fruits)
br_fruits$tx_name <- rownames(br_fruits)

cis_de_df <- cis_df[cis_df$tx_name %in% br_fruits$tx_name,]
head(cis_de_df)
dim(cis_de_df)
dim(cis_df)

trans_de_df <- trans_df[trans_df$tx_name %in% br_fruits$tx_name,]
head(trans_de_df)
dim(trans_de_df)
dim(trans_df)

cis_de_plot <- ggplot(cis_df)
cis_de_plot <- cis_de_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_de_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_de_plot

ggsave("cis_de_fruit_plot.png", width = 10, height = 10)

trans_de_plot <- ggplot(trans_df)
trans_de_plot <- trans_de_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_plot
ggsave("trans_de_fruit_plot.png", width = 10, height = 10)


br_leaf <- read.delim("parental_leaf_field_DE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_leaf)
dim(br_leaf)
br_leaf$tx_name <- rownames(br_leaf)


cis_de_leaf_df <- cis_df[cis_df$tx_name %in% br_leaf$tx_name,]
head(cis_de_leaf_df)
dim(cis_de_leaf_df)
dim(cis_df)

trans_de_leaf_df <- trans_df[trans_df$tx_name %in% br_leaf$tx_name,]
head(trans_de_leaf_df)
dim(trans_de_leaf_df)
dim(trans_df)

cis_de_leaf_plot <- ggplot(cis_df)
cis_de_leaf_plot <- cis_de_leaf_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_de_leaf_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_de_leaf_plot
ggsave("cis_de_leaf_plot.png", width = 10, height = 10)

trans_de_leaf_plot <- ggplot(trans_df)
trans_de_leaf_plot <- trans_de_leaf_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_leaf_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_leaf_plot
ggsave("trans_de_leaf_plot.png", width = 10, height = 10)

#############
# eQTL GO enrichment
#############
library(Biostrings)
library(GO.db)
library(goseq)

# infile and manipulate go annotation file
setwd("/Users/rjcmarkelz1/git.repos/brassica_genome_db/raw_data")
brass_go <- read.table("Brassica_rapa_v1.5_final_annot.wego", header = FALSE)
head(brass_go)
colnames(brass_go) <- c("Gene", "GO")
dim(brass_go)

# make list object to store data
brass_go_list <- strsplit(as.character(brass_go[,2]), split=",", fixed=T)
head(brass_go_list)
names(brass_go_list) <- as.character(brass_go[,1])
head(brass_go_list)
length(brass_go_list)

# also need to infile the gene length data as goseq uses this to estimate any biases based on gene
# lengths in the dataset. 
brass_gene_lengths <- read.table("Brassica_rapa_v1.5_final_gene_lengths", header = FALSE)
head(brass_gene_lengths)
colnames(brass_gene_lengths) <- c("Gene", "length")
dim(brass_gene_lengths)
str(brass_gene_lengths)
# [1] 43463     2
brass_gene_lengths$Gene <- as.character(brass_gene_lengths$Gene)
dim(brass_gene_lengths)
# [1] 43463     2

gene_names <- as.data.frame(high1$names)
head(gene_names)
names(gene_names) <- "tx_name"
keeps <- as.character(gene_names$tx_name)
# subset brass_gene_lengths for eQTL genes
# cis df
brass_gene_lengths_red <- brass_gene_lengths[brass_gene_lengths$Gene %in% keeps,]
str(brass_gene_lengths_red)
# 35039 total genes tested for eQTL

# of total genes tested need to get the ones that are significant for whatever
keeps_cis <- as.character(cis_df$tx_name)
head(keeps_cis)
cis_gene_lengths <- brass_gene_lengths_red[brass_gene_lengths_red$Gene %in% keeps_cis,]
not_cis_gene <- brass_gene_lengths_red[!brass_gene_lengths_red$Gene %in% keeps_cis,]

dim(cis_gene_lengths)
head(cis_gene_lengths)
# [1] 8907    2
dim(not_cis_gene)
head(not_cis_gene)
# [1] 26132     2

nullp_vector <- rep(c(1,0),c(nrow(cis_gene_lengths), nrow(not_cis_gene)))
head(nullp_vector)
tail(nullp_vector)
names(nullp_vector) <- c(cis_gene_lengths$Gene, not_cis_gene$Gene)

bias_data_vector <- c(cis_gene_lengths$length, not_cis_gene$length)
head(bias_data_vector)
# ?rep
# ?nullp

brass_nullp <- nullp(nullp_vector, genome = NULL, id = NULL, bias.data = bias_data_vector)
# # Warning message:
# # In pcls(G) : initial point very close to some inequality constraints

# ?goseq
# # need to decide whether or not to use the useuse_genes_without_cat=TRUE option. 
# # ~8000 genes not included in analysis that do not have a GO category classification
go_analysis_cis  <-  goseq(brass_nullp, gene2cat = brass_go_list, use_genes_without_cat = TRUE)
head(go_analysis_cis, 100)

setwd('~/git.repos/brassica_eqtl_v1.5/output/')
write.table(go_analysis_cis, "cis_eqtl_enrichment.csv", sep = ",", col.names = TRUE, row.names = TRUE)

##########
#Trans
##########

# of total genes tested need to get the ones that are significant for whatever
keeps_trans <- as.character(trans_df$tx_name)
head(keeps_trans)
length(keeps_trans)
# [1] 3749
trans_gene_lengths <- brass_gene_lengths_red[brass_gene_lengths_red$Gene %in% keeps_trans,]
not_trans_gene <- brass_gene_lengths_red[!brass_gene_lengths_red$Gene %in% keeps_trans,]

dim(trans_gene_lengths)
head(trans_gene_lengths)
# [1] 3596    2

dim(not_trans_gene)
head(not_trans_gene)
# [1] 31443     2

nullp_vector <- rep(c(1,0),c(nrow(trans_gene_lengths), nrow(not_trans_gene)))
head(nullp_vector)
tail(nullp_vector)
names(nullp_vector) <- c(trans_gene_lengths$Gene, not_trans_gene$Gene)

bias_data_vector <- c(trans_gene_lengths$length, not_trans_gene$length)
head(bias_data_vector)

brass_nullp <- nullp(nullp_vector, genome = NULL, id = NULL, bias.data = bias_data_vector)
# # Warning message:
# # In pcls(G) : initial point very close to some inequality constraints

# ?goseq
# # need to decide whether or not to use the useuse_genes_without_cat=TRUE option. 
# # ~8000 genes not included in analysis that do not have a GO category classification
go_analysis_trans  <-  goseq(brass_nullp, gene2cat = brass_go_list, use_genes_without_cat = TRUE)
head(go_analysis_trans, 100)

setwd('~/git.repos/brassica_eqtl_v1.5/output/')
write.table(go_analysis_trans, "trans_eqtl_enrichment.csv", sep = ",", col.names = TRUE, row.names = TRUE)

setwd("/Users/rjcmarkelz1/git.repos/brassica_eqtl_v1.5/data")
save.image(file = "un_eqtl_parent_field.RData", version = NULL, ascii = FALSE, safe = TRUE)

# end
