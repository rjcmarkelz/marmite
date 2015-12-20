###########
# Cody Markelz
# markelz@gmail.com
# Modified December 17, 2015
###########

# TODO make sure there are not two trans eQTL per chromosome per gene
# Go Enrichment?
# trans hotspots

library(dplyr)
library(data.table)
library(ggplot2)

# load dataset
load('~/git.repos/brassica_eqtl_v1.5/data/un_eqtl.RData')
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

?aggregate
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
setwd('~/git.repos/brassica_eqtl_v1.5/output/')
ggsave("cis_eqtl_plot.pdf", width = 10, height = 15)

# trans eQTL
trans_df <- subset(cistrans_df, cis_trans == "trans")
dim(trans_df)
trans_df <- trans_df[!grepl("^Sc", trans_df$tx_chrom),]
dim(trans_df)
# [1] 11520    13

head(trans_df)
tail(trans_df)

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
ggsave("cis_eqtl_plot.pdf", width = 4, height = 4)

head(trans_df)
head(cis_df)
str(cis_df)
str(trans_df)
ct_merge <- rbind(cis_df, trans_df)
head(ct_merge)
dim(ct_merge)
str(ct_merge)

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod, color = cis_trans), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
merge_plot

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot + geom_point(aes(x = qtl_pos, y = tx_start, color = lod ), size = 1.5) +
                        scale_y_reverse() +
                        facet_grid(tx_chrom ~ qtl_chrom) + theme_bw() + 
                        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_blank())
merge_plot


ggsave("cis_trans_eqtl_plot.pdf", width = 10, height = 10)

setwd('~/git.repos/brassica_eqtl_v1.5/data/')
br_flr <- read.delim("br_flowering_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_flr)
dim(br_flr)
colnames(br_flr)

# make trans plot with flowering time coordinates
# call each df seperately for each geom
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
ggsave("trans_flr_plot.pdf", width = 10, height = 10)

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
ggsave("cis_flr_plot.pdf", width = 10, height = 10)


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

ggsave("cis_de_fruit_plot.pdf", width = 10, height = 10)

trans_de_plot <- ggplot(trans_df)
trans_de_plot <- trans_de_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_plot
ggsave("trans_de_fruit_plot.pdf", width = 10, height = 10)


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
ggsave("cis_de_leaf_plot.pdf", width = 10, height = 10)

trans_de_leaf_plot <- ggplot(trans_df)
trans_de_leaf_plot <- trans_de_leaf_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_leaf_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_leaf_plot
ggsave("trans_de_leaf_plot.pdf", width = 10, height = 10)

# end
