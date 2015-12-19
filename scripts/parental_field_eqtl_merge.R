###########
# Cody Markelz
# markelz@gmail.com
# Modified December 17, 2015
###########

# TODO make sure there are not two trans eQTL per chromosome per gene

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

cis_plot <- ggplot(cis_df)
cis_plot <- cis_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
cis_plot

qplot(trans_df$lod, trans_df$qtl_pos, facet = trans_df$tx_chrom)



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

trans_plot <- ggplot(trans_df)
trans_plot <- trans_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
trans_plot

qplot(trans_df$lod, trans_df$qtl_pos, facet = trans_df$tx_chrom)


head(trans_df)
head(cis_df)
str(cis_df)
str(trans_df)
ct_merge <- rbind(cis_df, trans_df)
head(ct_merge)
dim(ct_merge)

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod, color = cis_trans), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
merge_plot

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = tx_start, color = cis_trans), size = 2) +
                        facet_grid(qtl_chrom ~ qtl_chrom)
                        theme(text = element_text(size = 20))
merge_plot

