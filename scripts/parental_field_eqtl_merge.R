###########
# Cody Markelz
# markelz@gmail.com
# Modified December 17, 2015
###########

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
dups <- cis_df[!duplicated(cis_df$tx_name),]

dim(dups)
head(dups)
tail(dups)
dups$tx_chrom

#####
# START HERE CIS eQTL here
#####



trans_df <- subset(cistrans_df, cis_trans == "trans")
trans_df <- trans_df[!grepl("^Sc", trans_df$tx_chrom),]
dim(trans_df)
head(trans_df)

# we really do not have the resolution to care about how EXACT the qtl is next
# to the start site for the cis so we should just choose largest lod
# that is the statistical information that we have for the population
test <- head(cis_df, n = 20)
test
?aggregate

test <- aggregate(lod ~ tx_name, data = test, max)
test



