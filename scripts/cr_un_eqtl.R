############
# the following block of code is to format the data for RQTL
# for difference mapping to look at genotype by environment interactions
############
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_blues <- read.delim("brassica_blues.csv", 
	                               header = TRUE, row.names = 1, sep = ",")


# start with a test file to make sure we can subtract columns for difference mapping
# 
test <- head(br_blues)[,1:20]
test

?grepl
test2 <- test[,grepl("*\\CR$",names(test))] - test[,grepl("*\\UN$",names(test))]
test2
dim(test)
dim(test2)
# works

#add 10 to DF so it is recentered at 10 and all values are above 0
#subtract columns from one another
br_blues_10 <- br_blues + 10
head(br_blues_10)[,1:20]
min(br_blues_10)

br_blues_sub <- br_blues_10[,grepl("*\\CR$",names(br_blues_10))] - br_blues_10[,grepl("*\\UN$",names(br_blues_10))]
head(br_blues_sub)[,1:20]
dim(br_blues_10)
dim(br_blues_sub)
# looks good

#transpose
br_blues_t <- as.data.frame(t(br_blues_sub))
head(br_blues_t)[,1:10]

#add ril numbers for RQTL
br_blues_t$id <- sub("(Br_group)(\\d+)(_)(CR)", "RIL_\\2", row.names(br_blues_t))
head(br_blues_t)[,35020:35023]

#transpose back
#takes a bit to transpose back in this direction
br_blues_final <- as.data.frame(t(br_blues_t))

head(br_blues_final)[,1:10]
tail(br_blues_final)[,1:10]

#save output
write.table(br_blues_final, "br_blues_RQTL.csv", col.names= FALSE, 
	         row.names = TRUE, sep = ",")



###########
###########
#RQTL MAPPING
###########
###########
library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

#takes a minute
brassica_genes <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="br_blues_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brassica_genes)

class(brassica_genes)[1] <- "riself"
brassica_genes <- jittermap(brassica_genes)
brassica_genes

brassica_genes <- est.rf(brassica_genes)
plot.rf(brassica_genes) 

#about a minute
brassica_genes <- calc.errorlod(brassica_genes, error.prob=0.001)

system.time(scanone.imp.1 <- scanone(brassica_genes, pheno.col = 1:35039,
 method = "imp", use="all.obs"))
#    user  system elapsed 
# 845.002  11.957 854.041 

# saved output in data directory






