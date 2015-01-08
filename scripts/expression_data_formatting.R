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

br_blues_total <- br_blues_10[,grepl("*\\UN$",names(br_blues_10))]
dim(br_blues_total)

#transpose
br_blues_t <- as.data.frame(t(br_blues_sub))
head(br_blues_t)[,1:10]

br_total_t <- as.data.frame(t(br_blues_total))
head(br_total_t)[,1:10]


#add ril numbers for RQTL
br_blues_t$id <- sub("(Br_group)(\\d+)(_)(CR)", "RIL_\\2", row.names(br_blues_t))
head(br_blues_t)[,35020:35023]

br_total_t$id <- sub("(Br_group)(\\d+)(_)(UN)", "RIL_\\2", row.names(br_total_t))
head(br_total_t)[,35020:35040]

#transpose back
#takes a bit to transpose back in this direction
br_blues_final <- as.data.frame(t(br_blues_t))
br_blues_total <- as.data.frame(t(br_total_t))

head(br_blues_final)[,1:10]
tail(br_blues_final)[,1:10]

head(br_blues_total)[,1:10]
tail(br_blues_total)[,1:10]

#save output
write.table(br_blues_final, "br_blues_RQTL.csv", col.names= FALSE, 
	         row.names = TRUE, sep = ",")

write.table(br_blues_total, "br_blues_total_RQTL.csv", col.names= FALSE, 
	         row.names = TRUE, sep = ",")