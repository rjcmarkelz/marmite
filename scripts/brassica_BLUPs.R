setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
br_blups <- read.delim("Brapa_Height_BLUPs.csv", 
	                               header = TRUE, row.names = 1, sep = ",")
br_blups2 <- read.delim("2011_phenology_heights_blups.csv", 
	                               header = TRUE, row.names = 1, sep = ",")

head(br_blups)
head(br_blups2)
#add ril numbers for RQTL
row.names(br_blups) <- sub("(\\d+)", "RIL_\\1", row.names(br_blups))
head(br_blups)

# difference mapping 
br_blups_diff <- br_blups[,grepl("*\\UN*",names(br_blups))] - br_blups[,grepl("*\\CR*",names(br_blups))] 
                   
colnames(br_blups_diff) <- c("2012_diff_dur", "2012_diff_hmax", "2012_diff_Inflect_DD", 
	                          "2012_diff_Inflect_size", "2012_diff_r", "2012_diff_STP")

head(br_blups_diff)
br_blups <- cbind(br_blups, br_blups_diff)
head(br_blups)


br_blups$id <- row.names(br_blups)
head(br_blups)

br_blups2$id <- row.names(br_blups2)
head(br_blups2)

dim(br_blups)
dim(br_blups2)
br_blups2
?merge
br_blups_total <- merge(br_blups2, br_blups, by = "id", all= TRUE)
head(br_blups_total)
dim(br_blups_total)



# just curious how automated blup caller was doing in R vs. SAS
# take a look at the lmer vs. the SAS output. Not exactly perfect.
# Adjusted R-squared:  0.9807, but do not know the exact SAS model
# also do not know if Rob used same QC steps
# 
plot(br_blups_total$plant_ht_6_final, br_blups_total$X2011_Hmax)
lmer_SAS <- lm(plant_ht_6_final ~ X2011_Hmax, data = br_blups_total)
summary(lmer_SAS)
plot(lmer_SAS)

# move id to last column
br_blups_total <- br_blups_total[,c(2:38,1)]
head(br_blups_total)

br_blups_final <- as.data.frame(t(br_blups_total))
head(br_blups_final)
tail(br_blups_final)

write.table(br_blups_final, "br_blups_RQTL.csv", col.names= FALSE, 
	         row.names = TRUE, sep = ",")


library(qtl)

br_phys <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="br_blups_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = c("-","NA"))

class(br_phys)[1] <- "riself"
br_phys <- jittermap(br_phys)
br_phys

#about a minute
br_phys <- calc.errorlod(br_phys, error.prob=0.001)
pheno_list <- pull.pheno(br_phys)
colnames(pheno_list)

so_flowering <- scanone(br_phys, pheno.col = 4,
 method = "imp", use="all.obs")
plot(so_flowering)

so_flowering <- scanone(br_phys, pheno.col = "germ_flr", method = "imp", use="all.obs")
plot(so_flowering)

so_X2012_UN_r <- scanone(br_phys, pheno.col = "X2012_UN_r", method = "imp", use="all.obs")
plot(so_X2012_UN_r)








