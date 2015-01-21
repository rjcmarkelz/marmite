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
#
br_phys
pheno_list <- pull.pheno(br_phys)
colnames(pheno_list)

so_flowering <- scanone(br_phys, pheno.col = 4,
 method = "imp", use="all.obs")
plot(so_flowering)

so_flowering <- scanone(br_phys, pheno.col = "germ_flr", method = "imp", use="all.obs")
plot(so_flowering, chr = "A10")

so_X2012_UN_r <- scanone(br_phys, pheno.col = "X2012_UN_r", method = "imp", use="all.obs")
plot(so_X2012_UN_r)

#potential overlapping eQTL hotspots
so_X2011_STP <- scanone(br_phys, pheno.col = "X2011_STP", method = "imp", use="all.obs")
plot(so_X2011_STP, chr = "A03")

so_X2012_CR_STP <- scanone(br_phys, pheno.col = "X2012_CR_STP", method = "imp", use="all.obs")
plot(so_X2012_CR_STP, chr = "A03")

so_X2012_UN_STP <- scanone(br_phys, pheno.col = "X2012_UN_STP", method = "imp", use="all.obs")
plot(so_X2012_UN_STP, chr = "A03")

so_X2012_diff_STP <- scanone(br_phys, pheno.col = "X2012_diff_STP", method = "imp", use="all.obs")
plot(so_X2012_diff_STP, chr = "A03")

save.image(file = "brassica_BLUPs.RData", version = NULL,
 ascii = FALSE, safe = TRUE)
so.perm2 <- scanone(br_phys, method = "imp", n.perm = 1000) 
summary(so.perm2)
so.perm295 <- summary(so.perm2)[1] #keep 95%
# LOD thresholds (1000 permutations)
#      lod
# 5%  3.23
# 10% 2.81

#old new compare
brassica_traits <- read.cross("csvsr", genfile ="old_map_rqtl_missing_RILS_removed.csv", 
	                       phefile="br_blups_RQTL.csv", genotypes=c("0","2"))

class(brassica_traits)[1] <- "riself"
brassica_traits <- jittermap(brassica_traits)
brassica_traits
traits_list <- pull.pheno(brassica_traits)
colnames(traits_list)

X2011_STP
scanone_2011_STP <- scanone(brassica_traits, pheno.col = "X2011_STP", method = "imp", use="all.obs")
scanone_2011_STP
plot(scanone_2011_STP, chr = 3)


so_flower_oldmap <- scanone(brassica_traits, pheno.col = 4, method = "imp", use="all.obs")
plot(so_flower_oldmap)
plot(so_flowering)


so.perm <- scanone(brassica_traits, method = "imp", n.perm = 1000) 
summary(so.perm)
so.perm95 <- summary(so.perm)[1] #keep 95%
# LOD thresholds (1000 permutations)
#      lod
# 5%  2.44
# 10% 2.13

so_flowering <- scanone(br_phys, pheno.col = 4,
 method = "imp", use="all.obs")
plot(so_flowering)



qtl_plot <- qtl_plot +        facet_grid(~ chr, scales = "free_x", space = "free_x") +
geom_rect(aes(fill = background.color), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 0.50, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        scale_y_continuous(expand = c(0, 0), limits = c((peak * -0.06), (peak * 1.02))) +
                        theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90),
                          axis.line=element_line(),
                          panel.margin = unit(0, "cm")) +
                        ggtitle("LOD Curves for QTLs") +
                        xlab("Position in cM") +
                        ylab("LOD Score") 


