
###########
# cody markelz 
# quick hacking to make some progress
# April 30, 2015
###########


#######
# 2010
#######
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
brass_2010 <- read.table("brassica_field_2010_blups.csv", sep = ",", header = TRUE)
head(brass_2010)
brass_2010$Line <- sub("(L_)(\\d+)", "RIL_\\2", brass_2010$Line)

#######
# 2012
#######
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
brass_2012 <- read.table("Brassica2012allmeans.csv", sep = ",", header = TRUE)
head(brass_2012)
brass_2012$id <- sub("(\\d+)", "RIL_\\1", brass_2012$id)

###############
# merge 2010 2012
###############
field_data <- merge(brass_2010, brass_2012, by.x = "Line", by.y = "id", all = TRUE)
dim(field_data)

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")
write.table(field_data, "2010_2012_field_data.csv", sep = ",")