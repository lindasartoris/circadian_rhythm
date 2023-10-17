# how many ants moved outside?
rm(list=ls())

usr = "lsartori"
WORKDIR <- paste("/media",usr,"LS_1/circadian_rhythm_2022_experiment",sep="/")
DATADIR <- paste(WORKDIR, "tracking", sep = "/")
metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2022_2023-07-03.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")


item_df <- split(metadata, metadata$REP_treat)
#par(mfrow = 7, 7)

for (i in 1:length(item_df)) {
  hist(item_df[[i]]$prop_time_outside,
       main = paste(unique(item_df[[i]]$REP_treat)),
       xlab = "x")
}

item_df[[1]]$prop_time_outside
hist(item_df[[1]]$prop_time_outside,
     main = paste(unique(item_df[[i]]$REP_treat)),
     xlab = "x")


#check histogram -> what is a good cut off?
#add extra column for moved_outside
#add extra column for N_exposed_alive