#script to change 
#Mean_ant_length_colonies.txt (originated from script 1_CAPSULE_ASSIGNMENT_Linda)
#to
#Mean_ant_length_per_TrackingSystem.txt (for EXP1_base_analysis)
usr <- "lsartori"
SCRIPTDIR <- paste("/media",usr,"Seagate Portable Drive/circadian_rhythm_2022_experiment/scripts/EXP1_base_analysis/EXP1_analysis_scripts",sep="/") # "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
BEH_FUNCTIONS <-  paste(SCRIPTDIR,"Behavioural_Inference",sep="/")
BODYLENGTH_FILE_L <- paste(BEH_FUNCTIONS,"Mean_ant_length_colonies.txt", sep = "/")
BODYLENGTH_FILE_A <- paste(BEH_FUNCTIONS,"ADRIANO_Mean_ant_length_per_TrackingSystem.txt", sep = "/")

Mean_ant_length_colonies <- read.table(BODYLENGTH_FILE_L, header = T, stringsAsFactors = F, sep = ",") ### body length information
str(Mean_ant_length_colonies)
Mean_ant_length_colonies_A <- read.table(BODYLENGTH_FILE_A, header = T, stringsAsFactors = F, sep = ",") ### body length information
str(Mean_ant_length_colonies_A)

#change column names
colnames(Mean_ant_length_colonies)  <- c("mean_worker_length_px", "mean_worker_length_mm", "myrmidon_file")

#get tracking system ("TS") from myrmidon

for (i in c(1:8)){
  Mean_ant_length_colonies$TS[i] <- unlist(strsplit(Mean_ant_length_colonies$myrmidon_file[i],split="_"))[1]
  
}
str(Mean_ant_length_colonies)

output_name <- file.path(paste0(BEH_FUNCTIONS, "/Mean_ant_length_per_TrackingSystem.txt"))
write.table(Mean_ant_length_colonies,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
