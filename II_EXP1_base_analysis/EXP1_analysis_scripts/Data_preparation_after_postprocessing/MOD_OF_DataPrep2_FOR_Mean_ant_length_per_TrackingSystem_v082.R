rm(list=ls())
gc()
mallinfo::malloc.trim(0L)

########################################################################
################# ORIENT ALL FILES (tracking_data) #######################
# Once DataPrep1_Clone-capsule-manual-to-manual_v082.R is run, check in fort-studio that capsules are present and normal and that Queen infos are overwritten (tag size, manual orientation, manual capsules).
# In auto_orientation_loop.R, load these 5 oriented and capsule provided files and orient all of the other files

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(data.table)
library(circular)
library(R.utils)
library(dplyr)

#### FUNCTIONS
#list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

# SWITCHES
#check if there are multiple detections caused by corrupted files/tag duplication/etc
CHECK_MULTIPLE_DETECT <- FALSE
#set to true if you want to obtain the Mean_ant_length_per_TrackingSystem.txt file
#set to false if you want to obtain the ant_length_per_Colony.txt file
LENGTH_SUMMARY        <- FALSE

###source C++ movement direction program
#sourceCpp("/media/cf19810/DISK4/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
#sourceCpp("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
#sourceCpp("/media/eg15396/DISK4/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
#sourceCpp("/home/bzniks/Downloads/PhD-exp1-data-analysis-main/scriptsR/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")

### directory of data and myrmidon files
#dir_data <- '/media/bzniks/DISK3/ADRIANO/EXPERIMENT_DATA'
#dir_data <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA"
dir_data <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA"
#dir_data  <- "/media/cf19810/Seagate Portable Drive/ADRIANO/EXPERIMENT_DATA_EXTRAPOLATED"


#### ACCESS FILES
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(dir_data, n = 1)

#select REP folders
files_list <- files_list[grep("REP",files_list)]

#### OPEN REPLICATE
EXP_list <- NULL
# replicate folder
for (REP.n in 1:length(files_list)) {
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = ".myrmidon")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  for (variable in REP.files) {
    if (!grepl("*CapsuleDef*",variable)) {
    #get substring in variable until R9BS_
    REP_treat_name <- sub("\\_.*", "", variable)
    treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
    
    #find first folder of the experiment
    REP.data_folder      <- list.files(path=REP.folder,pattern=glob2rx(paste0("*",treat_name,"*.0000*")))
    REP.data_folder <- REP.data_folder[!grepl(REP.data_folder,pattern="myrmidon")]
    
    INFO.folder <- paste(REP.folder,REP.data_folder,sep="/")
    INFO.file <- list.files(path=INFO.folder,pattern=glob2rx(paste0("artemis.INFO*")))
      
    INFO.file <-  paste(INFO.folder,INFO.file,sep="/")
    
    XX <- grep("Running on machine:", readLines(INFO.file), value = TRUE)
    TrackSys_name <- sub(".*Running on machine: ", "", XX) 
    
    EXP_list <- rbind(EXP_list ,data.frame(TrackSys_name,REP_treat_name,path_name=paste(REP.folder,variable,sep="/")))
    }
  }
}


# # Select the metadata-rich files
# Metadata_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "DeathRecord_NoOrient|ManOriented.myr")),]
# Metadata_list <- Metadata_list[which(!grepl(Metadata_list$path_name,pattern = "base")),]

#exclude files already auto-oriented
EXP_list <- EXP_list[which(!grepl(EXP_list$path_name,pattern = "AutoOrient")),]

# flag which ones are the base oriented files
EXP_list$OrientedCapsule = ifelse(grepl("-base.myrmidon",EXP_list$path_name),"true","false")
#EXP_list$OrientedCapsule = ifelse(grepl("ManOriented_CapsuleDef3.myrmidon",EXP_list$path_name),"true","false") # To use as base the CapsuleDef3!

### manually oriented ref file name
ref_orient_caps_file <- EXP_list[which(EXP_list$OrientedCapsule=="true"),]

# feed in the analysis of auto-orientation all tracking_data files.
#ToOrient_file <- EXP_list[which(EXP_list$OrientedCapsule=="false"),]
ToOrient_file <- EXP_list[which(grepl(EXP_list$path_name,pattern = "AntsCreated.myrmidon")),]

#check how many have base files (unoriented)
base_files <- EXP_list[which(grepl(EXP_list$path_name,pattern = "_base.myrmidon")),]

CHECKDATA <- data.frame(REP_treat_name=unique(EXP_list$REP_treat_name),base_oriented=FALSE, base=FALSE, AntsCreated=FALSE, MultipleDetect=NA)

# flag which ones are the base oriented files
EXP_list$OrientedCapsule <- ifelse(grepl("-base.myrmidon",EXP_list$path_name),"true","false")
CHECKDATA$base_oriented[grep(paste(EXP_list[which(EXP_list$OrientedCapsule %in% c("true")),"REP_treat_name"],collapse="|"),CHECKDATA$REP_treat_name)] <- "TRUE"
# flag which ones are Man_oriented files
EXP_list$ManOriented <- ifelse(grepl("*ManOriented.myrmidon",EXP_list$path_name),"true","false")
ref_ManOriented_file <- EXP_list[which(EXP_list$ManOriented=="true"),]

# flag which ones are the base files
EXP_list$base <- ifelse(grepl("_base.myrmidon",EXP_list$path_name),"true","false")
CHECKDATA$base[grep(paste(EXP_list[which(EXP_list$base %in% c("true")),"REP_treat_name"],collapse="|"),CHECKDATA$REP_treat_name)] <- "TRUE"
# flag which ones have the AntCreated files
EXP_list$AntsCreated <- ifelse(grepl("AntsCreated.myrmidon",EXP_list$path_name),"true","false")
CHECKDATA$AntsCreated[grep(paste(EXP_list[which(EXP_list$AntsCreated %in% c("true")),"REP_treat_name"],collapse="|"),CHECKDATA$REP_treat_name)] <- "TRUE"

### Loop through all the directories in the dir_folder

# loop through the unique Tracking systems
#TS <- "karla" #
LENGTH_OUTPUT <- NULL
to_keep_1 <- c(ls(),"to_keep_1","TS", "length_info")


for (TS in unique(ref_orient_caps_file$TrackSys_name)){
  print(paste("START: ant size for",TS, sep =" "))
  #################################################################################################################################################################################################################
  ###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
  ###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
  ###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
  #################################################################################################################################################################################################################
  #data_list         <- list ("/home/eg15396/Documents/Data/NTM/NTM_s30_auto_orient.myrmidon") ###here list all the myrmidon files containing oriented data
  if (LENGTH_SUMMARY){
    data_list         <- ref_orient_caps_file[which(ref_orient_caps_file$TrackSys_name==TS),"path_name"]
  }else{
    data_list         <- ref_ManOriented_file[which(ref_ManOriented_file$TrackSys_name==TS),"path_name"]
  }
  
  oriented_metadata <- NULL
  capsule_list <- list()
  for (myrmidon_file in data_list){
    experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
    print(paste("AntPose and Capsules extraction for",experiment_name, sep =" "))
    oriented_data <- fmExperimentOpen(myrmidon_file)
    oriented_ants <- oriented_data$ants
    capsule_names <- oriented_data$antShapeTypeNames
    for (ant in oriented_ants){
      ###extract ant length and capsules
      ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
      ant_length_mm <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_mm)
      if (LENGTH_SUMMARY){
      capsules      <- ant$capsules
      for (caps in 1:length(capsules)){
        capsule_name  <- capsule_names[[capsules[[caps]]$type]]
        capsule_coord <- capsules[[caps]]$capsule
        capsule_info <- data.frame(experiment = experiment_name,
                                   antID      = ant$ID,
                                   c1_ratio_x = capsule_coord$c1[1]/ant_length_px,
                                   c1_ratio_y = capsule_coord$c1[2]/ant_length_px,
                                   c2_ratio_x = capsule_coord$c2[1]/ant_length_px,
                                   c2_ratio_y = capsule_coord$c2[2]/ant_length_px,
                                   r1_ratio   = capsule_coord$r1[1]/ant_length_px,
                                   r2_ratio   = capsule_coord$r2[1]/ant_length_px
        )
        
        if (!capsule_name %in%names(capsule_list)){ ###if this is the first time we encounter this capsule, add it to capsule list...
          capsule_list <- c(capsule_list,list(capsule_info)) 
          if(length(names(capsule_list))==0){
            names(capsule_list) <- capsule_name
          }else{
            names(capsule_list)[length(capsule_list)] <- capsule_name
          }
        }else{###otherwise, add a line to the existing dataframe within capsule_list
          capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
        }
      }
      }#length summary
      
      ###extract offset btewen tag centre and ant centre
      for (id in ant$identifications){
        oriented_metadata <- rbind(oriented_metadata,data.frame(experiment       = experiment_name,
                                                                antID            = ant$ID,
                                                                tagIDdecimal     = id$tagValue,
                                                                angle            = id$antAngle,
                                                                x_tag_coord      = id$antPosition[1], 
                                                                y_tag_coord      = id$antPosition[2],
                                                                x_ant_coord      = id$antPosition[1]*cos(-id$antAngle) - id$antPosition[2]*sin(-id$antAngle),
                                                                y_ant_coord      = id$antPosition[1]*sin(-id$antAngle) + id$antPosition[2]*cos(-id$antAngle),
                                                                length_px        = ant_length_px,
                                                                length_mm        = ant_length_mm,
                                                                stringsAsFactors = F))
      }
    }
  }
  
  #the file may not include only oriented ants
  oriented_metadata <- oriented_metadata[which(!is.na(oriented_metadata$length_px)),]
  
  ### Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen #########
  ### So we need to remove the queen from the computation
  ### One way of doing so is to find and remove outliers in the ant length measure (provided there is enough variation between queen and worker size)
  interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ###apply outlier exclusion to oriented_metadata...
  oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  
  ###...and to capsule list
  if (LENGTH_SUMMARY){
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
  }
    ###Once queen(s) has(have) been removed, get the mean coordinates of the offset between tag centre and ant centre
    mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)
    # mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
    mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation
    ###Furthermore, get the average worker length from the data
    mean_worker_length_px <- mean(oriented_metadata$length_px)
    
    length_info <- data.frame(mean_worker_length_px=mean_worker_length_px,mean_worker_length_mm=mean_worker_length_mm, TS=TS,myrmidon_file=myrmidon_file)
  }else{
    length_info <- data.frame(worker_length_px=oriented_metadata$length_px, worker_length_mm=oriented_metadata$length_mm, TS=TS,myrmidon_file=oriented_metadata$experiment)
    length_info$myrmidon_file <- unlist(strsplit(length_info$myrmidon_file,split="/"))
    length_info$REP_treat <- sub("\\_.*", "", length_info$myrmidon_file)
  }
  

  # ###################################################################################
  # ############### list the experiment files to CHECK ###############################

  ### CHECK MULTIPLE DETECTIONS
  if (CHECK_MULTIPLE_DETECT) {
    ToOrient_data_list         <- ToOrient_file[which(ToOrient_file$TrackSys_name==TS),"path_name"]
    # 
    to_keep_2 <- c(ls(),"to_keep_2","ToOrient_myr_file")
    
  for (ToOrient_myr_file in ToOrient_data_list){
    REP_treat_name <- sub("\\_.*", "", basename(ToOrient_myr_file))
    
    #ToOrient_myr_file <- ToOrient_data_list[1] #temp
print(ToOrient_myr_file)
    # # if the _AutoOriented file file doesn't exist, then continue
    # if ( !file.exists(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented_withMetaData_NS.myrmidon"))) {
      ###get experiment and replicate name, and identify the corresponding metadata file
      ToOrient_exp_name <- unlist(strsplit(ToOrient_myr_file,split="/"))[length(unlist(strsplit(ToOrient_myr_file,split="/")))]
      ToOrient_Repl_name <- unlist(strsplit(ToOrient_exp_name,split="_"))[1]
      #Metadata_myr_file <- Metadata_list[which(Metadata_list$REP_treat_name==ToOrient_Repl_name),]$path_name
      #print(paste("Assign Orientation, Capsules and Metadata to",ToOrient_exp_name, sep =" "))

      ###open the myrmidon file to orient nd the metadata myrmidon file
      #Metadata_exp <- fmExperimentOpen(Metadata_myr_file)
      #Metadata_ants <- Metadata_exp$ants
      ToOrient_data <- fmExperimentOpen(ToOrient_myr_file)
      # fmQueryGetDataInformations(ToOrient_data)$details
      #N of multiple detections
      CHECKDATA["MultipleDetect"][which(CHECKDATA$REP_treat_name==REP_treat_name),] <-   sum(fmQueryComputeTagStatistics(ToOrient_data)$multipleSeen)
      
      to_keep_3 <- c(ls(),"to_keep_3","ant_INDEX")
    gc()
  } # IF EXISTS, SKIP PROCESSING
}

  LENGTH_OUTPUT <- rbind(LENGTH_OUTPUT,length_info)
  
  
  rm(list=ls()[which(!ls()%in%to_keep_1)])
  gc()
} # TS ID LOOP # closes: for (TS in unique(EXP_list$TrackSys_name))

if (CHECK_MULTIPLE_DETECT) {
CHECKDATA$ToCheck <- "no"
CHECKDATA["ToCheck"][which(CHECKDATA$REP_treat_name %in% c("R9SS","R9BS")),] <- "yes"
write.table(CHECKDATA,file=file.path(dir_data,"EXP1_REPS_data_check_12-09-22.txt"),append=F,col.names=T,row.names=F,quote=T,sep=",")
}

if (LENGTH_SUMMARY) {
  #summary used to adjust mean_ant size in behavioural inference
  write.table(LENGTH_OUTPUT,file=file.path(dir_data,"Mean_ant_length_per_TrackingSystem.txt"),append=F,col.names=T,row.names=F,quote=T,sep=",")
}else{
  #full report of collected measurments to check differences between colonies BIG and SMALL
  write.table(LENGTH_OUTPUT,file=file.path(dir_data,"ant_length_per_Colony.txt"),append=F,col.names=T,row.names=F,quote=T,sep=",")
}


