rm(list=ls())
gc()

########################################################################
########### CLONE CAPSULES FROM MANUAL TO MANUAL FILES #################
# CREATE THE MANUALLY ORIENTED BASE FILES
# Orient 1 large colony per tracking system used (5 total) by hand.
# pick 1 out of this 5 and, in FortStudio, create a capsule definition for a medium sized ant and replicate the shape for all of the ants of the colony.
# Copy this capsule for the remaining 4 colonies using Clone_capsule_manual_to_manual.R . The originals of these files have been stored as *.myrmidon.old

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

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

### directory of data and myrmidon files
#dir_data <- '/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA'
dir_data <- '/media/bzniks/DISK4/ADRIANO/EXPERIMENT_DATA'
#dir_data <- '/home/cf19810/Documents/TEMP/EXPERIMENT_DATA/'


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
  REP.files       <- list.files(REP.folder, pattern = "CapsuleDef2018.myr")
  #ensure they contain metadata
  REP.files       <- REP.files[grep(pattern = "withMetaData",REP.files)]
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  for (variable in REP.files) {
    
    #get substring in variable until R9BS_
    REP_treat_name <- sub("\\_.*", "", variable)
    treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
    
    #find first folder of the experiment
    REP.data_folder      <- list.files(path=REP.folder,pattern=glob2rx(paste0("*",treat_name,"*.0000*")))
    REP.data_folder <- REP.data_folder[!grepl(REP.data_folder,pattern="myrmidon")]
    
    INFO.file <- paste(REP.folder,REP.data_folder,"artemis.INFO",sep="/")
    
    XX <- grep("Running on machine:", readLines(INFO.file), value = TRUE)
    TrackSys_name <- sub(".*Running on machine: ", "", XX) 
    
    EXP_list <- rbind(EXP_list ,data.frame(TrackSys_name,REP_treat_name,path_name=paste(REP.folder,variable,sep="/")))
    
  }
}

### manually oriented myrmidon file name
#defined_capsule_file <- paste(dir_data,"REP1/R1SP_27-02-21_Oriented.myrmidon",sep='') #polyakov
#defined_capsule_file <- paste(dir_data,"REP3/CapsuleDef3_R3SP.myrmidon",sep='') #prideaux
#defined_capsule_file <- paste(dir_data,"/REP1/R1SP_27-02-21_AntsCreated_AutoOriented_withMetaData_NS.myrmidon",sep='')
defined_capsule_file <- paste(dir_data,"/REP1/R1SP_27-02-21_AntsCreated_AutoOriented_withMetaData_NS_NS_q_CapsuleDef2018.myrmidon",sep='')


# Select the queen capsule missing files
warning("select the appropriate file here!")
no_capsule_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "CapsuleDef2018.myr")),]

### FILES TO REVIEW
warning("make sure that the files do have a calculated size for the queen or otherwise these do not get processed")
# ##TEMPORARY
# #select only R9SS, R9BS, R9SP has these where missing the measurement
#TEMP <- c("R9SS","R9BS","R9SP")
# no_capsule_list <- no_capsule_list[,"path_name"]
#no_capsule_list <- no_capsule_list[grepl(paste(TEMP,collapse="|"),no_capsule_list)]



#################################################################################################################################################################################################################
###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
#################################################################################################################################################################################################################
#data_list         <- list ("/home/eg15396/Documents/Data/NTM/NTM_s30_auto_orient.myrmidon") ###here list all the myrmidon files containing oriented data
data_list         <- list (defined_capsule_file)
oriented_metadata <- NULL
capsule_list <- list()
for (myrmidon_file in data_list){
  experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
  oriented_data <- fmExperimentOpen(myrmidon_file)
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  
  # check who's the queen
  for (ant in oriented_ants){
    individual  <- ant$ID
    IsQueen <- NULL
    for (ROW in 1:nrow(ant$getValues("IsQueen"))) {
      IsQueen <- c(IsQueen,ant$getValues("IsQueen")[ROW,"values"]) #in case there are multiple rows of the variable (it is currently the case in files generated by DataPrep3_CopyMetadata)
      if (IsQueen==TRUE) { OrientQueenID <- ant$ID }
    }
  }
  
  
  for (ant in oriented_ants){
    #get info only for queen
    if (ant$ID==OrientQueenID){
    ###extract ant length and capsules
    ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
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
                                 r2_ratio   = capsule_coord$r2[1]/ant_length_px)
      
      
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
                                                              stringsAsFactors = F))
     }
    }
  }
}


for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}

no_capsule_list <- no_capsule_list[,"path_name"]
### Write capsule data in each manually oriented file
for (no_capsule_file in no_capsule_list) {
  
  no_capsule_name <- unlist(strsplit(no_capsule_file,split="/"))[length(unlist(strsplit(no_capsule_file,split="/")))]
  no_capsule_Repl_name <- unlist(strsplit(no_capsule_name,split="_"))[1]

# open tracking data which need new capsule
tracking_data <- fmExperimentOpen(no_capsule_file) 
ants <- tracking_data$ants

#Select queen only
# check who's the queen
for (no_cap_ant in ants){
  individual  <- no_cap_ant$ID
  IsQueen <- NULL
  for (ROW in 1:nrow(no_cap_ant$getValues("IsQueen"))) {
    IsQueen <- c(IsQueen,no_cap_ant$getValues("IsQueen")[ROW,"values"]) #in case there are multiple rows of the variable (it is currently the case in files generated by DataPrep3_CopyMetadata)
    if (IsQueen==TRUE) { 
      print(no_cap_ant)
      NoCapQueenID <- no_cap_ant$ID }
  }
}


  #use mean size of each manually oriented file that needs the capsule
  ant_length_px <- mean(fmQueryComputeMeasurementFor(tracking_data,antID=ants[[NoCapQueenID]]$ID)$length_px)
  
  ants[[NoCapQueenID]]$clearCapsules()
  
  #for (caps in 1:length(capsule_list)){
  ##assign capule numbers that match the order of the looped capsule names positions
  capsule_number <- 0
   for (capsule_name in unlist(tracking_data$antShapeTypeNames)) {
     capsule_number <- capsule_number +1
    # the file information
    #MAKE SURE THERE IS CAPSULE MATCHING, TO AVOID MIXING UP SHAPE INDEXES 
    #capsule_names <-   tracking_data$antShapeTypeNames[[which(tracking_data$antShapeTypeNames[[caps]] %in% names(capsule_list))]]
    
    capsule_ratios <- capsule_list[[capsule_name]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
    capsule_coords <- ant_length_px*capsule_ratios
  
    # capsule_name  <- capsule_names[[capsules[[caps]]$type]]
    ants[[NoCapQueenID]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )

}

#tracking_data$save(no_capsule_file) 
tracking_data$save(paste0(sub("\\..*", "", no_capsule_file),"_q.myrmidon"))
print(paste0("Adding capsule for QUEEN to",no_capsule_name))

#close experiment
rm(list=(c("tracking_data")))

}#no_capsule_file

#IMPORTANT:
# WOULD BE BETTER TO DO THIS STEP AUTOMATICALLY BUT IT HAS BEEN DONE MANUALLY AT THE MOMENT:
# SAVE THE FILES WITH A NEW NAME "NAME_TrackSystemName-base.myrmidon" 
# THIS IS THE BASE INPUT FOR THE FOLLOWING STEP IN auto_orientation_loop.R

##LASTLY!! Go to fort-studio to check things look all right
