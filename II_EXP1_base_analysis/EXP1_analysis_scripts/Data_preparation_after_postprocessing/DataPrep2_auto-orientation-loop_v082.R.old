rm(list=ls())
gc()

########################################################################
################# ORIENT ALL FILES (AntsCreated) #######################
# Once DataPrep1_Clone-capsule-manual-to-manual_v082.R is run, check in fort-studio that capsules are present and normal and that Queen infos are overwritten (tag size, manual orientation, manual capsules).
# In auto_orientation_loop.R, load these 5 oriented and capsule provided files and orient all of the other files

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(data.table)
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

###source C++ movement direction program
sourceCpp("/media/cf19810/DISK4/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
#sourceCpp("/media/eg15396/DISK4/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")

### directory of data and myrmidon files
dir_data <- '/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA'
#dir_data <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA"



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

#exclude files already auto-oriented
EXP_list <- EXP_list[which(!grepl(EXP_list$path_name,pattern = "AutoOrient")),]

# flag which ones are the base oriented files
EXP_list$OrientedCapsule = ifelse(grepl("-base.myrmidon",EXP_list$path_name),"true","false")

### manually oriented ref file name
ref_orient_caps_file <- EXP_list[which(EXP_list$OrientedCapsule=="true"),]

# feed in the analysis of auto-orientation all AntsCreated files.
#ToOrient_file <- EXP_list[which(EXP_list$OrientedCapsule=="false"),]
ToOrient_file <- EXP_list[which(grepl(EXP_list$path_name,pattern = "AntsCreated.myrmidon")),]

# Select the metadata-rich files for the Queen ID check
Metadata_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "DeathRecord_NoOrient|ManOriented")),]
Metadata_list <- Metadata_list[which(!grepl(Metadata_list$path_name,pattern = "base")),]

### Loop through all the directories in the dir_folder

# loop through the unique Tracking systems
#TS <- "westerby" #
for (TS in unique(ref_orient_caps_file$TrackSys_name)){
    print(paste("START: Auto-orient + capsule for",TS, sep =" "))
   
  #################################################################################################################################################################################################################
  ###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
  ###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
  ###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
  #################################################################################################################################################################################################################
  #data_list         <- list ("/home/eg15396/Documents/Data/NTM/NTM_s30_auto_orient.myrmidon") ###here list all the myrmidon files containing oriented data
  data_list         <- ref_orient_caps_file[which(ref_orient_caps_file$TrackSys_name==TS),"path_name"]
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
  
  ### Measures of mean ant length and offset between tag centre and ant centre will be heavily influenced by the queen #########
  ### So we need to remove the queen from the computation
  ### One way of doing so is to find and remove outliers in the ant length measure (provided there is enough variation between queen and worker size)
  interquartile_range <- quantile(oriented_metadata$length_px,probs=c(0.25,0.75))
  outlier_bounds      <- c(interquartile_range[1]-1.5*(interquartile_range[2]-interquartile_range[1]),interquartile_range[2]+1.5*(interquartile_range[2]-interquartile_range[1]))
  ###apply outlier exclusion to oriented_metadata...
  oriented_metadata <- oriented_metadata[which(oriented_metadata$length_px>=outlier_bounds[1]&oriented_metadata$length_px<=outlier_bounds[2]),]  
  ###...and to capsule list
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <-capsule_list[[caps]] [ as.character(interaction(capsule_list[[caps]] $experiment,capsule_list[[caps]] $antID))%in%as.character(interaction(oriented_metadata $experiment,oriented_metadata $antID)),]
  }
  
  
  ###Once queen(s) has(have) been removed, get the mean coordinates of the offset between tag centre and ant centre
  mean_x_ant_coord <- mean(oriented_metadata$x_ant_coord)
  # mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
  mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation
  ###Furthermore, get the average worker length from the data
  mean_worker_length_px <- mean(oriented_metadata$length_px)
  ###Finally, get information on each capsule
  for (caps in 1:length(capsule_list)){
    capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
  }
  
  ###################################################################################
  ############### open the experiment files to ORIENT ###############################
  ToOrient_data_list         <- ToOrient_file[which(ToOrient_file$TrackSys_name==TS),"path_name"]
 
for (ToOrient_myr_file in ToOrient_data_list){
    #ToOrient_myr_file <- ToOrient_data_list[1] #temp
  
  # if the _AutoOriented file file doesn't exist, then continue
  if ( !file.exists(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon"))) {
    
    ToOrient_exp_name <- unlist(strsplit(ToOrient_myr_file,split="/"))[length(unlist(strsplit(ToOrient_myr_file,split="/")))]
    print(paste("Assign Orientation and Capsule to",ToOrient_exp_name, sep =" "))
    tracking_data <- fmExperimentOpen(ToOrient_myr_file)
   
    ## The AntsCreated files 
    ToOrient_ants <- tracking_data$ants
    ToOrient_capsule_names <- tracking_data$antShapeTypeNames
    
# ### DELETE ALL POSSIBLE BASE INFO NOT REQUIRED FROM THE FILES TO ORIENT
#If ever processing files with already assigned measurements, make sure to remove the neasurements (capules, sposes, etc) before assigning new ones
#     # delete individuals' capsule data IF PRESENT
#     for (ant in ToOrient_ants){
#       ToOrient_ants[[ant$ID]]$clearCapsules()
#     }
#     
#     # delete the capsule shapes IF PRESENT
# if (length(ToOrient_capsule_names)>0) {
#     for (caps in 1:length(ToOrient_capsule_names)){
#       tracking_data$deleteAntShapeType(caps)
#     }
# }
#     
# 
#     #delete the ant pose (it doesn't seem to change the output of the script)
#     for (i in 1:length(ToOrient_ants)){
#       ##delete pose
#       for (id in ToOrient_ants[[i]]$identifications){
#         id$clearUserDefinedAntPose()
#       }
#     }
#       
#   
#       
    # # check print identifications
    # for (a in ToOrient_ants) {
    #  printf("Ant %s is identified by:\n", fmFormatAntID(a$ID))
    #  for (i in a$identifications){
    #    printf(" * %s\n", capture.output(i))
    #  }
    # }
    
    ###create capsule list
    ###CAUTION the relationship between shape name and id may be different from your manually oriented files.
    ###Hence for post-processing analyses it will be important to use the capsule names rather than IDs
    ###Also, it would be safer not to use manually annotated files in the analysis anyway for consistency 
    ### (it would not do to anaylse some colonies with precise manually annotated data and other with approximate automated data)
    ###So you need to create new automatically oriented myrmidon files for all your colonies including the manually-oriented ones
    for (caps in 1:length(capsule_list)){
      tracking_data$createAntShapeType(names(capsule_list)[caps])
    }
    
    ###get trajectory data to extract ant orientation (if tracking is too long perhaps only use 24 hours?)
    ## IT Crashes with times longer than 6 hours.....
    from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) ###experiment start time
    to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start + 12*3600  ) ###experiment start time
    
    ###Then compute trajectories
    max_gap <- fmHour(24*365)  ###this parameter is very important to - use a super large value to make sure you get only one trajectory per ant!!!
    positions <- fmQueryComputeAntTrajectories(tracking_data,start = from,end = to,maximumGap = max_gap,computeZones = TRUE)
    
    #hard-wire ant correspondence between trajectories_summary and trajectories
    positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
    names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
    
    ##define a max temporal gap for which you are happy to calculate a movement angle; e.g. 0.5 s
    max_time_gap <- 0.2 #0.5 ##LOWER THAN 0.2 BREAKS THE CODE
    
    ##define a minimum distance moved, as you don't want to use noise or small shifts in position in this calculation; e.g. 30 pix (to think about)
    min_dist_moved <- 15 #30
    
    # open Metadata experiment to get Queen info
    Metadata_exp <- fmExperimentOpen(Metadata_list[which(sub("\\_.*", "", basename(ToOrient_myr_file)) == Metadata_list$REP_treat_name),"path_name"] ) 
    Metadata_ants <- Metadata_exp$ants
    # check who's the queen
    for (ant in Metadata_ants){
      individual  <- ant$ID
      IsQueen <- NULL
      for (ROW in 1:nrow(ant$getValues("IsQueen"))) {
        IsQueen <- c(IsQueen,ant$getValues("IsQueen")[ROW,"values"]) #in case there are multiple rows of the variable (it is currently the case in files generated by DataPrep3_CopyMetadata)
        if (IsQueen==TRUE) { QueenID <- ant$ID }
      }
    }
    #check that it is correct (IsQueen value = TRUE)
    Metadata_ants[[QueenID]]$getValues("IsQueen")
    rm(list=c("Metadata_exp","Metadata_ants"))
    
    for (i in 1:length(ToOrient_ants)){
      #skip Queen
      if (!ToOrient_ants[[i]]$ID==QueenID) {
      #check that the ant trajectory exists
      aux <- positions$trajectories_summary$antID==ToOrient_ants[[i]]$ID
      if (sum(aux)==0) {
        print(paste('ANT',i,"HAS NO TRAJECTORY. SHE MAY BE DEAD OR UNTAGGED",sep=" "))
        next}
        #print(paste0("ANT ", ToOrient_ants[[i]]$ID))
      ####to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
      traj <- positions$trajectories [[   positions$trajectories_summary[which(aux),"antID_str"]    ]]
      
      
      ###### THIS MESSES UP WITH THE ORIENTATION................
      # ##### cut traj until metadata IsAlive is FALSE
      # IsAlive <-  ToOrient_ants[[i]]$getValues("IsAlive")
      # #check if there is a death time (IsAlive=FALSE)
      #   if(any(grepl(FALSE,IsAlive$values))){
      #     print(paste0("Ant's ",i," traj cut until death time"))
      #     death_time <-   IsAlive$times[!is.na(IsAlive$times)]
      #     #subtract timeStart from death time
      #     MAX_traj <- as.numeric(as.POSIXlt(death_time) - fmQueryGetDataInformations(tracking_data)$start, units="secs")
      #     #CUT TRAJ WHEN TIME > OF DEATH TIME
      #     traj <- traj[which(traj$time < MAX_traj),]
      # }else{
      #   #print("good little ant that didn't die")
      # }
  
      
      ###feed traj to c++ program
      traj <- cbind(traj,add_angles(traj,max_time_gap,min_dist_moved))
      
      
      ## get mean deviation angle between body and tag - the ant angle is equal to minus the Tag minus Movement angle output by C++ program
      AntAngle <- as.numeric(- mean(circular(na.omit(traj$Tag_minus_Movement_Angle),units="radians",zero=0)))
      ##now use trigonometry to calculate the pose, using AntAngle
      x_tag_coord <- mean_x_ant_coord*cos(AntAngle) - mean_y_ant_coord*sin(AntAngle)
      y_tag_coord <- mean_x_ant_coord*sin(AntAngle) + mean_y_ant_coord*cos(AntAngle)
      
      ##write this into ant metadata
      for (id in ToOrient_ants[[i]]$identifications){
        id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
      }
      
      ##also add this to trajectories_summary
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ToOrient_ants[[i]]$ID),"ant_angle"] <- AntAngle
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ToOrient_ants[[i]]$ID),"x_tag_coord"] <- x_tag_coord
      positions$trajectories_summary[which(positions$trajectories_summary$antID==ToOrient_ants[[i]]$ID),"y_tag_coord"] <- y_tag_coord
      
      ###finally, for each ant, add capsules using mean_ant_length and capsule_list
      for (caps in 1:length(capsule_list)){
        capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
        capsule_coords <- mean_worker_length_px*capsule_ratios
        
        
        ToOrient_ants[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
        
      }#ADD CAPSULES
    }# Skip queen
  }# LOOP ANTS
    
  tracking_data$save(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon"))
  print(paste("Saving",ToOrient_exp_name, "as AutoOriented.myrmidon",sep =" "))
  rm(list=(c("tracking_data"))) #remove experiment
 }else{print(paste0(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon")," already exists! Skip!"))} # ToOrient FILES LOOP #closes: for (ToOrient_myr_file in ToOrient_data_list){
  rm(list=(c("oriented_data"))) #remove experiment
} # IF EXISTS, SKIP PROCESSING
} # TS ID LOOP # closes: for (TS in unique(EXP_list$TrackSys_name))

cat("LOOP ENDED!! \n Go to fort-studio to check things look all right
\n AND OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)
")

##LASTLY!! Go to fort-studio to check things look all right
### AND TO OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)
