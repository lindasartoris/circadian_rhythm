rm(list=ls())
gc()

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
library(reader)

##define a max temporal gap for which you are happy to calculate a movement angle; e.g. 0.5 s
max_time_gap <- 0.2 #0.5 ##LOWER THAN 0.2 BREAKS THE CODE

##define a minimum distance moved, as you don't want to use noise or small shifts in position in this calculation; e.g. 30 pix (to think about)
min_dist_moved <- 15 #30

trajectory_length_hours <- 12 ###duration of trajectory (in hours) on which to base automated orientation

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
#sourceCpp("/media/cf19810/DISK4/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
sourceCpp("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/determine_angle_automatically/Get_Movement_Angle.cpp")
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
  # REP.n <- 1
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
    
    INFO.folder <- paste(REP.folder,REP.data_folder,sep="/")
    INFO.file <- list.files(path=INFO.folder,pattern=glob2rx(paste0("*artemis.INFO*")))
    
    INFO.file <-  paste(INFO.folder,INFO.file,sep="/")
    

    XX <- grep("Running on machine:", readLines(INFO.file), value = TRUE)
    TrackSys_name <- sub(".*Running on machine: ", "", XX) 
    
    EXP_list <- rbind(EXP_list ,data.frame(TrackSys_name,REP_treat_name,path_name=paste(REP.folder,variable,sep="/")))
    
  }
}


# Select the metadata-rich files
Metadata_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "DeathRecord_NoOrient|ManOriented.myr")),]
Metadata_list <- Metadata_list[which(!grepl(Metadata_list$path_name,pattern = "base")),]

#exclude files already auto-oriented
AlreadyDone <- EXP_list[which(grepl(EXP_list$path_name,pattern = "NS_NS_q.myrmidon")),]
EXP_list <- EXP_list[which(!grepl(EXP_list$path_name,pattern = "AutoOrient")),]

# flag which ones are the base oriented files
EXP_list$OrientedCapsule = ifelse(grepl("R3SP_13-03-21_AntsCreated_AutoOriented_withMetaData_CapDef3.myrmidon",EXP_list$path_name),"true","false")
#EXP_list$OrientedCapsule = ifelse(grepl("ManOriented_CapsuleDef3.myrmidon",EXP_list$path_name),"true","false") # To use as base the CapsuleDef3!

### manually oriented ref file name
ref_orient_caps_file <- EXP_list[which(EXP_list$OrientedCapsule=="true"),]

# feed in the analysis of auto-orientation all tracking_data files.
#ToOrient_file <- EXP_list[which(EXP_list$OrientedCapsule=="false"),]
ToOrient_file <- EXP_list[which(grepl(EXP_list$path_name,pattern = "AntsCreated.myrmidon")),]



### Loop through all the directories in the dir_folder

# loop through the unique Tracking systems
#TS <- "polyakov" #
to_keep_1 <- c(ls(),"to_keep_1","TS")
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
      if (oriented_ants[[29]]$ID==29) {
        
      }else{
      capsules      <- ant$capsules}
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
  ############### list the experiment files to ORIENT ###############################
  ToOrient_data_list         <- ToOrient_file[which(ToOrient_file$TrackSys_name==TS),"path_name"]
  
  to_keep_2 <- c(ls(),"to_keep_2","ToOrient_myr_file")
  for (ToOrient_myr_file in ToOrient_data_list){
    #ToOrient_myr_file <- ToOrient_data_list[9] #temp
    
    # if the _AutoOriented file file doesn't exist, then continue
    if ( !file.exists(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon"))) {
   # if ( !file.exists(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented_withMetaData_NS_NS_q.myrmidon"))) {
      ###get experiment and replicate name, and identify the corresponding metadata file
      ToOrient_exp_name <- unlist(strsplit(ToOrient_myr_file,split="/"))[length(unlist(strsplit(ToOrient_myr_file,split="/")))]
      ToOrient_Repl_name <- unlist(strsplit(ToOrient_exp_name,split="_"))[1]
      Metadata_myr_file <- Metadata_list[which(Metadata_list$REP_treat_name==ToOrient_Repl_name),]$path_name
      print(paste("Assign Orientation, Capsules and Metadata to",ToOrient_exp_name, sep =" "))
      
      
      ###open the myrmidon file to orient nd the metadata myrmidon file
      Metadata_exp <- fmExperimentOpen(Metadata_myr_file) 
      Metadata_ants <- Metadata_exp$ants
      ToOrient_data <- fmExperimentOpen(ToOrient_myr_file)
      
      # check who's the queen
      for (ant in Metadata_ants){
        individual  <- ant$ID
        IsQueen <- NULL
        for (ROW in 1:nrow(ant$getValues("IsQueen"))) {
          IsQueen <- c(IsQueen,ant$getValues("IsQueen")[ROW,"values"]) #in case there are multiple rows of the variable (it is currently the case in files generated by DataPrep3_CopyMetadata)
          if (IsQueen==TRUE) { QueenID <- ant$ID }
        }
      }

      
      
      ############# MAKE SURE ISALIVE METADATA KEY IS PRESENT #########################
      if (is.null(Metadata_exp$metaDataKeys$IsAlive)){
      Metadata_exp$setMetaDataKey("IsAlive",TRUE)
        }

      ############# ASSIGN METADATA KEYS to tracking data #########################
      list_keys <- list()
      #assign metadata keys
      for (KEY in   1:length(Metadata_exp$metaDataKeys)) {
        key <- names(Metadata_exp$metaDataKeys[KEY])
        defaultvalue <- unname(Metadata_exp$metaDataKeys[KEY][[1]])
        ToOrient_data$setMetaDataKey(key,defaultvalue)
        #check
        #ToOrient_data$metaDataKeys[KEY]
        list_keys <- c(list_keys,key)
      }

      # ############# ASSIGN ZONES #########################
      # #assign zones
      # for (space_index in 1:length(Metadata_exp$spaces)){
      #   list_zones <- NULL
      #   for (ZONE in   1:length(Metadata_exp$spaces[[space_index]]$zones)) {
      #     zone <- Metadata_exp$spaces[[space_index]]$zones[[ZONE]]$name
      #     #defaultvalue <- unname(Metadata_exp$metaDataKeys[ZONE][[1]])
      #     ToOrient_data$spaces[[space_index]]$createZone(zone)
      #     #check
      #     #ToOrient_data$spaces[[space_index]]$zones
      #     list_zones <- c(list_zones,zone)
      #   }
      #   
      #   #Assign SHAPE to ZONE (geometry to the nest zone)
      #   #extract metadata info per key
      #   for (ZONE.1 in   1:length(ToOrient_data$spaces[[space_index]]$zones)) {
      #     #for (ZONE_KEY in list_zones) {
      #     zone_definition <- Metadata_exp$spaces[[space_index]]$zones[[ZONE.1]]$definitions
      #     #assign shapes
      #     ToOrient_data$spaces[[space_index]]$zones[[ZONE.1]]$addDefinition(shapes=zone_definition[[1]][["shapes"]],start= fmTimeSinceEver(),end=fmTimeForever())
      #   }#ZONE.1
      # }
      
      ###create capsule list
      ###CAUTION the relationship between shape name and id may be different from your manually oriented files.
      ###Hence for post-processing analyses it will be important to use the capsule names rather than IDs
      ###Also, it would be safer not to use manually annotated files in the analysis anyway for consistency 
      ### (it would not do to anaylse some colonies with precise manually annotated data and other with approximate automated data)
      ###So you need to create new automatically oriented myrmidon files for all your colonies including the manually-oriented ones
      for (caps in 1:length(capsule_list)){
        ToOrient_data$createAntShapeType(names(capsule_list)[caps])
      }
      ToOrient_capsule_names <- ToOrient_data$antShapeTypeNames
      
      ###POSITIONS for all, cut trajs after
      FROM <- fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$start) ###experiment start time
      TO   <- fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$start + 12*3600  ) ###experiment start time
      max_gap <- fmHour(24*365)  ###this parameter is very important to - use a super large value to make sure you get only one trajectory per ant!!!
      positions <- fmQueryComputeAntTrajectories(ToOrient_data,start = FROM,end = TO,maximumGap = max_gap,computeZones = TRUE)
      
      
      ###copy individual metadata, orient ant and add capsules
      to_keep_3 <- c(ls(),"to_keep_3","ant_INDEX")
      for (ant_INDEX in 1:length(ToOrient_data$ants)){   ####LOOP OVER ANTS
        ###get corresponding ant from metadata
        metadata_ant_INDEX <- which(unlist(lapply(Metadata_ants,function(x)x$ID))==ToOrient_data$ants[[ant_INDEX]]$ID )
        
        ### delete all identification for that ant from tracking data, to clear
        for (identif in rev(1:length(ToOrient_data$ants[[ant_INDEX]]$identifications))){
          ToOrient_data$deleteIdentification( ToOrient_data$ants[[ant_INDEX]]$identifications[[identif]]  )
        }
        
        ### then create exactly the same identifications as in metadata
        for (identif in 1:length(Metadata_ants[[metadata_ant_INDEX]]$identifications)){
          ToOrient_data$addIdentification(antID= Metadata_ants[[metadata_ant_INDEX]]$identifications[[identif]]$targetAntID,
                                          tagID= Metadata_ants[[metadata_ant_INDEX]]$identifications[[identif]]$tagValue,
                                          start= Metadata_ants[[metadata_ant_INDEX]]$identifications[[identif]]$start,
                                          end=   Metadata_ants[[metadata_ant_INDEX]]$identifications[[identif]]$end)
          
        }
        
        #copy metadata info for each key
        for (METADATA_KEY in list_keys) {
          if (METADATA_KEY != "Exposed"){ ### not exposed
            for (ROW in 1:nrow(Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY))) {
              if ( is.na(Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)[ROW,"times"])){ 
                ToOrient_data$ants[[ant_INDEX]]$setValue(key=METADATA_KEY,value=Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeSinceEver()) #fmTimeSinceEver()
              }else{
                ToOrient_data$ants[[ant_INDEX]]$setValue(key=METADATA_KEY,value=Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeCreate(offset=Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)[ROW,"times"]))
              }#assign value
            }#ROW       
            
          }else{ ###exposed
            if (T%in%Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)$values){
              exp_index <- which(Metadata_ants[[metadata_ant_INDEX]]$getValues(METADATA_KEY)$values)
              ToOrient_data$ants[[ant_INDEX]]$setValue(key=METADATA_KEY,value=TRUE,time=fmTimeSinceEver()) #fmTimeSinceEver()
            }
          }
          

        }#METADATA_KEY
        
 
        ###orient ant and define capsule: loop over identifications
        ###SKIP FOR QUEEN SO THAT YOU CAN ORIENT THE QUEEN AND CREATE CAPSULES MANUALLY IN FORTSTUDIO
        if (ToOrient_data$ants[[ant_INDEX]]$ID!=QueenID){
          to_keep_4 <- c(ls(),"to_keep_4","identif")
          for (identif in 1:length(ToOrient_data$ants[[ant_INDEX]]$identifications)){
            
            ###get trajectory data to extract ant orientation (if tracking is too long perhaps only use 24 hours?)
            
            ###first define start time. Should be the maximum time between beginning of experiment and beginning of identification
            ### to do so we will build a from_to_choose_from_fmTime containing the different possible times in fmTime format
            ###                    and a from_to_choose_from_value containing the corresponding times in numeric format
            ###and choose the largest
            from_to_choose_from_fmTime <- list( 
              fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$start) ###start of experiment
              ,
              ToOrient_data$ants[[ant_INDEX]]$identifications[[1]]$start           ###start of identification
            )  
            from_to_choose_from_value  <- c(
              fmQueryGetDataInformations(ToOrient_data)$start 
              ,
              as.POSIXct(capture.output(print(ToOrient_data$ants[[ant_INDEX]]$identifications[[1]]$start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) 
            )
            ###and define from as the maximum time between the two
            from <- from_to_choose_from_fmTime[[which.max(from_to_choose_from_value)]]
            from_value <- as.POSIXct(capture.output(print(from)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
            
            ### Similarly, define end time for the trajectory. 
            ### Should be the minimum time between 1/ from + trajectory_length_hours hours ; 2/end of experiment; 3/ end of identification;  and 4/ death
            ###initialise to_to_choose_from_fmTime and to_to_choose_from_value with from + trajectory_length_hours hours, end of experiment and  end of identification
            
            to_to_choose_from_fmTime <- list( 
              fmTimeCreate(offset=from_value + trajectory_length_hours*3600 ) ###1/ from + trajectory_length_hours hours
              ,
              fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$end  )            ###2/ end of experiment
              ,
              ToOrient_data$ants[[ant_INDEX]]$identifications[[1]]$end                        ###3/ end of identification
            )  
            to_to_choose_from_value  <- c(
              from_value + trajectory_length_hours*3600
              ,
              fmQueryGetDataInformations(ToOrient_data)$end 
              ,
              as.POSIXct(capture.output(print(ToOrient_data$ants[[ant_INDEX]]$identifications[[1]]$end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
            )
            ##check if ant died and if so add death time
            if ((F%in%ToOrient_data$ants[[ant_INDEX]]$getValues("IsAlive")$values)){
              to_to_choose_from_fmTime <- c(to_to_choose_from_fmTime, fmTimeCreate(offset=ToOrient_data$ants[[ant_INDEX]]$getValues("IsAlive")[which(!ToOrient_data$ants[[ant_INDEX]]$getValues("IsAlive")$values),"times"]))
              to_to_choose_from_value  <- c(to_to_choose_from_value,ToOrient_data$ants[[ant_INDEX]]$getValues("IsAlive")[which(!ToOrient_data$ants[[ant_INDEX]]$getValues("IsAlive")$values),"times"])
              
            }
            to <- to_to_choose_from_fmTime[[which.min(to_to_choose_from_value)]]
            
            
            ###Then compute trajectory
            #positions <- fmQueryComputeAntTrajectories(ToOrient_data,start = from,end = to,maximumGap = max_gap,computeZones = TRUE, matcher = fmMatcherAntID(ToOrient_data$ants[[ant_INDEX]]$ID))
            ###this should return a single trajectory for the ant of interest and the identification of interest
            
            
            #check that the ant trajectory exists
            aux <- positions$trajectories_summary$antID==ToOrient_data$ants[[ant_INDEX]]$ID
            if (sum(aux)==0) {
              print(paste('ANT',identif,"HAS NO TRAJECTORY. SHE MAY BE DEAD OR UNTAGGED",sep=" "))
              next}
            #hard-wire ant correspondence between trajectories_summary and trajectories
            positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
            names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
            
            #print(paste0("ANT ", ToOrient_data$ants[[i]]$ID))
            ####to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
            traj <- positions$trajectories [[   positions$trajectories_summary[which(aux),"antID_str"]    ]]
            
            #create time column with time since positions + start_time
            positions$trajectories[[ToOrient_data$ants[[ant_INDEX]]$ID]]$time_UNIX <- positions$trajectories[[ToOrient_data$ants[[ant_INDEX]]$ID]]$time + fmQueryGetDataInformations(ToOrient_data)$start
            #use this time to isolate relevant time-span
            #Reconvert convert Zulu time to GMT... (TEMP FIX)
            from_UNIX  <- as.POSIXct(capture.output(print(from)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
            to_UNIX  <- as.POSIXct(capture.output(print(to)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
            
            traj <- traj [ which(positions$trajectories[[ToOrient_data$ants[[ant_INDEX]]$ID]]$time_UNIX >= from_UNIX &
                                 positions$trajectories[[ToOrient_data$ants[[ant_INDEX]]$ID]]$time_UNIX <= to_UNIX),]
            
            
            ###feed traj to c++ program
            traj <- cbind(traj,add_angles(traj,max_time_gap,min_dist_moved))
            
            ## get mean deviation angle between body and tag - the ant angle is equal to minus the Tag minus Movement angle output by C++ program
            AntAngle <- as.numeric(- mean(circular(na.omit(traj$Tag_minus_Movement_Angle),units="radians",zero=0)))
            ##now use trigonometry to calculate the pose, using AntAngle
            x_tag_coord <- mean_x_ant_coord*cos(AntAngle) - mean_y_ant_coord*sin(AntAngle)
            y_tag_coord <- mean_x_ant_coord*sin(AntAngle) + mean_y_ant_coord*cos(AntAngle)
            
            ##write this into ant metadata
            ToOrient_data$ants[[ant_INDEX]]$identifications[[identif]]$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
            
            rm(list=ls()[which(!ls()%in%to_keep_4)])
            gc()
          }## identif
          ###finally, once ant has been oriented for all identifications, add its capsules
          ##assign capule numbers that match the order of the looped capsule names positions
          capsule_number <- 0
          for (capsule_name in unlist(ToOrient_data$antShapeTypeNames)) {
            capsule_number <- capsule_number +1
            #MAKE SURE THERE IS CAPSULE MATCHING, TO AVOID MIXING UP SHAPE INDEXES
            capsule_ratios <- capsule_list[[capsule_number]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
            capsule_coords <- ant_length_px*capsule_ratios
            
            ToOrient_data$ants[[ant_INDEX]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
            
            }  #ADD CAPSULES #IF THIS RESULTS IN A ERROR, CHECK PREVIOUS VERSION ON GIT OR COMPARE WITH DataPrep4_Clone-capule-queens-only_v082.R
        
            #} ## if no trajectory (eg. dead from)
        }###if not queen
        rm(list=ls()[which(!ls()%in%to_keep_3)])
        gc()
      }# LOOP ANTS
      
      ToOrient_data$save(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon"))
      print(paste("Saving",ToOrient_exp_name, "as AutoOriented.myrmidon",sep =" "))
      rm(list=(c("ToOrient_data"))) #remove experiment
    }else{print(paste0(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon")," already exists! Skip!"))} # ToOrient FILES LOOP #closes: for (ToOrient_myr_file in ToOrient_data_list){
    
    rm (list = ls() [which(!ls()%in%to_keep_2)] )
    gc()
  } # IF EXISTS, SKIP PROCESSING
  
  rm(list=ls()[which(!ls()%in%to_keep_1)])
  gc()
} # TS ID LOOP # closes: for (TS in unique(EXP_list$TrackSys_name))

cat("LOOP ENDED!! \n Go to fort-studio to check things look all right
\n AND OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)
")

##LASTLY!! Go to fort-studio to check things look all right
### AND TO OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)
