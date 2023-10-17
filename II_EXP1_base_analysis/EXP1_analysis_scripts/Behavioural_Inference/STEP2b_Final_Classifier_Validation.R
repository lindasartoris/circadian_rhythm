#### THIS VERSION IS FORT 0.8.2 COMPATIBLE ####

# Script created by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson, with contributions by Enrico Gavagnign

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral



#clean start
rm(list=ls())
gc()

###parameter to set at start. To prepare data...
USER                       <- "Nathalie"#"Adriano"
training_set               <- "both" ##Adriano, Vasudha or both
BEH                        <- "G"
FRAME_RATE                 <- 8
proportion_training_events <- 0.8 ##PROPORTION OF ANNOTATED TIME WE WANT TO USE FOR TRAINING THE CLASSIFIER

###...and to read selected classifier.
measure         <- "Fbeta" ##""Fbeta" "CSI"                               ####which measure to do the selection on
priority        <- "overall" ##~trends overall
beta            <- 1      
if (measure=="Fbeta"){
  measure <- paste(measure,beta,sep="_")
}
if (BEH =="G"){
  interactions_of_interest <- list(c("head","body"))
  ###if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))
}

###############################################################################
###### GLOSSARY ###############################################################
###############################################################################

# MAN         : manual interactions/trajectories deriving from the hand annotated data (from the file annotations)
# AUTO        : automatically extracted interactions/trajectories deriving from fmQueryComputeAntInteractions
# REPLICATE   : nest ID, either "R3SP" or "R9SP" (SP: small pathogen)
# PERIOD      : the treatment period, either "pre" or "post" pathogen exposure
# REP_PER     : each of the 4 blocks analised, crossing the REPLICATE and PERIOD

###############################################################################
###### LOAD LIBRARIES AND FUNCTIONS #####################################
###############################################################################

####### navigate to folder containing myrmidon file
WORKDIR <- "/media/bzniks/Seagate\ Portable\ Drive/ADRIANO/Ants_behaviour_analysis_EXTRAPOLATED"
DATADIR <- "/media/bzniks/Seagate\ Portable\ Drive/ADRIANO/EXPERIMENT_DATA_EXTRAPOLATED"
SCRIPTDIR <- "~/Dropbox/SeniorLectureship_Bristol/Students_postdocs/PhD_students/2019 Adriano Wanderlingh/code/PhD-exp1-data-analysis-main/ScriptsR_FINAL"
SAVEOUTPUT <- "/media/bzniks/FiveTB"
BODYLENGTH_FILE <- paste(WORKDIR,"Data","/Mean_ant_length_per_TrackingSystem.txt",sep="/")
MachineLearningOutcome_DIR <- file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL")

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"BEH_Extract_movement_variables_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_PCA_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_self_defined_functions.R",sep="/"))
source(paste(SCRIPTDIR,"interaction_detection.R",sep="/"))
suppressMessages(source(paste(SCRIPTDIR,"BEH_libraries.R",sep="/")))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(SCRIPTDIR,"add_angles.cpp",sep="/"))
sourceCpp(paste(SCRIPTDIR,"merge_interactions.cpp",sep="/"))

duplicate_annotations <- function(annotation_subset,time_limit,total_duplicated_events){
  ###time_limit will work well if those two events don't overlap
  ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
  ###first list events to duplicate (if any)
  to_duplicate  <- annotation_subset[which(  (as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
  
  ###perform duplication if necessary
  if (nrow(to_duplicate)>=1){
    
    total_duplicated_events <- total_duplicated_events+nrow(to_duplicate)
    
    to_keep_as_is <- annotation_subset[which(!(as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
    to_duplicate_before             <-  to_duplicate
    to_duplicate_before$T_stop_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )
    to_duplicate_before$duration    <- as.numeric(to_duplicate_before$T_stop_UNIX - to_duplicate_before$T_start_UNIX)
    
    to_duplicate_after              <-  to_duplicate
    to_duplicate_after$T_start_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )+1/FRAME_RATE
    to_duplicate_after$duration     <- as.numeric(to_duplicate_after$T_stop_UNIX - to_duplicate_after$T_start_UNIX)
    
    annotation_subset <- rbind(to_keep_as_is,to_duplicate_before,to_duplicate_after)
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
  }
  
  return(list(annotation_subset=annotation_subset,total_duplicated_events=total_duplicated_events))
}

modify_columns <- function(annotations,BEH){
  ##rename some columns
  names(annotations)[which(grepl("rep",names(annotations)))]    <- "REPLICATE"
  names(annotations)[which(names(annotations)=="period"|names(annotations)=="treatment")] <- "PERIOD"
  
  ###specify behaviour
  if ("Behaviour"%in%names(annotations)){
    annotations$Behaviour     <- as.character(annotations$Behaviour)
  }else{
    annotations$Behaviour     <- BEH
  }
  
  ###modify class of some columns
  annotations$Actor         <- as.character(annotations$Actor)
  annotations$Receiver      <- as.character(annotations$Receiver)
  
  ###define new time columns called T_start_UNIX and T_stop_UNIX and delete T_start and T_stop
  annotations$T_start_UNIX  <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations$T_stop_UNIX   <- as.POSIXct(annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations$T_start <- NULL
  annotations$T_stop <- NULL
  return(annotations)
}

Vasudha_split <- function(annotations_Vasudha,focal_list_Vasudha,time_window_Vasudha,annotations_Adriano,BEH){
  set.seed(3)
  annotations_Vasudha <- within(annotations_Vasudha, full_Actor_ID    <- paste(REPLICATE,Actor   ,sep="_"))
  annotations_Vasudha <- within(annotations_Vasudha, full_Receiver_ID <- paste(REPLICATE,Receiver   ,sep="_"))
  
  full_ant_list_Vasudha <- NULL
  ###loop over unique combinations of REPLICATE and PERIOD
  for (REPLICATE in unique(time_window_Vasudha$REPLICATE)) 
  {
    ###############################################################################
    ###### OPEN EXPERIMENT INFORMATION ############################################
    ###############################################################################
    
    ## locate the ant info file for REPLICATE
    MyrmidonCapsuleFiles <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
    MyrmidonCapsuleFiles <- MyrmidonCapsuleFiles[which(grepl("myrmidon",MyrmidonCapsuleFiles))]
    
    if (length(MyrmidonCapsuleFiles)>0){
      MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl("CapsuleDef3",MyrmidonCapsuleFiles)]
    }else{
      MyrmidonCapsuleFiles <- list.files(path=file.path(DATADIR,gsub( "R", "REP"  ,     substr(REPLICATE,1,nchar(REPLICATE)-2))), pattern=REPLICATE, full.names=T)
      MyrmidonCapsuleFiles <- MyrmidonCapsuleFiles[which(grepl("myrmidon",MyrmidonCapsuleFiles))]
      MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl("CapsuleDef3",MyrmidonCapsuleFiles)]
    }
    e <- fmExperimentOpen(MyrmidonCapsuleFile)
    replicate_full_ant_list <- unlist(lapply(e$ants,function(x)x$ID))
    full_ant_list_Vasudha <- c(full_ant_list_Vasudha,paste(REPLICATE,replicate_full_ant_list,sep="_"))
    rm(list=c("e"))
    gc()
  }
  
  full_ant_list_Adriano <- unique(with(annotations_Adriano[which(!is.na(annotations_Adriano$Focal)),]   ,paste(REPLICATE,Focal   ,sep="_")))
  
  ###draw one colony that will go in full into the test dataset, and another which will go in full into the training dataset
  pre_colony_training     <- sample(unique(annotations_Vasudha$REPLICATE),1)
  pre_colony_test <- unique(annotations_Vasudha$REPLICATE)[which(unique(annotations_Vasudha$REPLICATE)!=pre_colony_training)]
  
  ###and exchange for the post window
  post_colony_test     <- pre_colony_training
  post_colony_training <- pre_colony_test
  
  ###define Vasudha_test and Vasudha_training datasets
  Vasudha_test <- rbind(
    annotations_Vasudha[which(annotations_Vasudha$PERIOD=="pre"&annotations_Vasudha$Behaviour==BEH&annotations_Vasudha$REPLICATE==pre_colony_test),]
    ,
    annotations_Vasudha[which(annotations_Vasudha$PERIOD=="post"&annotations_Vasudha$Behaviour==BEH&annotations_Vasudha$REPLICATE==post_colony_test),]
  )
  
  Vasudha_training <- rbind(
    annotations_Vasudha[which(annotations_Vasudha$PERIOD=="pre"&annotations_Vasudha$Behaviour==BEH&annotations_Vasudha$REPLICATE==pre_colony_training),]
    ,
    annotations_Vasudha[which(annotations_Vasudha$PERIOD=="post"&annotations_Vasudha$Behaviour==BEH&annotations_Vasudha$REPLICATE==post_colony_training),]
  )
  
  ###define time_window_Vasudha_training and time_window_Vasudha_test
  time_window_Vasudha_training <- time_window_Vasudha[which((time_window_Vasudha$PERIOD=="pre"&time_window_Vasudha$REPLICATE==pre_colony_training)|(time_window_Vasudha$PERIOD=="post"&time_window_Vasudha$REPLICATE==post_colony_training)),]
  time_window_Vasudha_test     <- time_window_Vasudha[which((time_window_Vasudha$PERIOD=="pre"&time_window_Vasudha$REPLICATE==pre_colony_test)|(time_window_Vasudha$PERIOD=="post"&time_window_Vasudha$REPLICATE==post_colony_test)),]
  
  ###define all_training_ants training and non_training_ants - that is all ants in the colonies!
  all_training_ants <- full_ant_list_Vasudha[which( (grepl(post_colony_training,full_ant_list_Vasudha)) | (grepl(pre_colony_training,full_ant_list_Vasudha)))]
  non_training_ants <- full_ant_list_Vasudha[which( (grepl(post_colony_test,full_ant_list_Vasudha)) | (grepl(pre_colony_test,full_ant_list_Vasudha)))]
  
  ###add focal information
  receiver_focal <- with(Vasudha_training,full_Receiver_ID%in%all_training_ants)
  Vasudha_training[which(receiver_focal),"Focal"] <- Vasudha_training[which(receiver_focal),"Receiver"] 
  Vasudha_training[which(!receiver_focal),"Focal"] <- Vasudha_training[which(!receiver_focal),"Actor"] 
  
  return(list(Vasudha_training=Vasudha_training,Vasudha_test=Vasudha_test,all_training_ants=all_training_ants,non_training_ants=non_training_ants,time_window_Vasudha_training=time_window_Vasudha_training,time_window_Vasudha_test=time_window_Vasudha_test))
}


test_training_split_by_number <- function(annotations_all,time_window_all,seed){
  set.seed(seed)
  ##############################################################################
  ######### SPLIT ANNOTATIONS DATASET INTO TEST AND TRAINING  ##################
  ##############################################################################
  ### Be careful about this - in the previous you had randomly allocated behaviours to test or training throughout the period
  ### This means true Hits were wrongly identified as misses when looking at automatic interactions, because the time span of automatic interaction detection was unchanged
  ### Instead you need to define contiguous periods of time that contain half the events, for each colony/period
  print("Splitting manual annotations into test and training chunks...")
  ###first list nb of events for each behaviour, each replicate and each period
  nb_events <- aggregate ( Actor ~ Behaviour + PERIOD + REPLICATE, FUN=length, data=annotations_all)
  names(nb_events)[which(names(nb_events)=="Actor")] <- "Nb"
  ###narrow down to behaviour of interest only
  nb_events <- nb_events[which(nb_events$Behaviour==BEH),]
  
  ###initialise new test and training objects
  annotations_training <- NULL
  annotations_test     <- NULL
  time_window_training <- NULL
  time_window_test     <- NULL
  
  ###then loop over nb_events
  
  
  total_duplicated_events <- 0
  for (i in 1:nrow(nb_events)){
    ###subset annotations_all to period/replicate of interest
    PERIOD    <- nb_events[i,"PERIOD"]
    REPLICATE <- nb_events[i,"REPLICATE"]
    annotation_subset <- annotations_all[which(annotations_all$Behaviour==BEH & annotations_all$PERIOD==PERIOD & annotations_all$REPLICATE==REPLICATE),]
    
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to half the events for this period  
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to 3/4 of the events for this period  
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
    
    ### now we split annotation_subset in two chunks before and after time_limit,
    ### randomly allocate each time chunk to test or training dataset,
    ### and store the time windows for those time chunks
    if (runif(1,0,1)<0.5){
      
      time_limit <- min (c(as.numeric(annotation_subset[floor(proportion_training_events*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+proportion_training_events*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(proportion_training_events*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(proportion_training_events*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
      if (length(time_limit)==0){
        time_limit <- min (c(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+0.5*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(0.5*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
        
      }
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      
      
      
      annotations_training <- rbind( annotations_training , annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),]  )
      annotations_test     <- rbind( annotations_test     , annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),])
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                       time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      )
      )
      
      time_window_test     <- rbind( time_window_test     , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit + 1/FRAME_RATE,origin="1970-01-01", tz="GMT") - 1/FRAME_RATE, 
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
    }else{
      time_limit <- min (c(as.numeric(annotation_subset[floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+(1-proportion_training_events)*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor((1-proportion_training_events)*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
      if (length(time_limit)==0){
        time_limit <- min (c(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"]),as.numeric(annotation_subset[1+0.5*floor(nb_events[i,"Nb"]),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(0.5*nb_events[i,"Nb"]),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(0.5*nb_events[i,"Nb"]),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
        
      }
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      annotations_training <- rbind( annotations_training, annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),]  )
      annotations_test     <- rbind( annotations_test,annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),])
      
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit+1/FRAME_RATE,origin="1970-01-01", tz="GMT")- 1/FRAME_RATE,
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
      
      
      time_window_test     <- rbind( time_window_test , data.frame(PERIOD=PERIOD, 
                                                                   REPLICATE=REPLICATE,
                                                                   time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                   time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      ))
      
    }
    
  }
  
  ###finally, recalculate duration for all annotations_all objects and define new time columns express in seconds
  for (annotation_object in c("annotations_all","annotations_training","annotations_test")){
    annot <- get(annotation_object)
    #convert Zulu time to GMT
    annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
    # #transform zulu time in GMT
    # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    #assign time in sec to avoid issues on time management and matching
    annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
    annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
    
    assign(annotation_object,annot)
  }
  print(paste("Had to duplicate",total_duplicated_events,"events"))
  return(list(annotations_all=annotations_all,annotations_training=annotations_training,annotations_test=annotations_test,time_window_all=time_window_all,time_window_training=time_window_training,time_window_test=time_window_test))
}

test_training_split_by_time <- function(annotations_all,time_window_all){
  if (training_set=="both"){
    set.seed(7)
  }else{
    set.seed(9)
  }
  ##############################################################################
  ######### SPLIT ANNOTATIONS DATASET INTO TEST AND TRAINING  ##################
  ##############################################################################
  ### Be careful about this - in the previous you had randomly allocated behaviours to test or training throughout the period
  ### This means true Hits were wrongly identified as misses when looking at automatic interactions, because the time span of automatic interaction detection was unchanged
  ### Instead you need to define contiguous periods of time that contain half the events, for each colony/period
  print("Splitting manual annotations into test and training chunks...")
  
  ###initialise new test and training objects
  annotations_training <- NULL
  annotations_test     <- NULL
  time_window_training <- NULL
  time_window_test     <- NULL
  
  ###then loop over time_window_all
  total_duplicated_events <- 0
  for (i in 1:nrow(time_window_all)){
    ###subset annotations_all to period/replicate of interest
    PERIOD    <- time_window_all[i,"PERIOD"]
    REPLICATE <- time_window_all[i,"REPLICATE"]
    annotation_subset <- annotations_all[which(annotations_all$Behaviour==BEH & annotations_all$PERIOD==PERIOD & annotations_all$REPLICATE==REPLICATE),]
    
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to half the events for this period  
    ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to 3/4 of the events for this period  
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
    
    ### now we split annotation_subset in two chunks before and after time_limit,
    ### randomly allocate each time chunk to test or training dataset,
    ### and store the time windows for those time chunks
    time_chunk <- (as.numeric(time_window_all[i,"time_stop"])-as.numeric(time_window_all[i,"time_start"]))*proportion_training_events
    time_chunk <- round(time_chunk*FRAME_RATE)/FRAME_RATE
    if (runif(1,0,1)<0.5){ ###training = first_chunk
      time_limit <- as.numeric(time_window_all[i,"time_start"])+time_chunk
      
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      
      annotations_training <- rbind( annotations_training , annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),]  )
      annotations_test     <- rbind( annotations_test     , annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),])
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                       time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      )
      )
      
      time_window_test     <- rbind( time_window_test     , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit + 1/FRAME_RATE,origin="1970-01-01", tz="GMT") - 1/FRAME_RATE, 
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
    }else{
      time_limit <- as.numeric(time_window_all[i,"time_stop"])-time_chunk
      
      ###time_limit will work well if those two events don't overlap
      ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
      ###first list events to duplicate (if any)
      annotation_subset <- duplicate_annotations (annotation_subset,time_limit,total_duplicated_events)
      total_duplicated_events <- annotation_subset[["total_duplicated_events"]]
      annotation_subset       <- annotation_subset[["annotation_subset"]]
      
      annotations_training <- rbind( annotations_training, annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),]  )
      annotations_test     <- rbind( annotations_test,annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),])
      
      time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                       REPLICATE=REPLICATE,
                                                                       time_start = as.POSIXct(time_limit+1/FRAME_RATE,origin="1970-01-01", tz="GMT")- 1/FRAME_RATE,
                                                                       time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
      
      
      time_window_test     <- rbind( time_window_test , data.frame(PERIOD=PERIOD, 
                                                                   REPLICATE=REPLICATE,
                                                                   time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                   time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
      ))
      
    }
    
  }
  
  ###finally, recalculate duration for all annotations_all objects and define new time columns express in seconds
  for (annotation_object in c("annotations_all","annotations_training","annotations_test")){
    annot <- get(annotation_object)
    #convert Zulu time to GMT
    annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
    # #transform zulu time in GMT
    # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
    #assign time in sec to avoid issues on time management and matching
    annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
    annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
    
    assign(annotation_object,annot)
  }
  return(list(annotations_all=annotations_all,annotations_training=annotations_training,annotations_test=annotations_test,time_window_all=time_window_all,time_window_training=time_window_training,time_window_test=time_window_test))
}

###############################################################################
###### PARAMETERS #############################################################
###############################################################################
###body length information
all_body_lengths <-read.table(BODYLENGTH_FILE,header=T,stringsAsFactors = F,sep=",")

#plotting limits used for the coordinates plotting
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

#general
N_DECIMALS                  <- 3 ## when assigning time in seconds, the number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH                 <- TRUE  ## fuzzy matching between data frames in collisions detection

#trajectories cutting gap, relevant for fmQueryComputeAntTrajectories
max_gap                     <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)

###NATH_FLAG: don't just use max and min, give some extra 
###NATH_FLAG: these parameters should be adjusted for each box as ant length may depend on focus /  camera distance
# AntDistanceSmallerThan      <- 300 #for higher accuracy, recalculate it from: max(interaction_MANUAL$straightline_dist_px,na.rm = T)
# AntDistanceGreaterThan      <- 70 #for higher accuracy, recalculate it from: min(interaction_MANUAL$straightline_dist_px,na.rm = T)
# ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
# minimumGap                  <- fmSecond(0.5) ## THIS OPTION DOES NOT WORK IN INTERACTIONS SO DISABLED! for a given pair of interacting ants, when interaction is interrupted by more than minimumGap, interaction will check whether ants have moved since - and if so, will create new interaction

DISAGREEMENT_THRESH <- 0.5
###Fixed parameter

###############################################################################
###READ WHOLE ANNOTATIONS DATASET #############################################
###############################################################################
print("Loading manual annotations...")
#the current annotation file FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv underwent cross-validation by Adriano
annotations_Vasudha   <- modify_columns(read.csv(paste(WORKDIR,"/Data","/R3SP_R9SP_All_data_FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv",sep = ""), sep = ","),BEH)
annotations_Adriano   <- modify_columns(read.csv(paste(WORKDIR,"Data","Grooming_Classifier_CrossVal_ANNOTATIONS.csv",sep = "/"), sep = ","),BEH)
annotations_Charlotte <- modify_columns(read.csv(paste(WORKDIR,"Data","Grooming_Classifier_CrossVal_Charlotte2023_ANNOTATIONS.csv",sep = "/"), sep = ","),BEH)
# add Charlotte annotations to Adriano table
annotations_Adriano <- annotations_Adriano[which(names(annotations_Adriano)%in%names(annotations_Charlotte))]
annotations_Adriano <- rbind(annotations_Adriano,annotations_Charlotte[names(annotations_Adriano)])

###read or extract metadata - Vasudha
focal_list_Vasudha  <- read.table(paste(WORKDIR,"/Data","/Exposed_nurses_R3SP_R9SP.txt",sep = ""), sep = ",",header=T,stringsAsFactors = F)
time_window_Vasudha        <- merge(aggregate(T_start_UNIX ~ PERIOD + REPLICATE, FUN=min, data=annotations_Vasudha),aggregate(T_stop_UNIX  ~ PERIOD + REPLICATE, FUN=max, data=annotations_Vasudha))
names(time_window_Vasudha) <- c("PERIOD","REPLICATE","time_start","time_stop")

###read metadata - Adriano
metadata_info_Adriano    <- read.csv(paste(WORKDIR,"Data","Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU.csv",sep = "/"), sep = ",")
focal_list_Adriano       <- unique(with(metadata_info_Adriano[which(!is.na(metadata_info_Adriano$antID)),]   ,paste(REP_treat,antID   ,sep="_")))
metadata_info_Charlotte  <- read.csv(paste(WORKDIR,"Data","Grooming_Classifier_CrossVal_Charlotte2023_RETURN_EXP_TIME_ZULU.csv",sep = "/"), sep = ",")
focal_list_Charlotte     <- unique(with(metadata_info_Charlotte[which(!is.na(metadata_info_Charlotte$antID)),]   ,paste(REP_treat,antID   ,sep="_")))
# add Charlotte metadata to Adriano metadata
focal_list_Adriano       <- unique(c(focal_list_Adriano,focal_list_Charlotte))

time_window_Adriano        <- data.frame(
  PERIOD = "post"
  ,REPLICATE = metadata_info_Adriano$REP_treat
  ,time_start = metadata_info_Adriano$ReturnExposed_time
)
time_window_Adriano$time_start <- as.POSIXct(time_window_Adriano$time_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
###add two minutes to start
time_window_Adriano$time_start <- time_window_Adriano$time_start + 2*60
time_window_Adriano$time_stop  <- time_window_Adriano$time_start + 30*60


if (training_set=="Adriano"){
   annotations_all_training  <- annotations_Adriano
  time_window_all_training   <- time_window_Adriano
  FOCAL_training             <- T
  all_training_ants <- focal_list_Adriano
  
  
  annotations_non_training   <- annotations_Vasudha
  time_window_non_training   <- time_window_Vasudha
  FOCAL_non_training         <- F
  non_training_ants          <- NULL
  
}else if (training_set=="Vasudha"){
  annotations_all_training   <- annotations_Vasudha
  time_window_all_training   <- time_window_Vasudha
  FOCAL_training             <- F
  all_training_ants          <- NULL
  
  annotations_non_training   <- annotations_Adriano
  time_window_non_training   <- time_window_Adriano
  FOCAL_non_training         <- T
  non_training_ants <- unique(with(annotations_Adriano[which(!is.na(annotations_Adriano$Focal)),]   ,paste(REPLICATE,Focal   ,sep="_")))
  
}else if (training_set=="both"){
  
  ##subset a random number of ants in annotations_Vasudha
  split_annotations_Vasudha    <- Vasudha_split(annotations_Vasudha,focal_list_Vasudha,time_window_Vasudha,annotations_Adriano,BEH)
  
  ###assign non-training dataset
  annotations_non_training   <- split_annotations_Vasudha[["Vasudha_test"]]
  time_window_non_training   <- split_annotations_Vasudha[["time_window_Vasudha_test"]]
  FOCAL_non_training         <- T
  non_training_ants          <- split_annotations_Vasudha[["non_training_ants"]]
  
  
  ##DO NOT add Charlotte focals to non_training_ants
  # non_training_ants <- c(non_training_ants, focal_list_Charlotte)
  
  ##DO NOT merge Charlotte and Vasudha non_training
  # common_names             <- names(annotations_non_training)[which(names(annotations_non_training)%in%names(annotations_Charlotte))]
  # annotations_non_training <- rbind(annotations_Charlotte[common_names],annotations_non_training[common_names])
  # time_window_non_training <- rbind(time_window_Adriano,time_window_non_training)
  
  
  ###define Vasudha training dataset
  annotations_Vasudha_all_training <- split_annotations_Vasudha[["Vasudha_training"]]
  time_window_Vasudha_all_training <- split_annotations_Vasudha[["time_window_Vasudha_training"]]
  FOCAL_training                   <- T
  all_training_ants                <- split_annotations_Vasudha[["all_training_ants"]]
  
  ##add Adriano focals to all_training_ant_list
  all_training_ants <- c(all_training_ants, unique(with(annotations_Adriano[which(!is.na(annotations_Adriano$Focal)),]   ,paste(REPLICATE,Focal   ,sep="_"))))
  
  ##merge Adriano and Vasudha training
  common_names <- names(annotations_Vasudha_all_training)[which(names(annotations_Vasudha_all_training)%in%names(annotations_Adriano))]
  annotations_all_training <- rbind(annotations_Adriano[common_names],annotations_Vasudha_all_training[common_names])
  time_window_all_training <- rbind(time_window_Adriano,time_window_Vasudha_all_training)
}

###make sure annotations don't overrun 30 minute window
for (suffix in c("all_training","non_training")){
  annot     <- get(paste("annotations",suffix,sep="_"))
  time_wind <- get(paste("time_window",suffix,sep="_"))
  
  for (i in 1:nrow(annot)){
    annot[i,"T_stop_UNIX"] <- min(annot[i,"T_stop_UNIX"],time_wind[which(time_wind$REPLICATE==annot[i,"REPLICATE"]&time_wind$PERIOD==annot[i,"PERIOD"]),"time_stop"]-1/FRAME_RATE)
    annot[i,"T_start_UNIX"] <- max(annot[i,"T_start_UNIX"],time_wind[which(time_wind$REPLICATE==annot[i,"REPLICATE"]&time_wind$PERIOD==annot[i,"PERIOD"]),"time_start"]+1/FRAME_RATE)
  }
  annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
  annot               <- annot[which(annot$duration>0),] 
  
  ###finally, subset annot for the behaviour of interest
  annot <- annot[which(annot$Behaviour==BEH),]
  
  assign(paste("annotations",suffix,sep="_"),annot)
}

for (annotation_object in c("annotations_non_training")){
  annot <- get(annotation_object)
  #convert Zulu time to GMT
  annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
  # #transform zulu time in GMT
  # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  #assign time in sec to avoid issues on time management and matching
  annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
  annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
  
  assign(annotation_object,annot)
}

###split into test and training sets
random_split <- test_training_split_by_number(annotations_all_training,time_window_all_training,seed=20)
# random_split <- test_training_split_by_time(annotations_all_training,time_window_all_training)

annotations_all_training      <- random_split[["annotations_all"]]
annotations_training_training <- random_split[["annotations_training"]]
annotations_training_test     <- random_split[["annotations_test"]]
time_window_all_training      <- random_split[["time_window_all"]]
time_window_training_training <- random_split[["time_window_training"]]
time_window_training_test     <- random_split[["time_window_test"]]



###############################################################################
###### EXTRACT CHOSEN PARAMETERS FROM CHOSEN CLASSIFIER########################
###############################################################################

###define name of general output table containing quality scores
chosen_file_name <- file.path(MachineLearningOutcome_DIR,paste("quality_scores_",measure,"_priority",priority,"_CHOSEN.txt",sep=""))
chosen <- read.table(chosen_file_name,header=T,stringsAsFactors = F)

# #####Arguments to loop over - to comment out when running the loop
subDir                      <- paste0("Loop_ID_",chosen[,"Loop_ID"])
CAPSULE_FILE                <- chosen[,"CAPSULE_FILE"]
DT_dist_THRESHOLD_BL        <- chosen[,"DT_dist_THRESHOLD_BL"]
MAX_INTERACTION_GAP         <- chosen[,"MAX_INTERACTION_GAP"]
DISAGREEMENT_THRESH         <- chosen[,"DISAGREEMENT_THRESH"]
trim_length_sec             <- chosen[,"trim_length_sec"]
DT_frame_THRESHOLD          <- chosen[,"DT_frame_THRESHOLD"]
beta                        <- chosen[,"beta"]
###Load BN_list
BN_list <- dget (  file.path(MachineLearningOutcome_DIR,subDir,"BN_object_list.dat") )
###Load_classifier
classifier        <- list(readRDS (file.path(MachineLearningOutcome_DIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep=""))))
names(classifier) <- chosen[,"classifier"]

#############################################################################################
##################### EXTRACT VARIABLES FOR MANUAL AND AUTOMATIC INTERACTIONS ###############
#############################################################################################
###for these particular parameters, evaluate how successful the interaction detection parameters are at detecting candidate frames
###extract interactions and manual annotations for the non-training dataset
non_training           <- extraction_loop(chunk="non_training"
                                          ,extract_movement_variables=T
                                          ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
                                          ,all_body_lengths=all_body_lengths
                                          ,focal=FOCAL_non_training
                                          ,focal_list = non_training_ants
)

###extract interactions and manual annotations for the training subset of the training dataset
training_training      <- extraction_loop(chunk="training_training"
                            ,extract_movement_variables=T
                            ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
                            ,all_body_lengths=all_body_lengths
                            ,focal=FOCAL_training
                            ,focal_list = all_training_ants
)

###extract interactions and manual annotations for the test subset of the training dataset
training_test           <- extraction_loop(chunk="training_test"
                                          ,extract_movement_variables=T
                                          ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
                                          ,all_body_lengths=all_body_lengths
                                          ,focal=FOCAL_training
                                          ,focal_list = all_training_ants
)




# ###repeat training_test, but without subsetting focals
# training_test_all_ants <- extraction_loop(chunk="training_test"
#                                            ,extract_movement_variables=T
#                                            ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
#                                            ,all_body_lengths=all_body_lengths
#                                            ,focal=F
#                                            ,focal_list = NULL
# )
# 

###now, make predictions for all these objects
for (object in c("non_training","training_training","training_test")){#,"training_training","training_test","training_test_all_ants",random_split <- test_training_split_by_time(annotations_all_training,time_window_all_training)
  dataset <- get(object)
  
  candidate_groomings <- dataset[["summary_AUTO"]]
  candidate_groomings$predicted_Hit <- NULL
  #~~~~~PREDICT GROONING USING CLASSIFIED
  print("Predicting grooming...")
  prediction_start <- Sys.time()
  candidate_groomings["predicted_Hit"] <- predict_class  (summary_AUTO =    candidate_groomings
                                                          ,BN_list      =  BN_list
                                                          ,classifier   = classifier
  )
  print("Grooming predicted.")
  
  # update object
  dataset[["summary_AUTO"]] <- candidate_groomings
  assign(object,dataset)
}

# ###BIG QUESTION 1 - do predictions depend on distribution of false/true events in the dataset?
# large_prediction_set <- training_test_all_ants$summary_AUTO
# small_prediction_set <- training_test$summary_AUTO
# ###For this, compare training_test and training_test_all_ants
# nrow(large_prediction_set) ###training_test_all_ants has 17101 interactions
# nrow(small_prediction_set)  ###training_test has 714 interactions
# ###create a unique interaction identity so we can cross-compare predictions on the two sets
# large_prediction_set <- within(large_prediction_set, unique_interaction_id <- paste(REPLICATE,PERIOD,ant1,ant2,frame_start,sep="_"))
# small_prediction_set <- within(small_prediction_set, unique_interaction_id <- paste(REPLICATE,PERIOD,ant1,ant2,frame_start,sep="_"))
# ###check that small_prediction_set is an exact subset of large_prediction_set using unique_interaction_id
# all(small_prediction_set$unique_interaction_id%in%large_prediction_set$unique_interaction_id)
# ###now, subset large_prediction_set to match exactly the interactions contained in the small_prediction_set
# large_prediction_set_SUBSET <- large_prediction_set[which(large_prediction_set$unique_interaction_id%in%small_prediction_set$unique_interaction_id),]
# nrow(large_prediction_set_SUBSET)
# ###reorder large_prediction_set_SUBSET so it is the same order as small_prediction_set
# large_prediction_set_SUBSET <- large_prediction_set_SUBSET[match(small_prediction_set$unique_interaction_id,large_prediction_set_SUBSET$unique_interaction_id),]
# all(large_prediction_set_SUBSET$unique_interaction_id==small_prediction_set$unique_interaction_id)
# ###and finally, checks if hits are the same!
# all(na.omit(large_prediction_set_SUBSET$predicted_Hit==small_prediction_set$predicted_Hit))
# head(small_prediction_set[which(is.na(small_prediction_set$predicted_Hit)),])
  

###BIG QUESTION 2 - how well do the predictions do?
###get overall match for all these objects
quality_scores_REP_PER  <- NULL
for (object in c("non_training","training_training","training_test")){#
  dataset <- get(object)
  candidate_groomings <- dataset[["summary_AUTO"]]
  summary_AUTO                 <- candidate_groomings[which(candidate_groomings$predicted_Hit==1),]
  summary_MANUAL               <- dataset[["summary_MANUAL"]]
  list_IF_Frames               <- dataset[["list_IF_Frames"]]
  list_replicate_full_ant_list <- dataset[["list_replicate_full_ant_list"]]
  list_replicate_focal_list    <- dataset[["list_replicate_focal_list"]]

  ###loop over all data / per replicate / per period
  for (coverage in c("all","per_rep")){
    if (coverage=="all"){
      if (object=="non_training"){
        dataset_name <- "Vasudha"
      }else{
        dataset_name <- "Vasudha_Charlotte_and_Adriano"
      }
      
          aut_man_agreement <- auto_manual_agreement (    summary_AUTO   = summary_AUTO
                                                        , summary_MANUAL = summary_MANUAL
                                                        , list_IF_Frames = list_IF_Frames
                                                        , list_replicate_full_ant_list = list_replicate_full_ant_list
                                                        , list_replicate_focal_list = list_replicate_focal_list
        )

      ###overall test
      quality_scores_overall <-round(quality_scores(aut_man_agreement[["true_false_positive_negatives"]],beta),digits=3 )
      quality_scores_REP_PER <- rbind(quality_scores_REP_PER,data.frame(what=object,replicate="all",dataset=dataset_name,period="all",as.list(quality_scores_overall),as.list(aut_man_agreement[["true_false_positive_negatives"]])))
      
    }else{
      for (REP in unique(candidate_groomings$REPLICATE)) {
        if (REP%in%c("R3SP","R9SP")){
          dataset_name <- "Vasudha"
          
        }else {
          dataset_name <- "Charlotte_plus_Adriano"
          
        }
        
        
        for (PER in unique(candidate_groomings[which(candidate_groomings$REPLICATE==REP),"PERIOD"])){
          if (nrow(summary_AUTO[which(summary_AUTO$REPLICATE==REP&summary_AUTO$PERIOD==PER),])<1) {
            #temporary solution for REPs without candidate grooming
            quality_scores_REP_PER <- rbind(quality_scores_REP_PER,data.frame(what=object,replicate=REP,dataset=dataset_name,period=PER,CSI=NA,Fbeta=NA,precision=NA,sensitivity=NA,true_negatives=NA,true_positives=NA,false_negatives=NA,false_positives=NA))
          }else{
             aut_man_agreement_REP_PER <- auto_manual_agreement (  summary_AUTO   = summary_AUTO  [which(summary_AUTO$REPLICATE==REP&summary_AUTO$PERIOD==PER),]
                                                                , summary_MANUAL = summary_MANUAL[which(summary_MANUAL$REPLICATE==REP&summary_MANUAL$PERIOD==PER),]
                                                                , list_IF_Frames = list_IF_Frames
                                                                , list_replicate_full_ant_list = list_replicate_full_ant_list
                                                                ,  list_replicate_focal_list = list_replicate_focal_list
            )
 
             qual_scores <-round(quality_scores(aut_man_agreement_REP_PER[["true_false_positive_negatives"]],beta),digits=3 )
            quality_scores_REP_PER <- rbind(quality_scores_REP_PER,data.frame(what=object,replicate=REP,dataset=dataset_name,period=PER,as.list(qual_scores),as.list(aut_man_agreement_REP_PER[["true_false_positive_negatives"]])))
          }
          
        }
      }
    }
    
    
  }
  
  
}

###explore how variable those predictions are
quality_scores_REP_PER <- within(  quality_scores_REP_PER,  treatment <- substr(replicate,nchar(replicate)-1,nchar(replicate)))
quality_scores_REP_PER$treatment <- factor(quality_scores_REP_PER$treatment )
# quality_scores_REP_PER$dataset   <- "Adriano"
# quality_scores_REP_PER[which(quality_scores_REP_PER$replicate%in%c("R3SP","R9SP")|quality_scores_REP_PER$what=="non_training"),"dataset"] <- "Vasudha"
##overall score, pooled per treatment
quality_scores_summary <- aggregate (cbind(true_negatives,true_positives,false_negatives,false_positives) ~ what  + period + treatment , FUN=sum, data=quality_scores_REP_PER)
quality_scores_summary <- cbind(quality_scores_summary,t(apply(quality_scores_summary[c("true_negatives","true_positives","false_negatives","false_positives")],1,FUN=quality_scores,beta=beta)))

quality_scores_summary <- quality_scores_summary[order(quality_scores_summary$what,quality_scores_summary$period,quality_scores_summary$treatment),]

###comparisons
print(aggregate(cbind(precision,sensitivity)~what + dataset + period + treatment,function(x)cbind(mean(x)),data=quality_scores_REP_PER[which(quality_scores_REP_PER$what!="training_training"&quality_scores_REP_PER$period=="post"&quality_scores_REP_PER$dataset=="Charlotte_plus_Adriano"),]))

model_sensitivity <- lm(    sensitivity^3 ~ treatment, data=quality_scores_REP_PER[which(quality_scores_REP_PER$what!="training_training"&quality_scores_REP_PER$period=="post"&quality_scores_REP_PER$dataset=="Charlotte_plus_Adriano"),])
shapiro.test(residuals(model_sensitivity))    
anova(model_sensitivity)

model_precision <- lm(    precision^2~ treatment, data=quality_scores_REP_PER[which(quality_scores_REP_PER$what!="training_training"&quality_scores_REP_PER$period=="post"&quality_scores_REP_PER$dataset=="Charlotte_plus_Adriano"),])
shapiro.test(residuals(model_precision))    
anova(model_precision)

write.table(quality_scores_summary,file=file.path(MachineLearningOutcome_DIR, paste("quality_scores_",measure,"_",priority,"_summary.txt",sep="")),col.names=T,row.names=F,append = F,quote=F)
