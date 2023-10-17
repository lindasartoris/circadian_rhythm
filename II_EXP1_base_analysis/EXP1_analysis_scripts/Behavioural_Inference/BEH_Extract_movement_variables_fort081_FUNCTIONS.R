extraction_loop <- function(chunk,extract_movement_variables=T, selected_variables = NULL,all_body_lengths=NULL,focal="all",focal_list=NULL){
  ### initialise variables
  interaction_MANUAL      <- NULL
  summary_MANUAL          <- NULL
  interaction_AUTO        <- NULL
  summary_AUTO            <- NULL
  list_IF_Frames          <- list()
  list_replicate_full_ant_list <- list()
  list_replicate_focal_list <- list()
  ### get appropriate time_window and annotations objects
  time_window_ori <- get(paste("time_window_",chunk,sep=""))
  annotations_ori <- get(paste("annotations_",chunk,sep=""))
  
  ###loop over unique combinations of REPLICATE and PERIOD
  for (REPLICATE in unique(time_window_ori$REPLICATE)) 
  {
    ###############################################################################
    ###### OPEN EXPERIMENT INFORMATION ############################################
    ###############################################################################
    
    ## locate the ant info file for REPLICATE
    MyrmidonCapsuleFiles <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
    MyrmidonCapsuleFiles <- MyrmidonCapsuleFiles[which(grepl("myrmidon",MyrmidonCapsuleFiles))]
    
    if (length(MyrmidonCapsuleFiles)>0){
      MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl(CAPSULE_FILE,MyrmidonCapsuleFiles)]
      # if (length(MyrmidonCapsuleFiles)>0){
      #   MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl(gsub("Capsule","Cap",CAPSULE_FILE),MyrmidonCapsuleFiles)]
      # }
    }else{
      MyrmidonCapsuleFiles <- list.files(path=file.path(DATADIR,gsub( "R", "REP"  ,     substr(REPLICATE,1,nchar(REPLICATE)-2))), pattern=REPLICATE, full.names=T)
      MyrmidonCapsuleFiles <- MyrmidonCapsuleFiles[which(grepl("myrmidon",MyrmidonCapsuleFiles))]
      MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl(CAPSULE_FILE,MyrmidonCapsuleFiles)]
      # if (length(MyrmidonCapsuleFiles)>0){
      #   MyrmidonCapsuleFile <- MyrmidonCapsuleFiles[grepl(gsub("Capsule","Cap",CAPSULE_FILE),MyrmidonCapsuleFiles)]
      # }
    }
    
     print(MyrmidonCapsuleFile)
    e <- fmExperimentOpen(MyrmidonCapsuleFile)
    replicate_full_ant_list <- unlist(lapply(e$ants,function(x)x$ID))
    
    if (is.null(focal_list)){
      replicate_focal_list <- NULL
    }else{
      replicate_focal_list <- sort(as.numeric(gsub(paste(REPLICATE,"_",sep=""),"",focal_list[which(grepl(REPLICATE,focal_list))])))
    }
    
    body_lengths <- get_body_lengths(e,all_body_lengths)
    
    ################################################################################
    ########### START PERIOD LOOP ##################################################
    ################################################################################
    
    for (PERIOD in unique(time_window_ori[which(time_window_ori$REPLICATE==REPLICATE),       "PERIOD"]))
    {
      print(paste(REPLICATE,PERIOD))
      time_window <- time_window_ori[which(time_window_ori$REPLICATE==REPLICATE&time_window_ori$PERIOD==PERIOD),      ]
      annotations <- annotations_ori[which(annotations_ori$REPLICATE==REPLICATE&annotations_ori$PERIOD==PERIOD),      ]
      
      ## set experiment time window 
      time_start <- fmTimeCreate(time_window[which(time_window$PERIOD==PERIOD & time_window$REPLICATE==REPLICATE),"time_start"])
      time_stop  <- fmTimeCreate(time_window[which(time_window$PERIOD==PERIOD & time_window$REPLICATE==REPLICATE),"time_stop"])
      
      # print(time_start)
      # print(time_stop)
      ###############################################################################
      ###### IDENTIFY FRAMES ########################################################
      ###############################################################################
      
      #Because of various issues raised, including the ones reported here https://github.com/formicidae-tracker/myrmidon/issues/240 ,
      #the analysis is performed not using UNIX_time but using frames. Every new queried file will then show a part including FRAME assignment.
      #Frames are then used for any trajectory cutting later on in the sub-scripts
      IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_stop,showProgress = FALSE)
      IF_frames           <- IdentifyFrames$frames
      # Assign a frame to each time since start and use it as baseline for all matching and computation
      IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))
      
      # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
      IF_frames$time_sec <- round(as.numeric(IF_frames$time),N_DECIMALS)
      
      ###NATH_FLAG
      # deduce FRAME_RATE from IF_frames
      # FRAME_RATE <- rev(c(1:30))[match.closest(x=median(diff(IF_frames$time_sec)),table=rev(1/c(1:30)))]
      
      #assign frame numbering to annotations 
      ###NATH_FLAG: why do you not create a new object that subsets the annotations for the replicate and period you want?
      ### here you have an object with lots of irrelevant rows
      
      
      annotations$frame_start <- match.closest(x = annotations$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
      annotations$frame_stop  <- match.closest(x = annotations$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
      
      # Creating a new zeroed-time since the start of the exp by  summing the cumulated differences between each pair of consecutive frames (to account for the time mismatch)
      IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
      IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)
      
      ###############################################################################
      ###### READING TRAJECTORIES ###################################################
      ###############################################################################
      #COMPUTE TRAJECTORIES  
      positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE,showProgress = FALSE) #set true to obtain the zone of the ant
      #add body lengths to trajectories
      positions <- add_body_length_to_traj(positions,body_lengths)
      
      positions$trajectories_summary$frame_num <- NA
      #assign starting frame number
      positions$trajectories_summary["frame_num"] <- lapply(positions$trajectories_summary["start"], function(x) IF_frames$frame_num[match(x, IF_frames$time)])
      
      ## immediately after computing your positions object:
      positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
      names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
      
      trajectories_summary <- positions$trajectories_summary
      
      ## Add ant_x names and times to the positions to convert from "FRAME since the start of the experiment", to FRAMES
      for (A in positions$trajectories_summary$antID)
      {
        AntID                         <- paste("ant_",A,sep="") ###NATH_FLAG: why not loop directly on antID_str which you just created?
        First_Obs_Frame               <- positions$trajectories_summary$frame_num [which(positions$trajectories_summary$antID_str==AntID)]
        IF_frames$new_zero_diff  <- IF_frames$cum_diff - IF_frames[IF_frames$frame_num==First_Obs_Frame,"cum_diff"] #subtracting the $frames zeroed-time  corresponding to the $start time from the zeroed-time column itself (New-Zeroed-time)
        # print(paste("Adding first obs FRAME", First_Obs_Frame, "to the time-zeroed trajectory of ant", AntID))
        #assign corresponding frame N when the New-Zeroed-time and $time correspond, closest.match 0.05 (well inside the 0.125 frame length in sec)
        positions$trajectories[[AntID]]$frame <- match.closest(x = positions$trajectories[[AntID]]$time, table = IF_frames$new_zero_diff, tolerance = 0.05)
        IF_frames$new_zero_diff <- NA
      }
      AntID <- NULL; First_Obs_Frame <- NULL
      # the older version used UNIX_time
      # ## Add ant_x names and times to the positions to convert from time since the start of the experiment, to UNIX time
      # for (A in positions$trajectories_summary$antID)
      #   {
      #   AntID <- paste("ant_",A,sep="") 
      #   First_Obs_Time <- as.POSIXct(positions$trajectories_summary$start [which(positions$trajectories_summary$antID_str==AntID)], tz="GMT",origin = "1970-01-01 00:00:00") ## find the first time after the user defined time_start_ISO that this ant was seen
      #   print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
      #   positions$trajectories[[AntID]] $UNIX_time <- as.POSIXct(positions$trajectories[[AntID]]$time, tz="GMT",origin = "1970-01-01 00:00:00")  + First_Obs_Time ##convert back to UNIX time  
      # }
      
      ######################################################################################
      ###### EXTRACT ANT TRAJECTORIES FROM MANUALLY-ANNOTATED DATA AND CALC PARAMETERS #####
      ######################################################################################
      ####First extract variables for manual interactions
      ## subset all hand-labelled bahavs for this behaviour type in this REPLICATE AND PERIOD
      # if (extract_movement_variables){
      #   extracted_variables_MANUAL <- extract_from_object(annotations
      #                                                     ,IF_frames
      #                                                     ,positions
      #                                                     ,BEH
      #                                                     ,REPLICATE 
      #                                                     ,PERIOD
      #                                                     ,selected_variables
      #   )
      #   
      #   summary_MANUAL     <- rbind(summary_MANUAL    ,extracted_variables_MANUAL[["summary_variables"]])
      # }else{
        ##ADD ant1 ant2
        annot_per_rep      <- annotations[which(annotations$Behaviour==BEH&annotations$REPLICATE==REPLICATE&annotations$PERIOD==PERIOD),]
        annot_per_rep$ant1 <- annot_per_rep$Actor; annot_per_rep$ant2 <- annot_per_rep$Receiver;annot_per_rep$pair <- apply(annot_per_rep[,c("ant1","ant2")],1,function(x){paste(sort(as.numeric(x)),collapse = "_") })
        summary_MANUAL     <- rbind(summary_MANUAL,annot_per_rep)   
      # }
      # interaction_MANUAL <- rbind(interaction_MANUAL,extracted_variables_MANUAL[["stored_trajectories"]])
      
      #########################################################################################
      ###### READING AUTOMATIC INTERACTIONS ###################################################
      #########################################################################################
      #Get info on the capules shapes and use the relevant ones in the ant interaction query
      capsules  <- e$antShapeTypeNames
      names(capsules) <- as.character( 1:length(capsules))
      
      ALL_CAPS_MATCHERS <- list()
      for (interaction_of_interest in interactions_of_interest){
        caps_matcher <- c()
        for (caps in interaction_of_interest[c(1:2)]){
          caps_matcher <- c(caps_matcher,as.numeric(names(capsules)[[which(capsules==caps)]]) )
        }
        ALL_CAPS_MATCHERS <- c(ALL_CAPS_MATCHERS,list(caps_matcher))
      }
      # head_id <- as.numeric(names(capsules)[[which(capsules=="head")]])
      # body_id <- as.numeric(names(capsules)[[which(capsules=="body")]])
      
      if (focal==T){
        interacts_AUTO_REP_PER <- interaction_detection (e=e
                                                         ,start=time_start
                                                         ,end=time_stop
                                                         ,max_time_gap = MAX_INTERACTION_GAP
                                                         ,max_distance_moved = 2*mean(body_lengths$body_length,na.rm=T)
                                                         ,capsule_matcher=ALL_CAPS_MATCHERS
                                                         ,desired_ants_OR = replicate_focal_list
                                                         ,IF_frames=IF_frames
        )
        
      }else{
        interacts_AUTO_REP_PER <- interaction_detection (e=e
                                                         ,start=time_start
                                                         ,end=time_stop
                                                         ,max_time_gap = MAX_INTERACTION_GAP
                                                         ,max_distance_moved = 2*mean(body_lengths$body_length,na.rm=T)
                                                         ,capsule_matcher=ALL_CAPS_MATCHERS
                                                         ,IF_frames=IF_frames
        )
      }
      
      # matcherCapType <- fmMatcherInteractionType(body_id,head_id)
      # matcherCapTypeAntDists <- fmMatcherAnd(list(fmMatcherInteractionType(body_id,head_id)
      #                                             # ,
      #                                             # fmMatcherAntDistanceSmallerThan(AntDistanceSmallerThan)
      #                                             # ,
      #                                             # fmMatcherAntDistanceGreaterThan(AntDistanceGreaterThan)
      #                                             # ,
      #                                             # fmMatcherAntDisplacement(mean(body_lengths$body_length), minimumGap) #check every minimumGap seconds if ant has displaced more than mean body length - if so, will create new interaction
      # ))
      # # AN EXTRA 2 THAT CAN BE USED ARE #fmMatcherAntAngleSmallerThan(), fmMatcherAntAngleGreaterThan()
      # 
      # ###NATH_FLAG: stick to object position to extract trajectories from automatic data
      # #note: adding AntDistances seems to reduce the false positives rate but has no effect on false negatives
      # interacts_AUTO_REP_PER <- fmQueryComputeAntInteractions(e,
      #                                                         start=time_start,
      #                                                         end=time_stop,
      #                                                         maximumGap = fmSecond(MAX_INTERACTION_GAP) , ## WHEN A PAIR DISENGAGE, HOW LONG IS THE INTERVAL?
      #                                                         reportFullTrajectories = F
      #                                                         ,showProgress = FALSE
      #                                                         ,
      #                                                         matcher = matcherCapTypeAntDists
      # )
      # 
      
      
      
      ## Add ant_x names and times to the interacts_AUTO_REP_PER to convert from FRAME since the start of the experiment, to FRAMES
      ##creates a ID string for each ant in $interactions
      interacts_AUTO_REP_PER$ant1ID_str            <- paste("ant_",interacts_AUTO_REP_PER$ant1,sep="")
      interacts_AUTO_REP_PER$ant2ID_str            <- paste("ant_",interacts_AUTO_REP_PER$ant2,sep="")
      # Assign interaction pair
      interacts_AUTO_REP_PER$pair <- paste(interacts_AUTO_REP_PER$ant1, interacts_AUTO_REP_PER$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
      
      ##convert times to different formats
      interacts_AUTO_REP_PER$T_start_UNIX <- as.POSIXct(interacts_AUTO_REP_PER$start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
      interacts_AUTO_REP_PER$T_stop_UNIX  <- as.POSIXct(interacts_AUTO_REP_PER$end,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
      #assign time in sec to avoid issues on time management and matching
      interacts_AUTO_REP_PER$T_start_sec <- round(as.numeric(interacts_AUTO_REP_PER$T_start_UNIX),N_DECIMALS)
      interacts_AUTO_REP_PER$T_stop_sec <- round(as.numeric(interacts_AUTO_REP_PER$T_stop_UNIX),N_DECIMALS)
      
      #assign start and end frame number
      interacts_AUTO_REP_PER$frame_start <- match.closest(x = interacts_AUTO_REP_PER$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
      interacts_AUTO_REP_PER$frame_stop  <- match.closest(x = interacts_AUTO_REP_PER$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
      
      #calc duration (including start frame)
      interacts_AUTO_REP_PER$duration <- interacts_AUTO_REP_PER$T_stop_sec - interacts_AUTO_REP_PER$T_start_sec + 1/FRAME_RATE
      #hist(interacts_AUTO_REP_PER$interactions$Duration,breaks = 30)
      
      ######################################################################################
      ###### CALC PARAMETERS FOR AUTOMATIC INTERACTIONS ####################################
      ######################################################################################
      if (extract_movement_variables){
        
        extracted_variables_AUTO <- extract_from_object(interacts_AUTO_REP_PER
                                                        ,IF_frames
                                                        ,positions
                                                        ,BEH
                                                        ,REPLICATE 
                                                        ,PERIOD
                                                        ,selected_variables
        )
        
        if(!is.null(extracted_variables_AUTO[["summary_variables"]])){
          summary_AUTO     <- rbind(summary_AUTO    , extracted_variables_AUTO[["summary_variables"]])
        }else{
          inter_list <- data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,Behaviour=BEH,interacts_AUTO_REP_PER,stringsAsFactors = F)
          inter_list[names(summary_AUTO)[which(!names(summary_AUTO)%in%names(inter_list))]] <- NA
          summary_AUTO     <- rbind(summary_AUTO    ,  inter_list[names(summary_AUTO)])
          
        }
        # interaction_AUTO <- rbind(interaction_AUTO, extracted_variables_AUTO[["stored_trajectories"]])
      }else{
        summary_AUTO     <- rbind(summary_AUTO    , data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,Behaviour=BEH,interacts_AUTO_REP_PER,stringsAsFactors = F))
       }
      
      
      ###NATH MARKER
      ###############################################################################
      ###### COLLISIONS #############################################################
      ###############################################################################
      
      #SAVE AND KEEP THE VARIOUS PRODUCED DATAFRAMES (outside of the REPLICATE and PERIOD)
      # the following two objects are used in the CSI score determination
      IF_frames <- list(IF_frames); names(IF_frames) <- paste(c("IF_frames",REPLICATE,PERIOD),collapse="_")
      list_IF_Frames <- c(list_IF_Frames,IF_frames)
      rm(list=ls()[which(ls()%in%c("extracted_variables_AUTO","extracted_variables_MANUAL","IF_frames"))])
    }##PERIOD
    replicate_full_ant_list <- list(replicate_full_ant_list); names(replicate_full_ant_list) <- paste(c("replicate_list",REPLICATE),collapse="_")
    list_replicate_full_ant_list <- c(list_replicate_full_ant_list,replicate_full_ant_list)
    
    replicate_focal_list <- list(replicate_focal_list); names(replicate_focal_list) <- paste(c("replicate_list",REPLICATE),collapse="_")
    list_replicate_focal_list <- c(list_replicate_focal_list,replicate_focal_list)
    
  }##REPLICATE
  
  return(list(  
    # interaction_MANUAL=interaction_MANUAL,
    summary_MANUAL=summary_MANUAL,
    # interaction_AUTO=interaction_AUTO,
    summary_AUTO=summary_AUTO,
    list_IF_Frames=list_IF_Frames,
    list_replicate_full_ant_list=list_replicate_full_ant_list,
    list_replicate_focal_list=list_replicate_focal_list
    
    )
  )
}


extract_from_object <- function (  object, IF_frames, positions, BEH , REPLICATE, PERIOD,selected_variables=NULL  ){

  
  if (!"unique_interaction_id"%in%names(object)){ object$unique_interaction_id <- 1:nrow(object)}
  if ("REPLICATE" %in% names(object)){ object <- object[which(object$REPLICATE == REPLICATE),]}
  if ("PERIOD"    %in% names(object)){ object <- object[which(object$PERIOD    == PERIOD),]}
  if ("Behaviour" %in% names(object)){ object <- object[which(object$Behaviour == BEH),]}
  
  
  name_list_summary <- NULL
  name_list_traj    <- NULL
  
  summary_variables  <- NULL
  stored_trajectories <- NULL
  
  if (is.null(selected_variables)){
    selected_variables_pair       <- NULL
    selected_variables_individual <- NULL 
  }else{
    selected_variables_pair       <- selected_variables[which(!grepl("ACT|REC|duration",selected_variables))]
    selected_variables_individual <- unique(c("transfmean_speed_BLpersec",gsub("_ACT","",gsub("_REC","",selected_variables[which(grepl("ACT|REC",selected_variables))]))   ))
  }
  
  
  ## loop through each event in object
  to_keep_annotations <- c("to_keep_annotations",ls(),"ROW","object")
  
  ###sort object by decreasing order of event duration to ensure that we don't start with a trajectory that has too few items to be processed
  object <- object[order(-object$duration),]
  
  for (ROW in 1:nrow(object))
  { 
    rm(list=ls()[which(!ls()%in%to_keep_annotations)])
    gc()
    
    if ((ROW/1000)  == round(ROW/1000)){print(paste("Extracting movement variables for interaction",ROW,"out of",nrow(object)))}
    #######get time information about event
    duration_sec   <-  object[ROW,"T_stop_sec"] - object[ROW,"T_start_sec"] + 1/FRAME_RATE 
    frame_start   <-  object[ROW,"frame_start"]
    frame_stop     <-  object[ROW,"frame_stop"]
    duration_frames <-  frame_stop - frame_start + 1
    
    ####################################################################
    ###for manual annotations: create ant1 and ant2 columns
    if (!"ant1"%in%names(object)){
      object$ant1 <- object$Actor
      object$ant2 <- object$Receiver
    }
    
    
    ## extract actor, receiver IDs & start & end times from the hand-annotated data
    ANT1 <- object$ant1[ROW]    ;     ANT1_Name <- paste("ant",ANT1,sep="_")
    ANT2 <- object$ant2[ROW] ;     ANT2_Name <- paste("ant",ANT2,sep="_")
    # if (is.na(ANT2)){
    #   print(paste("Behaviour:",BEH,"number",ROW,"Ant1:",ANT1_Name))
    # }else{
    #   print(paste("Behaviour:",BEH,"number",ROW,"Ant1:",ANT1_Name,"Ant2:",ANT2_Name))
    # }
    
    ###now for each ant extract individual movement characteristics
    for (AntType in c("ANT1","ANT2")){
      summary_individual_traj <- NULL; traj <- NULL
      ###extract trajectory of desired ant
      traj <- positions$trajectories[[get(paste(AntType,"_Name",sep=""))]]
      
      ###crop traj to desired times, add frame and time information, calculate time interval and distance moved whichg will be used for filtering
      traj <- process_traj(traj,IF_frames, frame_start, frame_stop)
      
      ###then run script to extract individual variables
      summary_individual_traj <- suppressWarnings(extract_individual_movement_variables(traj,frame_start, frame_stop,selected_variables_individual))
      traj <- summary_individual_traj[["traj"]];summary_individual_traj <- summary_individual_traj[["summary_individual_traj"]]
      if (is.null(name_list_summary)&!is.null(summary_individual_traj)){
        name_list_summary <- names(summary_individual_traj)
      }
      
      if (is.null(summary_individual_traj)){
        summary_individual_traj        <-as.data.frame(matrix(nrow=1,ncol=length(name_list_summary)))
        names(summary_individual_traj) <- name_list_summary
      }
      
      ###prepare to store:
      names(traj)[which(names(traj)!="frame")] <- paste(AntType,names(traj)[which(names(traj)!="frame")],sep=".")
      if (ncol(summary_individual_traj)>0){
        names(summary_individual_traj)           <- paste(names(summary_individual_traj),AntType,sep="_")
      }else{
        summary_individual_traj <- NULL
      }
      assign(paste("traj",AntType,sep="_"),traj); rm(list=c("traj"))
      assign(paste("summary_individual_traj",AntType,sep="_"),summary_individual_traj);rm(list=c("summary_individual_traj"))
    }

    ####Reassign ANT1 and ANT2 as actors and Receivers
    if ("Actor" %in% names(object)){
      ACT <- c("ANT1","ANT2")[which(c(ANT1,ANT2)==object$Actor[ROW])]
      REC <- c("ANT1","ANT2")[which(c(ANT1,ANT2)==object$Receiver[ROW])]
    }else{
      if (is.na(ANT1)&!is.na(ANT2)){
        ACT <- "ANT2"
        REC <- "ANT1"
      }else if (!is.na(ANT1)&is.na(ANT2)){
        ACT <- "ANT1"
        REC <- "ANT2"
      }else{
        ACT <- c("ANT1","ANT2") [which.max(c(summary_individual_traj_ANT1$transfmean_speed_BLpersec_ANT1,summary_individual_traj_ANT2$transfmean_speed_BLpersec_ANT2))]
        REC <- c("ANT1","ANT2") [which( c("ANT1","ANT2")!=ACT)]
      }
    }
    
    ###do the rest only if ACT and REC have been reassigned (i.e., exclude rare interactions where there is not information to even get the speed of each)
    if (length(ACT)>0){
      for (AntType in c("ACT","REC")){
        traj                    <- get(paste("traj",get(AntType),sep="_"));rm(list=paste("traj",get(AntType),sep="_"))
        summary_individual_traj <- get(paste("summary_individual_traj",get(AntType),sep="_"));rm(list=paste("summary_individual_traj",get(AntType),sep="_"))
        
        names(traj)                    <- gsub(get(AntType),AntType,names(traj))
        names(summary_individual_traj) <- gsub(get(AntType),AntType,names(summary_individual_traj))
        
        assign(paste("traj",AntType,sep="_"),traj); rm(list=c("traj"))
        assign(paste("summary_individual_traj",AntType,sep="_"),summary_individual_traj); rm(list=c("summary_individual_traj"))
        
        assign(paste(AntType,"_Name",sep=""),get(paste(get(AntType),"_Name",sep=""))); rm(list=paste(get(AntType),"_Name",sep=""))
      }
      
      ###prepare summary output line including both ants - containing individual variables
      summary_ROW <- data.frame(REPLICATE, 
                                PERIOD, 
                                BEH=BEH, 
                                ROW,
                                unique_interaction_id=object[ROW,"unique_interaction_id"],
                                Act_Name=ACT_Name, 
                                Rec_Name=REC_Name,
                                ant1 = as.numeric(gsub("ant_","", ACT_Name)), 
                                ant2=  as.numeric(gsub("ant_","", REC_Name)), 
                                frame_start               , 
                                frame_stop, 
                                duration_sec,
                                summary_individual_traj_ACT,
                                summary_individual_traj_REC)
      summary_ROW$ROW <- as.factor(summary_ROW$ROW )
      
      ###then pair variables
      if (is.null(selected_variables_pair)|length(selected_variables_pair)>0){
        summary_BOTH <- extract_pair_movement_variables(traj_ACT,traj_REC)
        # traj_BOTH <- summary_BOTH[["traj_BOTH"]]; 
        summary_BOTH <- summary_BOTH[["summary_both"]]
        
        ##add pair information to summary_ROW
        summary_ROW <- cbind(summary_ROW,summary_BOTH)
      }

      
      ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
      summary_ROW$pair <- apply(summary_ROW[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
      
      #create object to hold trajectories
      # interacts_ROW <- data.frame(REPLICATE=REPLICATE,
      #                             PERIOD=PERIOD,
      #                             ROW=ROW,
      #                             BEH=BEH,
      #                             Act_Name=ACT_Name,Rec_Name=REC_Name,
      #                             ant1 = as.numeric(gsub("ant_","", ACT_Name)), 
      #                             ant2=  as.numeric(gsub("ant_","", REC_Name)), 
      #                             traj_BOTH,
      #                             stringsAsFactors = F)
      # interacts_ROW$ROW <- as.factor(interacts_ROW$ROW )
      # interacts_ROW$pair <- apply(interacts_ROW[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
      
      ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
      # stored_trajectories <- rbind(stored_trajectories, interacts_ROW)
      summary_variables   <- rbind(summary_variables,     summary_ROW)
    }
    
  }##ROW
  # return(list(summary_variables=summary_variables,stored_trajectories=stored_trajectories))
  return(list(summary_variables=summary_variables))
  
}

process_traj <- function(traj,IF_frames, frame_start, frame_stop){
  ###crop traj to the desired time
  traj <- traj [ which(traj$frame >= frame_start & traj$frame <= frame_stop),]
  
  if (is.null(traj) | (length(traj$frame)==0)){ ###if traj is NULL or has no lines
    traj <- data.frame(x=numeric(),y=numeric(),angle=numeric(),zone=numeric(),frame=numeric(),time_sec=numeric(),time_interval=numeric(),dt_FRAME=numeric(),distance_moved_px=numeric())
  }else{
    ## remove 'time' column as it is confusing - it's not a common time
    traj$time <- NULL
    
    ###add a time_sec column containing common tme from IF_frames
    traj$time_sec      <- IF_frames[match (traj$frame,IF_frames$frame_num),"time_sec"]
    
    ###calculate successive time intervals and frame intervals
    traj               <- traj[order(traj$time_sec),]
    traj$time_interval <- c(diff(traj$time_sec),NA) ###time difference between t and t+1
    traj$dt_FRAME      <- c(diff(traj$frame)   ,NA) ## time difference in frames between time t and time t+1
    
    ###calculate distance moved between successive fixes
    traj$distance_moved_px <-  c(with(traj, (sqrt(diff(x)^2 + diff(y)^2))),NA) # euclidean distance
    
  }
  return(traj)
  
}

extract_individual_movement_variables <- function(traj,frame_start, frame_stop,selected_variables_individual=NULL){
  if (is.null(selected_variables_individual)){
    selected_variables_individual <- 
      c(
      "sum_moved_distance_BL"                  
      ,"transfmean_speed_BLpersec"             
      ,"transfSD_speed_BLpersec"               
      ,"transfmean_abs_accel_BLpersec2"       
      ,"transfSD_abs_accel_BLpersec2"         
      ,"transfmean_abs_jerk_BLPerSec3"        
      ,"transfSD_abs_jerk_BLPerSec3"          
      ,"transfmean_abs_TurnAngle"             
      ,"transfSD_abs_TurnAngle"               
      ,"stDev_turnAngle"                      
      ,"transfmean_abs_ang_Velocity_Movement" 
      ,"transfSD_abs_ang_Velocity_Movement"   
      ,"stDev_ang_Velocity_Movement"          
      ,"transfmean_abs_Body_Rotation"         
      ,"transfSD_abs_Body_Rotation"           
      ,"stDev_Body_Rotation"                  
      ,"transfmean_abs_ang_Velocity_Body"     
      ,"transfSD_abs_ang_Velocity_Body"       
      ,"stDev_ang_Velocity_Body"              
      ,"StDev_Body_angle"                     
      ,"StDev_Movement_angle"                 
      ,"transfmean_abs_Movement_Body_angle_diff"
      ,"transfSD_abs_Movement_Body_angle_diff"  
      ,"stDev_Movement_Body_angle_diff"         
      ,"transfmean_abs_Movement_Body_Inclination_angle"
      ,"transfSD_abs_Movement_Body_Inclination_angle"  
      ,"stDev_Movement_Body_Inclination_angle"         
      ,"root_mean_square_deviation_BL2"         
      ,"chull_area_BL2"                            
      ,"prop_time_undetected"   
    )
    
  }
  ##########################################################
  ## INDIVIDUAL TRAJECTORY frame-by-frame MEASURES #########
  ##########################################################
  ### Run cpp function to add salient angles:
  ### TurnAngle: change in movement direction between (t-1 to t) step and (t to t+1) step; from minus pi to pi
  ### Movement_direction: direction of movement from t to t+1, from minus pi to pi
  ### Movement_Body_angle_diff: difference between movement direction from t to t+1 and body angle at time t, from minus pi to pi
  ### Movement_Body_Inclination_angle: (acute) angle between (undirected) inclination of movement and (undirected) inclination of body; from 0 to pi/2
  ### Body_Rotation: change of body orientation between time t and t+1 
  traj <- cbind(traj, add_angles (traj))  
  
  #convert distance moved from px to BL
  traj$distance_moved_BL <- traj$distance_moved_px/traj$body_length_pixels
  
  #instantaneous speed, acceleration, and jerk(diff in accelerations) 
  traj$speed_BLPerSec  <- traj$distance_moved_BL  / traj$time_interval
  
  traj$accel_BLPerSec2 <- c(diff(traj$speed_BLPerSec) ,NA) / traj$time_interval
  traj$jerk_BLPerSec3  <- c(diff(traj$accel_BLPerSec2),NA) / traj$time_interval
  
  # angular velocity: ω = (α₂ - α₁) / t = Δα / t
  traj$ang_Velocity_Body     <- traj$Body_Rotation        / traj$time_interval
  traj$ang_Velocity_Movement <- traj$TurnAngle            / traj$time_interval
  
  ##########################################
  ## SUMMARY MEASURES ######################
  ##########################################
  ## Apply THRESHOLDs to exclude gaps in which the individuals were not detected for very long (FRAME)
  ## remove movement variables when traj$dt_FRAME > DT_frame_THRESHOLD
  # dt_FRAME corresponds to frame between T and T+1. This should be applied to all movement characteristics
  traj[which(traj$dt_FRAME>DT_frame_THRESHOLD), which(!names(traj)%in%c("x","y","angle","zone","body_length_pixels","frame","time_sec","time_interval","dt_FRAME"))] <- NA
  
  ## Apply THRESHOLDs to exclude jitter in the individuals' movement (DISTANCE)
  # remove movement variables when traj$distance_moved < DT_dist_THRESHOLD_BL
  
  traj[which(traj$distance_moved_BL<DT_dist_THRESHOLD_BL), which(!names(traj)%in%c("x","y","angle","zone","body_length_pixels","frame","time_sec","time_interval","dt_FRAME"))] <- NA
  
  ###now if traj has enough rows, compute summary variables
  if (nrow(traj)>1){
    summary_individual_traj        <- data.frame (matrix(ncol = length(selected_variables_individual), nrow = 1))
    names(summary_individual_traj) <-selected_variables_individual
    for ( func in selected_variables_individual){
      if (func !="prop_time_undetected"){
        summary_individual_traj[1,func] <-get(func)(traj)
      }else{
        summary_individual_traj[1,func] <-get(func)(traj,frame_start,frame_stop)
      }
    }
   }else{summary_individual_traj <- NULL}
  
  return(list(traj=traj,summary_individual_traj=summary_individual_traj))
}

extract_pair_movement_variables       <- function(traj_ACT,traj_REC){
  traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("frame"))  ## 21 Jan 2022: Changed to sue raw unix seconds to allow simpler matching below (posix formats cause problems with the matching...)
  #straight line - euclidean distance
  traj_BOTH$straightline_pair_dist_BL    <-  sqrt(((traj_BOTH$ACT.x-traj_BOTH$REC.x)/traj_BOTH$ACT.body_length_pixels)^2+((traj_BOTH$ACT.y-traj_BOTH$REC.y)/traj_BOTH$ACT.body_length_pixels)^2)
  traj_BOTH$pair_Body_Orient_diff_abs    <-  abs((traj_BOTH$ACT.angle  - traj_BOTH$REC.angle) %% pi)
  traj_BOTH$pair_Body_Inclination_diff   <-  unlist(lapply(traj_BOTH$pair_Body_Orient_diff , FUN=inclination_angle))###NATH_FLAG:do we want inclination here?
  summary_both <- data.frame(mean_strghtline_pair_dist_BL    = mean(traj_BOTH$straightline_pair_dist_BL , na.rm=T),
                             mean_pair_body_orient_diff_abs  = mean(traj_BOTH$pair_Body_Orient_diff_abs , na.rm=T),
                             mean_pair_body_inclination_diff = mean(traj_BOTH$pair_Body_Inclination_diff, na.rm=T) 
  )
  # return(list(traj_BOTH=traj_BOTH,summary_both=summary_both))
  return(list(summary_both=summary_both))
}

add_body_length_to_traj <- function(positions,body_lengths){
  for (traj_segment in 1:nrow(positions$trajectories_summary)){
    space_body_length <- body_lengths[which(body_lengths$space_ID==positions$trajectories_summary[traj_segment,"space"]),"body_length"]
    if (length(space_body_length)==1){
      positions$trajectories[[traj_segment]]$body_length_pixels <- space_body_length
    }
  }
  return(positions)
}

chull_area_BL2 <- function(traj){
  # Convex Hull
  box.coords <- traj[, c("x", "y","body_length_pixels")]
  box.coords$x <- box.coords$x/box.coords$body_length_pixels
  box.coords$y <- box.coords$y/box.coords$body_length_pixels
  box.coords <- as.matrix(box.coords[c("x","y")])
  
  box.hpts <- chull(x = traj$x/traj$body_length_pixels, y = traj$y/traj$body_length_pixels) # calculate convex hull for x and y columns
  box.hpts <- c(box.hpts, box.hpts[1]) #add first coord at the bottom to close area
  box.chull.coords <- box.coords[box.hpts,]
  chull.poly <- Polygon(box.chull.coords, hole=F); 
  return(chull.poly@area)
}

sum_moved_distance_BL <- function(traj){return(sum(traj$distance_moved_BL   ,  na.rm = T))}

get_body_lengths <- function(e,all_body_lengths) {
  ###check how many spaces there are
  space_list <- e$spaces
  body_lengths <- NULL
  
  for (space_ID in 1:length(space_list)){
    body_lengths <- rbind(body_lengths,data.frame(space_ID        = space_list[[space_ID]]$ID,
                                                  tracking_system = space_list[[space_ID]]$name,
                                                  body_length     = all_body_lengths[which(all_body_lengths$TS==space_list[[space_ID]]$name),"mean_worker_length_px"]
    ))
  }
  return(body_lengths)
}

transfmean_speed_BLpersec             <- function(traj){return( exp(mean  (log(na.omit(traj$speed_BLPerSec)                           ))))}
transfSD_speed_BLpersec               <- function(traj){return( exp(sd  (log(na.omit(traj$speed_BLPerSec)                           ))))}


transfmean_abs_accel_BLpersec2        <- function(traj){return( exp(mean  (log(na.omit(abs(traj$accel_BLPerSec2))                           ))))}
transfSD_abs_accel_BLpersec2          <- function(traj){return( exp(sd  (log(na.omit(abs(traj$accel_BLPerSec2))                           ))))}


transfmean_abs_jerk_BLPerSec3         <- function(traj){return( exp(mean  (log(na.omit(abs(traj$jerk_BLPerSec3))                           ))))}
transfSD_abs_jerk_BLPerSec3           <- function(traj){return( exp(sd  (log(na.omit(abs(traj$jerk_BLPerSec3))                           ))))}

transfmean_abs_TurnAngle              <- function(traj){return( (mean  (sqrt(na.omit(abs(traj$TurnAngle))                           )))^2)}
transfSD_abs_TurnAngle                <- function(traj){return(  (sd  (sqrt(na.omit(abs(traj$TurnAngle))                           )))^2)}
stDev_turnAngle                       <- function(traj){return( angular.deviation(traj$TurnAngle              ,  na.rm = T))}

transfmean_abs_ang_Velocity_Movement  <- function(traj){return( (mean  (sqrt(na.omit(abs(traj$ang_Velocity_Movement))                           )))^2)}
transfSD_abs_ang_Velocity_Movement    <- function(traj){return(  (sd  (sqrt(na.omit(abs(traj$ang_Velocity_Movement))                           )))^2)}
stDev_ang_Velocity_Movement           <- function(traj){return( angular.deviation(traj$ang_Velocity_Movement              ,  na.rm = T))}

transfmean_abs_Body_Rotation          <- function(traj){return( (mean  (sqrt(na.omit(abs(traj$Body_Rotation))                           )))^2)}
transfSD_abs_Body_Rotation            <- function(traj){return(  (sd  (sqrt(na.omit(abs(traj$Body_Rotation))                           )))^2)}
stDev_Body_Rotation                   <- function(traj){return( angular.deviation(traj$Body_Rotation              ,  na.rm = T))}


transfmean_abs_ang_Velocity_Body      <- function(traj){return( (mean  (sqrt(na.omit(abs(traj$ang_Velocity_Body))                           )))^2)}
transfSD_abs_ang_Velocity_Body        <- function(traj){return(  (sd  (sqrt(na.omit(abs(traj$ang_Velocity_Body))                           )))^2)}
stDev_ang_Velocity_Body               <- function(traj){return( angular.deviation(traj$ang_Velocity_Body              ,  na.rm = T))}

StDev_Body_angle                      <- function(traj){return( angular.deviation(traj$angle                  ,  na.rm = T))}
StDev_Movement_angle                  <- function(traj){return( angular.deviation(traj$Movement_direction     ,  na.rm = T))}

transfmean_abs_Movement_Body_angle_diff    <- function(traj){return( (mean  (sqrt(na.omit(abs(traj$Movement_Body_angle_diff))                           )))^2)}
transfSD_abs_Movement_Body_angle_diff      <- function(traj){return(  (sd  (sqrt(na.omit(abs(traj$Movement_Body_angle_diff))                           )))^2)}
stDev_Movement_Body_angle_diff           <- function(traj){return( angular.deviation(traj$Movement_Body_angle_diff              ,  na.rm = T))}

transfmean_abs_Movement_Body_Inclination_angle    <- function(traj){return( (mean  (sqrt(na.omit(traj$Movement_Body_Inclination_angle)                           )))^2)}
transfSD_abs_Movement_Body_Inclination_angle      <- function(traj){return(  (sd  (sqrt(na.omit(traj$Movement_Body_Inclination_angle)                           )))^2)}
stDev_Movement_Body_Inclination_angle           <- function(traj){return( angular.deviation(traj$Movement_Body_Inclination_angle              ,  na.rm = T))}

root_mean_square_deviation_BL2         <- function(traj){return( sqrt((sum(((traj$x/traj$body_length_pixels)-mean((traj$x/traj$body_length_pixels),  na.rm = T))^2+((traj$y/traj$body_length_pixels)-mean((traj$y/traj$body_length_pixels),  na.rm = T))^2))/(length(na.omit((traj$x/traj$body_length_pixels))))))} ### Root Mean Square Deviation,  as defined in Ulrich et al 2018)}
prop_time_undetected                  <- function(traj,frame_start,frame_stop){return( 1-(sum(!is.na(traj$x))/(frame_stop-frame_start+1)))}
