#########################################################
################## GET ANT TASK #########################
#########################################################
 
#### THIS VERSION IS FORT 0.8.1 COMPATIBLE ####
 
# Script created by Adriano Wanderlingh, with contributions by  Nathalie Stroeymeyt
# AntTasks_Linda was created by NS for Linda Sartoris


#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore:
# https://github.com/AdrianoWanderlingh/Ant_Tracking


#### THIS SCRIPT IS USED INSIDE:
# #### - Extract_Metadata_v082.R
# 
# AntTasks <- function(e){
#   print("Computing AntTasks based on 48h time-window before exposure")
#   
#   #required packages
#   require(FortMyrmidon)
#   require(mallinfo)
#   
#   #PARAMS for function
#   #hour_chunk_start represents the start time (end time defined by TimeWindow)
#   hour_chunk_start <- sort(seq(39, 75, by = 12),decreasing = TRUE ) # c(75,63,51,39)
#   TimeWindow <- 12
#   
# 
#   #define zones
#   #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
#   zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
#   zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
#   foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
#   nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
#   print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
#   
#   #Get complete list of Ants
#   AntID_list <- NULL
#   for (ant in   1: length(e$ants)) {
#     AntID_list <- c(AntID_list,e$ants[[ant]]$ID)
#     }
# 
#   
#   ######L: use Nathalie's function to extract trajectories from start to end #####
#   # source trajectory_extraction.R
#   source(paste(BEH_FUNCTIONS,"trajectory_extraction.R",sep="/"))
#   #postions <- extract_trajectories ()
#   
#   
#   
#   
#   
#   # positions_summaries_list <- list()
#   # loop_N <- 0
#   # 
#   # for (HOUR_start in hour_chunk_start) {
#   #   loop_N <- loop_N + 1  
#   #   print(paste0("Computing chunk ",loop_N," of 4"))
#   #   
#   #   ## get 2 12Hours window for the Task calculation
#   #   ## calcualte the task BEFORE the EXPOSURE
#   #   time_start <- fmQueryGetDataInformations(e)$end - HOUR_start*3600
#   #   #time_start <- fmQueryGetDataInformations(e)$start + 33*3600 ####first time in tracking plus 21 hours, to skip acclimation time + 12 HOURS
#   #   time_start <- fmTimeCreate(offset=time_start)
#   #   #time_start <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$time_stop - 24*3600)####last time in tracking minus 24 hours
#   #   stop <- fmQueryGetDataInformations(e)$end - (HOUR_start-TimeWindow)*3600
#   #   #stop  <- fmQueryGetDataInformations(e)$start + 45*3600 ####pre-tracking period 
#   #   time_stop   <- fmTimeCreate(offset=stop)
#   #   #time_stop  <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$end) ####last time in tracking
#   #   ###QUERY 3: fmQueryComputeAntTrajectories()
#   #   positions                 <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
#     positions_summaries       <- positions$trajectories_summary
#     positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
#     positions_list            <- positions$trajectories
#     
#     #for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
#     # if so the ant is a forager, if not the ant is a nurse
#     positions_summaries$AntTask <- NA
#     
#     ##before going back and forth between positions_summaries and positions_list:
#     ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
#     nrow(positions_summaries)==length(positions_list)
#     ####2/ always make sure that positions_summaries is ordered correctly, using the index column
#     positions_summaries <- positions_summaries[order(positions_summaries$index),]
#     ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
#     for ( ant_index in 1:length(positions_list)) {
#       positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
#       positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
#       
#       if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
#         positions_summaries[ant_index,"AntTask"] <- "forager"
#       }else{
#         positions_summaries[ant_index,"AntTask"] <- "nurse"
#       }
#     }
#     #match antID and tagID (live tracking gives tagID). 
#     IDs <- e$identificationsAt(fmTimeNow()) #this skips dead ants
#     IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
#     IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
#     positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
#     #positions_summaries1 <- positions_summaries
#     #positions_summaries_LOOP <- aggregate(cbind(positions_summaries$nb_frames_outside,positions_summaries$nb_frames_inside), by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries_LOOP) [match("Group.1",colnames(positions_summaries_LOOP))] <- "antID"
#     positions_summaries_LOOP <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
#     #add names that help in merging later on
#     colnames(positions_summaries_LOOP) [match("nb_frames_outside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_outside",loop_N)
#     colnames(positions_summaries_LOOP) [match("nb_frames_inside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_inside",loop_N)
#   #   
#   #   
#   #   #positions_summaries_list <- append(positions_summaries_list, positions_summaries_LOOP)
#   #   positions_summaries_list[[loop_N]] <-  positions_summaries_LOOP
#   #   
#   #   
#   #   rm(list=ls()[which(!ls()%in%c("positions_summaries_list","e","foraging_zone","nest_zone","AntID_list","hour_chunk_start","TimeWindow","loop_N"))]) #close experiment
#   #   gc()
#   #   mallinfo::malloc.trim(0L)
#   #   
#   # }
#   
#   #merge all data frames together
#   
#   # # non-good looking recursive merging
#   # positions_summaries_mergA <- merge(positions_summaries_list[[1]][c("antID","nb_frames_outside1","nb_frames_inside1")], # , "tag_hex_ID"
#   #                                    positions_summaries_list[[2]][c("antID","nb_frames_outside2","nb_frames_inside2")], # , "tag_hex_ID"
#   #                                    all.x=T,all.y=T)
#   # 
#   # positions_summaries_mergB <-  merge(positions_summaries_list[[3]][c("antID","nb_frames_outside3","nb_frames_inside3")], # , "tag_hex_ID"
#   #                                     positions_summaries_list[[4]][c("antID","nb_frames_outside4","nb_frames_inside4")], # , "tag_hex_ID"
#   #                                     all.x=T,all.y=T)
#   
#   # positions_summaries <- Reduce(function(x, y) merge(x, y, all=TRUE), positions_summaries_list)
#   
#   positions_summaries <- as.data.frame(sapply(positions_summaries,as.numeric))
#   
#   positions_summaries[is.na(positions_summaries)] <- 0
#   
#   #positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2)/(positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2+positions_summaries$nb_frames_inside1+positions_summaries$nb_frames_inside2)
#   
#   #sum inside & outside
#   positions_SUMS<- data.frame(antID=positions_summaries$antID,tag_hex_ID=positions_summaries$tag_hex_ID, outside=rowSums(positions_summaries[, grep("outside", colnames(positions_summaries))]),
#                               inside=rowSums(positions_summaries[, grep("inside", colnames(positions_summaries))]))
#   
#   positions_SUMS$prop_time_outside <- positions_SUMS$outside/(positions_SUMS$outside+positions_SUMS$inside)
#   
#   positions_SUMS[which(positions_SUMS$prop_time_outside<=0.01),"AntTask"] <- "nurse"
#   positions_SUMS[which(positions_SUMS$prop_time_outside>0.01),"AntTask"] <- "forager"
#   
#   AntTasks <- data.frame(antID=positions_SUMS[,"antID"],
#                          tag_hex_ID=positions_SUMS[,"tag_hex_ID"],
#                          AntTask= positions_SUMS[,"AntTask"],
#                          prop_time_outside= positions_SUMS[,"prop_time_outside"])
#   
#   print("AntTasks computed")
#   
#   # #add missing ants
#   missing_ants <- subset(AntID_list, !(AntID_list %in% AntTasks$antID))
#   missing_ants_table <- data.frame()
#   
#   for (MISSING in missing_ants) {
#   for (id in length(e$ants[[MISSING]]$identifications)) {
#    #print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
#     missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING, tag_hex_ID= e$ants[[MISSING]]$identifications[[id]]$tagValue ,AntTask=NA, prop_time_outside= NA))
#   }}
#   
#   AntTasks <- rbind(AntTasks, missing_ants_table)
#   # add misssing ants  as NURSE by DEFAULT
#   # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"
#   
#   AntTasks$AntTask_num <- NA
#   AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
#   AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
#   AntTasks <- AntTasks[order(AntTasks$antID),]
#   
#   rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","e","foraging_zone","AntID_list","AntTasks"))]) #close experiment
#   gc()
#   mallinfo::malloc.trim(0L)
#   
#   ##RETURN OUTPUT
#   # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
#   return(AntTasks)
#   
# } #, ZoneUsage= TRUE


AntTasks_Linda <- function(e){
  print("Computing AntTasks based on 24h time-window before exposure")
  
  #required packages

  
  #PARAMS for function
  #define zones
  #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
  zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
  zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
  foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
  
  #Get complete list of Ants
  ####L: in AntID_list add a line for tag hex ID ####
  AntID_list <- NULL
  for (ant in   1: length(e$ants)) {
    AntID_list <- rbind(AntID_list,
                        data.frame (
                          antID = e$ants[[ant]]$ID
                          ,
                          tag_hex_ID = e$ants[[ant]]$identifications[[1]]$tagValue
                          ,
                          alive = e$ants[[ant]]$getValue("IsAlive",fmTimeNow())
                          ,
                          treated = e$ants[[ant]]$getValue("Exposed",fmTimeNow())
                        ))
  }
  
  #### L: in metadata_colonies time is no longer in Zulu format but already GMT, so change the following lines: ####
  ###extract start and end of pre-treatment period from metadata to define the period on which ant tasks will be determined
  #L: this would be the time before I started with the forager collection
  time_start <- metadata_colonies[    which   (     grepl(   unlist(  strsplit(REP.FILES,split="/"))[1] ,  metadata_colonies$treat_TS_colony   )      )     ,   "exp_start"     ]
  time_end <- metadata_colonies[    which   (     grepl(   unlist(  strsplit(REP.FILES,split="/"))[1] ,  metadata_colonies$treat_TS_colony   )      )     ,   "middle_end"     ]
  ###turn time_start and time_end from character strings copied from fort-studio to time objects that R understands
  #time_start <- as.POSIXct(time_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  #time_end   <- as.POSIXct(time_end, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  #L: when time stamps were calculated in R use this:
  time_start <- as.POSIXct(time_start, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
  time_end   <- as.POSIXct(time_end, format = "%Y-%m-%d %H:%M:%S",  origin="1970-01-01", tz="GMT" )
  ###turn time_start and time_end from R time objects to FortMyrmidon time objects 
  time_start <- fmTimeCreate(time_start)
  time_end  <- fmTimeCreate(time_end)
  
  IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_end,showProgress = FALSE)
  IF_frames           <- IdentifyFrames$frames
  IF_frames           <- IF_frames[order(IF_frames$time), ] ## this was added by Daniel because of error in trajectory_extraction.R :'vec' must be sorted non-decreasingly and not contain NAs #### i.e. for NS_guillam_c08
  # Assign a frame to each time since start and use it as baseline for all matching and computation
  IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))
  
  # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
  IF_frames$time_sec <- round(as.numeric(IF_frames$time),3)
  
  
  
  ######L: use Nathalie's function to extract trajectories from start to end #####
positions <- extract_trajectories(e = e,
                                  start = time_start,
                                  end = time_end,
                                  maximumGap = fmHour(24*365) ,
                                  IF_frames = IF_frames) 
  positions_summaries       <- positions$trajectories_summary
  positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
  positions_list            <- positions$trajectories
  
  #for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
  # if so the ant is a forager, if not the ant is a nurse
  
  ##before going back and forth between positions_summaries and positions_list:
  ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
  nrow(positions_summaries)==length(positions_list)
  ####2/ always make sure that positions_summaries is ordered correctly, using the index column
  positions_summaries <- positions_summaries[order(positions_summaries$index),]
  ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
  for ( ant_index in 1:length(positions_list)) {
    positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
    positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
  }
  positions_summaries$nb_frames_inside[which(is.na(positions_summaries$nb_frames_inside))] <-0
  positions_summaries$nb_frames_outside[which(is.na(positions_summaries$nb_frames_outside))] <-0
  positions_summaries$prop_time_outside <- positions_summaries$nb_frames_outside/(positions_summaries$nb_frames_outside+positions_summaries$nb_frames_inside)
  #L: add a column of ants that moved outside (to check how many ants were outside ~constantly) ? ####
  # positions_summaries$moved_outside <- ifelse(positions_summaries$prop_time_outside >= 0.9, "yes", "no")
  
  positions_summaries[which(positions_summaries$prop_time_outside<=0.01),"AntTask"] <- "nurse"
  positions_summaries[which(positions_summaries$prop_time_outside>0.01),"AntTask"]  <- "forager"
  positions_summaries[which(positions_summaries$antID==1),"AntTask"]                <- "queen"
  
  merged_data <- merge(positions_summaries, AntID_list, by = "antID", all.x = TRUE)
  positions_summaries$tag_hex_ID <- merged_data$tag_hex_ID
  
  AntTasks <- data.frame(antID=positions_summaries[,"antID"],
                         tag_hex_ID=positions_summaries[,"tag_hex_ID"],
                         AntTask= positions_summaries[,"AntTask"],
                         prop_time_outside= positions_summaries[,"prop_time_outside"])
  
  print("AntTasks computed")
  
  # #add missing ants
  missing_ants <- subset(AntID_list$antID, !(AntID_list$antID %in% AntTasks$antID))
  missing_ants_table <- data.frame()
  
  for (MISSING in missing_ants) {
    for (id in length(e$ants[[MISSING]]$identifications)) {
      #print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
      missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING, tag_hex_ID= e$ants[[MISSING]]$identifications[[id]]$tagValue ,AntTask=NA, prop_time_outside= NA))
    }}
  
  AntTasks <- rbind(AntTasks, missing_ants_table)
  # add misssing ants  as NURSE by DEFAULT
  # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"
  
  AntTasks$AntTask_num <- NA
  AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
  AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
  AntTasks <- AntTasks[order(AntTasks$antID),]
  
  AntTasks <- merge(AntTasks,AntID_list,all.x=T)
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","e","foraging_zone","AntID_list","AntTasks"))]) #close experiment
  gc()
  mallinfo::malloc.trim(0L)
  
  ##RETURN OUTPUT
  # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
  return(AntTasks)
  
} #, ZoneUsage= TRUE
