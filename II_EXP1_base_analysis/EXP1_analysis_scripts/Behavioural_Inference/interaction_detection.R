chunk_size_for_merge <- 1000

collapse_collisions <- function(collisions){
  ###bring collisions metadata and zone/x-y/time info to same dataframe
  ###first, merge collisions$collisions with collisions$frames
  collisions$collisions <- cbind(collisions$collisions,collisions$frames[collisions$collisions$frames_row_index,])
  
  ###second, add information from positions to collisions$collisions
  collisions$collisions <- cbind(collisions$collisions,t(apply (collisions$collisions,1,get_positional_data,collisions=collisions)))
  
  return(collisions$collisions)
}

get_positional_data <- function (x,frames_row_index="frames_row_index",ant_pair = c("ant1","ant2"),collisions){
  pos_data <- unlist(collisions$positions[[as.numeric(x[frames_row_index])]][match(as.numeric(x[ant_pair]),as.numeric(collisions$positions[[as.numeric(x[frames_row_index])]]$antID)),])
  return(pos_data[which(!names(pos_data)%in%c("antID1","antID2"))])
}

filter_capsules     <- function(collisions,capsule_matcher){
  ###first, list types of interactions which we want to keep
  to_include <- c()
  for (caps in 1:length(capsule_matcher)){
    to_include <- c(to_include,c(paste(capsule_matcher[[caps]],collapse="-"),paste(rev(capsule_matcher[[caps]]),collapse="-")))
  }
  to_include <- sort(unique(to_include))
  
  ###use to_include and the content of collisions$types to list the interactions we ARE NOT INTERESTED IN
  all_types <- unique(unlist(strsplit((collisions$types),split=",")))  ###all types of interactions listed in collisions
  to_remove <- all_types [which(!all_types%in%to_include)]
  
  ##use the to_remove object to erase from collision$types the collisions we are not interested in, as well as the commas
  collisions_types <- collisions$types
  for (type in c(to_remove,",")){
    collisions_types <- gsub(type,"",collisions_types)
  }
  lines_to_remove <- which(nchar(collisions_types)==0)
  
  
  ##finally,remove lines with empty collisions$types
  collisions      <- collisions[-lines_to_remove,]
  
  return(collisions)
}

filter_ants         <- function(collisions,desired_ants,logical_link){
  if (logical_link == "or"){
    collisions <- collisions[which(collisions$ant1%in%desired_ants | collisions$ant2%in%desired_ants),]
  }
  if (logical_link == "and"){
    collisions <- collisions[which(collisions$ant1%in%desired_ants & collisions$ant2%in%desired_ants),]   
  }
 return(collisions) 
}

filter_distance     <- function(collisions,distance,logical_link){
  if (logical_link == "smaller"){
    collisions <- collisions[which(sqrt((collisions$x1-collisions$x2)^2+(collisions$y1-collisions$y2)^2)<=distance),]
  }
  if (logical_link == "greater"){
    collisions <- collisions[which(sqrt((collisions$x1-collisions$x2)^2+(collisions$y1-collisions$y2)^2)>=distance),]
  }
  return(collisions) 
}

filter_angle     <- function(collisions,angle_threshold,logical_link){
  ant_angle_diff  <-  (collisions$angle1-collisions$angle2) - (2*pi)*round((collisions$angle1-collisions$angle2)/(2*pi))
  angle_threshold <- (angle_threshold) - (2*pi)*round((angle_threshold)/(2*pi))
    
  if (logical_link == "smaller"){
    collisions <- collisions[which(abs(ant_angle_diff)<=angle_threshold),]
  }
  if (logical_link == "greater"){
    collisions <- collisions[which(abs(ant_angle_diff)>=angle_threshold),]  }
  return(collisions) 
}


interaction_detection <- function (e
                                   ,start
                                   ,end
                                   , max_time_gap
                                   , max_distance_moved
                                   , capsule_matcher=NULL
                                   , desired_ants_OR = NULL
                                   , desired_ants_AND = NULL
                                   , distance_smaller_than = NULL
                                   , distance_greater_than = NULL
                                   , angle_smaller_than = NULL
                                   , angle_greater_than = NULL
                                   , IF_frames 
){
  ###e: experiment
  ###time start: fmTime
  ###time_stop: fmTime
  ###capsule_matcher: list of vectors - each vector containing 2 capsule IDs that we want to see collide
  ###desired_ants_OR: vector of (any number) of ant IDs that we are interested in (any interaction involving at least one of these ants)
  ###desired_ants_AND: vector of (any number) of ant IDs that we are interested in (any interaction involving two ants from this list)
  ###distance_smaller_than: only collisions in which the distance between the two ants tag is less than distance_smaller_than will be kept
  ###distance_greater_than: only collisions in which the distance between the two ants tag is more than distance_greater_than will be kept
  ###angle_smaller_than: only collisions in which the angle between the two ants is less than angle_smaller_than will be kept
  ###angle_greater_than: only collisions in which the angle between the two ants is more than angle_greater_than will be kept
  
  
  ###to avoid crashes, perform queries by chunks of 12 hours
  absolute_end   <- end ###store the absolute end time so we know when to end
  
  ###initialise an object that will store all successive collisions
  all_collisions <-  NULL
  
  
  ###initialise argument for while
  extraction_complete <- F
  
  ###while extraction is not complete, keep querying for 12 hours segments
  while (!extraction_complete){
    ###get IF_frames for start and end
    frame_start <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    frame_end   <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    
    ###test if remaining time is superior to 12 hours
    if (as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) -   as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) > 13 * 60 * 60){
      ###if superior, define new end time as equal to start + 12 hours
      frame_end <- IF_frames[ match.closest(x = IF_frames[which(IF_frames$frame_num==frame_start),"time"]+13*60*60,table = as.numeric(IF_frames$time)),"frame_num"]
      end <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==frame_end) ,"time"   ])
    }else{
      ###otherwise, keep end time and declare extraction to be complete
      extraction_complete <- T
    }
    
    ###extract trajectory
    print(paste("Querying collisions between",capture.output(print(start)),"and",capture.output(print(end)),"..."))

    #######################################
    ###query collisions ###############
    #######################################
    collisions          <- fmQueryCollideFrames(e, start=start, end=end, showProgress = FALSE)
    
    #######################################
    ###filter on capsules #################
    #######################################
    if (!is.null(capsule_matcher)){
      collisions$collisions <- filter_capsules (collisions$collisions,capsule_matcher)
    }
    
    #######################################
    ###filter on ant ids #################
    #######################################
    if (!is.null(desired_ants_OR)){
      collisions$collisions <- filter_ants (collisions$collisions,desired_ants_OR,"or")
    }
    if (!is.null(desired_ants_AND)){
      collisions$collisions <- filter_ants (collisions$collisions,desired_ants_AND,"and")
    }
    
    #################################################################################################
    ###collpase collisions into single dataframe (necessary for next filtering step) ################
    #################################################################################################
    print("Collapsing collisions...")
    collisions <- collapse_collisions(collisions)
    print("Collisions collapsed.")
    
    #################################################################################################################################################
    ###remove duplicates between collisions and all_coliisions (as there are a few frames overlapping between successive queries ) ##################
    #################################################################################################################################################
    if (!is.null(all_collisions)){
      collisions <- collisions[which(as.numeric(collisions$time)>max(as.numeric(all_collisions$time),na.rm=T)),]
    }
    
    ##########################################################################################################################################
    ###apply correct offset to frames_row_index to correct frame number (automatically restarts at 1 in fmQueryCollideFrames) ################
    ##########################################################################################################################################
    collisions$frame_number <- IF_frames[match.closest(x = collisions[,"time"],table = as.numeric(IF_frames$time)),"frame_num"]
    
    
    all_collisions          <-  rbind(all_collisions,collisions)
    
    ###update start and end
    if (!extraction_complete){
      start <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==frame_end-5) ,"time"   ])          ### new start time = end time of last query - a few frames to ensure some overlap to avoid losing data (problems with rounding leading to mismatches between timestamps!)
      end   <- absolute_end ### new end time = absolute_end
    }
    
    ###clear memory
    rm(list=c("collisions"))
    gc()
    mallinfo::malloc.trim(0L)
  }
  
  collisions <- all_collisions; rm(list=c("all_collisions"));gc();mallinfo::malloc.trim(0L)
  #######################################
  ###filter on ant distances ############
  #######################################
  if (!is.null(distance_smaller_than)){
    collisions <- filter_distance (collisions,distance_smaller_than,"smaller")
  }
  if (!is.null(distance_greater_than)){
    collisions <- filter_distance (collisions,distance_greater_than,"greater")
  }
  
  #######################################
  ###filter on ant angles ###############
  #######################################
  if (!is.null(angle_smaller_than)){
    collisions <- filter_angle (collisions,angle_smaller_than,"smaller")
  }
  if (!is.null(angle_greater_than)){
    collisions <- filter_angle (collisions,angle_greater_than,"greater")
  }
  
  ###################################################################
  ###assemble successive collisions into interactions ###############
  ###################################################################
  ###prepare collisions to be loaded into cpp function
  collisions$time_second <- round(as.numeric(collisions$time),3)
 
   ##SLOW!!! :-( collisions$pair <- apply(collisions[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
  collisions <- within(collisions,pair <- paste(ant1,ant2,sep="_")) ###FAST :-)
  collisions$pair <- match(collisions$pair,sort(unique(collisions$pair)))-1 ###necessary for c++
  pair_list <- sort(unique(collisions$pair))
  
  print("About to start mergeing collisions into interaction list...")
  if (length(pair_list)>chunk_size_for_merge){ ###if too many pairs, merge interactions by chunks
    interactions <- NULL; 
    n_chunks <- ceiling(length(pair_list)/chunk_size_for_merge)
    for (chunk in 1:n_chunks){
      pair_list_chunk  <- pair_list [   ( 1 + (chunk_size_for_merge * (chunk-1)) ) : (min(chunk *chunk_size_for_merge, length(pair_list)     )  )] ###selects the next chunk_size_for_merge pairs of ants
      collisions_chunk <- collisions[which(collisions$pair %in%pair_list_chunk),]                                  ###subsets collision table to only contain those pairs
      ###the merge_interactions.cpp programs expects successive pair IDs starting from 0, so we need to update the content of both collisions_chunks and pair_list_chunk
      starting_pair    <- min(pair_list_chunk)                     ###first get the value of the first pair ID                                                 
      collisions_chunk$pair <- collisions_chunk$pair - starting_pair     ###subtract        starting_pair from       collisions_chunk$pair                           
      pair_list_chunk <- pair_list_chunk - starting_pair ###subtract        starting_pair from      pair_list_chunk            
      ####merge interactions for that chunk
      interactions <- rbind(interactions,merge_interactions(collisions_chunk, pair_list_chunk, max_distance_moved=max_distance_moved, max_time_gap=max_time_gap))
      ###clear memory 
      rm(list=c("pair_list_chunk","collisions_chunk","starting_pair"));mallinfo::malloc.trim(0L)
     }
  }else{
    interactions       <- merge_interactions(collisions, pair_list, max_distance_moved=max_distance_moved, max_time_gap=max_time_gap)
  }
  interactions       <- interactions[order(interactions$end,interactions$start,interactions$ant1,interactions$ant2),] ### reorder the table, necessary if interaction mergeing was done by chunks
  interactions$start <- as.POSIXct(interactions$start,  origin="1970-01-01", tz="GMT" )
  interactions$end   <- as.POSIXct(interactions$end,  origin="1970-01-01", tz="GMT" )

  return(interactions)
}

