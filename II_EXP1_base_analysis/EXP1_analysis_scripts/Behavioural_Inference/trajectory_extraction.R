#### THIS SCRIPT IS USED INSIDE:
#### - AntTasks_v082.R
#### - SpaceUse_v082.R

merge_trajectory_segments <- function(positions,IF_frames){
  new_positions <- list(trajectories_summary=data.frame(),trajectories=list())

  print(paste("Mergeing trajectory segments for ",length(unique (positions$trajectories_summary$antID)),"ants..."))
  antID_str <-c()
  for (antID in sort(unique (positions$trajectories_summary$antID))){
    antID_str <- c(antID_str,paste("ant_",antID,sep="") )
    if (round(which(antID==sort(unique (positions$trajectories_summary$antID)))/10)==which(antID==sort(unique (positions$trajectories_summary$antID)))/10){
      print(paste("Mergeing trajectory segments for ",which(antID==sort(unique (positions$trajectories_summary$antID)))," out of ",length(unique (positions$trajectories_summary$antID)),"ants..."))
    }
    
    ant_indices <- which(positions$trajectories_summary$antID==antID)
    if (length(ant_indices)==1){ # LS: this "if statement" should not be commented out! -> otherwise ants with only one trajectory the time will be NA
      new_positions[["trajectories_summary"]]      <- rbind(new_positions[["trajectories_summary"]],positions$trajectories_summary[ant_indices,])
      new_positions[["trajectories"]] <- c(new_positions[["trajectories"]], list( positions$trajectories[[ant_indices]]      ))
    }else{
      ###the new summary object should contain only one line per ant - and that line should correspond to the first trajectory segment in line
      ant_summary       <- positions$trajectories_summary[ant_indices,]      ### extract all the lines corresponding to that ant in positions$trajectories_summary
      first_segment_idx <- which.min(ant_summary$start)                      ### find which one starts first 
      ant_summary$frame_num <- IF_frames[match.closest(x=ant_summary$start , table=as.numeric(IF_frames$time)),"frame_num"]
      new_positions[["trajectories_summary"]]      <- rbind(new_positions[["trajectories_summary"]],ant_summary[first_segment_idx,]) ### copy that line into new_positions[["trajectories_summary"]] 
      
      ###create a empty trajectory data.frame for that ant within the list
      new_positions[["trajectories"]] <- c(new_positions[["trajectories"]], list(data.frame()))

      ##now loop over the trajectory segments, in correct temporal order order
      for (ant_index in ant_indices[order(ant_summary$start)]){
        positions$trajectories[[ant_index]]$dt_diff                       <- c(NA,diff(positions$trajectories[[ant_index]]$time))
        positions$trajectories[[ant_index]]$frame_diff                    <- round(positions$trajectories[[ant_index]]$dt_diff*FRAME_RATE) 
        positions$trajectories[[ant_index]][1,"frame"]                    <- IF_frames[match.closest(x=ant_summary[which(ant_index==ant_indices[order(ant_summary$start)]),"start"] , table=as.numeric(IF_frames$time)),"frame_num"]
        positions$trajectories[[ant_index]][2:nrow(positions$trajectories[[ant_index]]),"frame"] <- positions$trajectories[[ant_index]][1,"frame"] + cumsum(positions$trajectories[[ant_index]][2:nrow(positions$trajectories[[ant_index]]),"frame_diff"])
        positions$trajectories[[ant_index]]$time                          <- as.numeric(IF_frames[positions$trajectories[[ant_index]]$frame,"time" ]) - ant_summary[first_segment_idx,"start"]
        new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]       <- rbind(new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]],positions$trajectories[[ant_index]])
        # free memory by removing the ant trajectory segment
        positions$trajectories[[ant_index]] <- data.frame()
        gc()
        mallinfo::malloc.trim(0L)
      }
      new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]    <- new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]   [which(!duplicated(new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]   $frame)),]
      new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]    <- new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]   [which(!names(new_positions[["trajectories"]] [[length(new_positions[["trajectories"]] )]]   )%in%c("dt_diff","frame_diff"))]
     }
    
    rm(list=ls()[which(ls()%in%c("ant_indices","ant_summary","first_segment_idx","ant_trajectory","ant_index"))])
    gc()
    mallinfo::malloc.trim(0L)
  }
  new_positions[["trajectories_summary"]]$antID_str <-  antID_str
  names(new_positions[["trajectories"]])            <-  antID_str
  print("All ant trajectories merged.")
  
  print("Clearing memory...")
  rm(list="positions")
  gc()
  mallinfo::malloc.trim(0L)
  Sys.sleep(5)
  
  print("About to return output from merge_trajectory_segments...")
  
  return(new_positions)
}

extract_trajectories <- function (e, start ,end , maximumGap ,computeZones=TRUE,showProgress = FALSE, IF_frames){
  ###to avoid crashes, perform queries by chunks of 12 hours
  absolute_end   <- end ###store the absolute end time so we know when to end
  
  ###initialise an object that will store all successive trajectory segments
  trajectories_summary <- data.frame()
  trajectories        <- list()
  
  ###initialise argument for while
  extraction_complete  <- F
  # first_extracted_frame <- NA
  # last_extracted_frame <-NA
  ###while extraction is not complete, keep querying for 12 hours segments
  while (!extraction_complete){
    ###get IF_frames for start and end
    frame_start <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    ##L: "Error in findInterval(x, table, rightmost.closed = FALSE, all.inside = TRUE) :'vec' must be sorted non-decreasingly and not contain NAs)  ####
    frame_end   <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    
    ###test if remaining time is superior to 12 hours
    if (as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) -   as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) > 6 * 60 * 60){
      ###if superior, define new end time as equal to start + 12 hours
      frame_end <- IF_frames[ match.closest(x = IF_frames[which(IF_frames$frame_num==frame_start),"time"]+6*60*60,table = as.numeric(IF_frames$time)),"frame_num"]
      end <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==frame_end) ,"time"   ])
    }else{
      ###otherwise, keep end time and declare extraction to be complete
      extraction_complete <- T
    }

    ###extract trajectory
    print(paste("Querying trajectories between",capture.output(print(start)),"and",capture.output(print(end)),"..."))
    positions <- fmQueryComputeAntTrajectories(e,start = start,end = end,maximumGap = maximumGap,computeZones = computeZones,showProgress = showProgress) #set true to obtain the zone of the ant
    
    ###get last extracted frame from the trajectories
    # print(paste("The last extracted frame in the previous segment was",last_extracted_frame))
    # first_extracted_frame   <-IF_frames[ match.closest( min(as.numeric(positions$trajectories_summary$start)),as.numeric(IF_frames$time)),"frame_num"]
    last_extracted_frame    <-  IF_frames[ match.closest( max(as.numeric(positions$trajectories_summary$start) + unlist(lapply(positions$trajectories,function(x)x[nrow(x),"time"]))),as.numeric(IF_frames$time)),"frame_num"]
    # print(paste("The first extracted frame in the new segment is",first_extracted_frame))
    ###add to overall objects 
    trajectories_summary <- rbind(trajectories_summary,positions$trajectories_summary)
    trajectories         <- c(trajectories,positions$trajectories)
    
    ###update start and end
    if (!extraction_complete){
      start <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==last_extracted_frame - 5 ) ,"time"   ]       ) ###extract from a few frames before to make sure no data is lost, as rounding times leads to uncertainties
      end   <- absolute_end ### new end time = absolute_end
     }
    
    ###clear memory
    rm(list=c("positions"))
    gc()
    mallinfo::malloc.trim(0L)
  }
  
  ###create new positions object that contains all the trajectory segments
  positions = list(trajectories_summary=trajectories_summary, "trajectories"=trajectories)
  
  ###clear memory 
  rm(list=c("trajectories_summary","trajectories"))
  gc()
  mallinfo::malloc.trim(0L)
  
  ###merge trajectory for each ant so we end up with a single trajectory per ant
  positions <- merge_trajectory_segments (positions,IF_frames)
  print("Output from merge_trajectory_segments received.")
  
  print("About to return output from extract_trajectories...")
  
  return(positions)
  
}