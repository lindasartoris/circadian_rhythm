
################## GET ZONES usage ###################################
SpaceUse <- function(e, start, end){
  
  #required packages
  require(lubridate)
  require(FortMyrmidon)
  require(mallinfo)
  #require(parallel)
  require(doParallel)
  require(data.table)
  
  # parameters (from Sceince2018)
  # warning("the parameters (from Sceince2018) may need moficiation, discuss with Nathalie" ) # L: Nathalie has approved the settings
  options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)
  max_time <- 1.5 # max time in sec, for which it is acceptable to calculate a distance moved (e.g. 1.5 sec: allows for up to 2 missing points between two successive detections and still calculates distance moved)
  min_segment_duration_sec <- 120; Lmin <- FRAME_RATE*min_segment_duration_sec #minimum duration of a bout, in frames
  max_gap <- 15 ##anything larger than that needs to be subdivided into two different bouts
  dist_threshold <- 25 # avg is 50px (from the fort experiment calculator) #previous value: 15, everything below that threshold (half a tag length) is considered inactive
  power <- 0.125 ###power used before clustering
  CORE_N <- 10
  
  ###  Functions
  insertRow <- function(existingDF, newrow, r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
  }
  
  #define zones
  #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
  zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
  zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
  
  
  foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
  
  #Get complete list of Ants
  AntID_list <- NULL
  for (ant in   1: length(e$ants)) {
    AntID_list <- c(AntID_list,e$ants[[ant]]$ID)
  }
  
  print("Querying frames...")
  IdentifyFrames      <- fmQueryIdentifyFrames(e,start, end,showProgress = FALSE)
  IF_frames           <- IdentifyFrames$frames
  # LS: was added because in some cases time was not sorted properly (only a few rows would be off)
  IF_frames           <- IF_frames[order(IF_frames$time), ] 
  
  rm(list=c("IdentifyFrames"))
  gc()
  mallinfo::malloc.trim(0L)
  
  # Assign a frame to each time since start and use it as baseline for all matching and computation
  IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))
  
  # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
  IF_frames$time_sec <- round(as.numeric(IF_frames$time),N_DECIMALS)
  IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
  IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)
  
  
  
  #   print(paste(" Time interval from hour",HOUR_start,"to",(HOUR_start-TimeWindow),"until e end",sep = " "))
  #positions                 <- fmQueryComputeAntTrajectories(e,start,end,maximumGap = fmHour(24*365),computeZones = TRUE,singleThreaded=FALSE)
  positions <- extract_trajectories(e, start , end, maximumGap = fmHour(24*365),computeZones = T,showProgress = F,IF_frames=IF_frames)
  positions_summaries       <- positions$trajectories_summary
  # positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
  positions_list            <- positions$trajectories
  
  
  # AW: warning("order the trajectories to ensure that the following calculations are correct" )
  # LS: addition to the script to check/sort trajectories because of the warning ####
  
  for (ant_idx in   1: length(e$ants)){
    
    ant <- paste("ant", e$ants[[ant_idx]]$ID , sep="_")

    # Check if the column is sorted in ascending order
    is_sorted <- !is.unsorted(positions_list[[ant]]$time) # omit NAs because otherwise !is.unsorted() does not work
    
    tryCatch({ # because not all ants from 1 to 1:length(e$ants) might exist in positions_list so "is_sorted" would neither be T or F and it gives an error # might no longer be necessary as we are using e$ants[[ant_idx]]$ID now
      if(is_sorted==FALSE)
      {warning(paste(REP_TREAT,"- trajectories for", ant, "were not sorted initially"))
        positions_list[[ant]] <- positions_list[[ant]] [ order(positions_list[[ant]]$time),]
      }else{
       # print(paste(ant,"trajectories are sorted"))
      }
    }, error = function(err){
      print(paste(ant, "did not appear in positions_list"))
    })
  }
  
  # loop through each dataframe and check for duplicated rows
  # LS: in my case some of the duplicated rows were just the repeated NAs 
  for (i in seq_along(positions_list)) {
    if (any(duplicated(positions_list[[i]]$time))) {
      print(i)
      remove_times <- positions_list[[i]][duplicated(positions_list[[i]]$time),"time"]
      positions_list[[i]] <- positions_list[[i]][which(!positions_list[[i]]$time %in% remove_times),]
      print(paste(length(remove_times),"duplicates found for antID ", names(positions_list)[i], sep=" "))
    }
  }
  
  
  ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
  if(!(nrow(positions_summaries)==length(positions_list)))stop("number of rows of positions_summaries is not equal to the length of positions_list")
  ####2/ always make sure that positions_summaries is ordered correctly, using the index column
  # positions_summaries <- positions_summaries[order(positions_summaries$index),]
  ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
  
  ### parallelise process
  # register a parallel backend with x cores
  registerDoParallel(CORE_N)
  
  to_keep <- unique(c(to_keep,ls(),"ant_index")) # "traj_file","to_keep_ori"
  
  SpaceUse_loop_start_time <- Sys.time()
  
  
  # parallel loop to calculate summary statistics for each dataframe in the list
  results <- foreach(antID = names(positions_list)) %dopar% { #, .combine = "list" # LS: changed to "antID = names(positions_list)" to keep identities instead of iterating by ant index (you should never do that)
    print(paste("currently computing for ant in antID", antID))
    if (length(positions_list[[antID]]$time)>3){ # LS: added this bc the calculations go wrong when there are only 3 or fewer trajectories
      
      #############################################  
      ####### analyse trajectory ANT BY ANT ####### 
      #SpaceUse_ant <- function(positions_list_ant_index) {
      #for ( ant_index in 1:length(positions_list)) {
      #print(paste0("ant_index: ",ant_index))
      #cat("\rant_index:", ant_index, "of", length(positions_list))
      ####8_process_trajectory_files.R#####
      # Modified from Stroeymeyt et al. 2018
      #### Defines activity bouts and calculates summary statistics for each trajectory 
      #### Outputs: modified trajectory file and modifed individual behaviour file
      
      ####get trajectory list #####
      #input_traj <- paste(data_path,"intermediary_analysis_steps/trajectories",sep="/")
      #setwd(input_traj)
      #traj_list <- positions_list[[ant_index]]
      #### reduce traj list to only PreTreatment folders
      #### for each ant we will successively consider the Pre-treatment, then the Post-treatment trajectory
      #traj_list <- traj_list[grepl("Pre",traj_list)]
      
      traj <- NULL
      
      ##################################
      #for (traj_folder in traj_list){
      # root_name <- unlist(strsplit(traj_folder,split="/"))[grepl("colony",unlist(strsplit(traj_folder,split="/")))]
      # print(paste(unlist(strsplit(root_name,split="_"))[!grepl("Treatment",unlist(strsplit(root_name,split="_")))],collapse="_"))
      # colony    <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
      # tagfile   <- tag_list[grepl(colony,tag_list)]
      # splitfiles   <- split_list[grepl(colony,split_list)]
      # 
      # ####read-in tag list to define list of live ants
      # tag        <- read.tag(tagfile)$tag; names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
      # alive      <- paste("ant_",tag[which(tag$final_status=="alive"),"tag"],".txt",sep="")
      # 
      # splitinfo_Pre <- read.table(splitfiles[grepl("Pre",splitfiles)],header=T, stringsAsFactors = F)
      # splitinfo_Post <- read.table(splitfiles[grepl("Post",splitfiles)],header=T, stringsAsFactors = F)
      # 
      # ###Navigate in folder
      # setwd(traj_folder)
      # traj_files <- list.files()
      # 
      # to_keep <- unique(c(to_keep,ls(),"traj_file","to_keep_ori"))
      ###Process each trajectory file
      #for (traj_file in traj_files){
      #print(gsub("_"," ",gsub("\\.txt","",traj_file)))
      #traj_file_Pre  <- paste(traj_folder,traj_file,sep="/")
      #traj_file_Post <- gsub("Pre","Post",traj_file_Pre)
      ####delete files corresponding to dead ants
      # if (!traj_file %in% alive){
      #   print("Deleting file (ant died before the end)")
      #   file.remove(traj_file_Pre)
      #   file.remove(traj_file_Post)
      # }else{# if (!traj_file %in% alive)
      # print("Modifying trajectory...")
      ###check if analysis has already been done
      #output <- system(paste( "head -1 ",traj_file_Post),intern=T)
      #if (!grepl("bout_id",output)){
      bout_counter <- 0
      ####do the analysis for each of the Pre- and Post- files 
      #for (suffix in c("_Pre","_Post")){
      
      #output <- system(paste( "head -1 ",get(paste("traj_file",suffix,sep=""))),intern=T)
      #if (!grepl("bout_id",output)){
      #warning("positions_list_ant_index also afterwards in the function!")
      traj <- positions_list[[antID]] #read.table(get(paste("traj_file",suffix,sep="")),header=T,sep=",",comment.char="%")
      #traj <- positions_list_ant_index
      #names(traj)[names(traj)=="X.frame"] <- "frame"
      #####Remove occasional duplicates
      #traj <- traj[!duplicated(traj$frame),]
      ##########################################################      
      ######1. calculate distance moved ########################
      ##########################################################
      ####Prepare two tables offset by one line to facilitate calculation for calculations
      traj_1 <- insertRow(traj,rep(NA,ncol(traj)),1)
      traj_2 <- insertRow(traj,rep(NA,ncol(traj)),(nrow(traj)+1))[1:(nrow(traj)+1),]
      names(traj_1) <- paste(names(traj_1),"_previousframe",sep="")
      ####merge them
      full_trajectory <- data.frame(traj_1,traj_2)
      ####calculate time difference between two successive coordinates
      full_trajectory$time_diff <- full_trajectory$time-full_trajectory$time_previousframe
      ####round time_diff to nearest (1/FRAME_RATE)
      full_trajectory$time_diff <- round(full_trajectory$time_diff/(1/FRAME_RATE))*(1/FRAME_RATE);full_trajectory$time_diff[full_trajectory$time_diff==0]<-(1/FRAME_RATE)
      #### calculate distance moved between two successive frames
      full_trajectory$total_dist <- (sqrt((full_trajectory$x-full_trajectory$x_previousframe)^2+(full_trajectory$y-full_trajectory$y_previousframe)^2))
      #### calculate a distance per frame (i.e. speed) for intervals that are less than max_time
      #AW: replaced as too slow
      #full_trajectory[which(full_trajectory$time_diff<=max_time),"dist_per_frame"] <- (sqrt((full_trajectory[which(full_trajectory$time_diff<=max_time),"x"]-full_trajectory[which(full_trajectory$time_diff<=max_time),"x_previousframe"])^2+(full_trajectory[which(full_trajectory$time_diff<=max_time),"y"]-full_trajectory[which(full_trajectory$time_diff<=max_time),"y_previousframe"])^2))/((full_trajectory[which(full_trajectory$time_diff<=max_time),"time_diff"]/(1/FRAME_RATE)))
      # Create logical vector to subset data
      subset_index <- full_trajectory$time_diff <= max_time
      #first row can't be calculated
      subset_index[is.na(subset_index)] <- FALSE
      # Calculate distance per frame using logical vector
      full_trajectory$dist_per_frame[subset_index] <- sqrt((full_trajectory$x[subset_index] - full_trajectory$x_previousframe[subset_index])^2 + (full_trajectory$y[subset_index] - full_trajectory$y_previousframe[subset_index])^2) / (full_trajectory$time_diff[subset_index] / (1/FRAME_RATE))
      #### round
      full_trajectory$dist_per_frame <- round(100*full_trajectory$dist_per_frame)/100
      full_trajectory$total_dist <- round(100*full_trajectory$total_dist)/100
      #### now use this new data to create a new traj file
      traj <- full_trajectory[(1:(nrow(full_trajectory)-1)),c("time","x","y","zone","dist_per_frame","total_dist")] #"frame"
      
      #####calculate turn angles
      trajectory <- as.ltraj(traj[c("x","y")],date=as.POSIXct(traj$time,origin="1970-01-01"),id=gsub("ant_","",antID) , typeII=T,slsp="missing") # ,id=gsub("\\.txt","",gsub("ant_","",traj_file))
      turn_angles <- abs(trajectory[[1]]$rel.angle)
      #####modify it slightly, because I defined distance as distance moved between current frame and the previous; whereas turn angles are defined between current frame and the next
      turn_angles <- c(NA,turn_angles[1:(length(turn_angles)-1)])
      traj$turn_angle <- turn_angles
      #####no movement = no turn. Therefore fill in NA turn angles when dist=0
      traj[which(is.na(traj$turn_angle)&traj$dist_per_frame==0),"turn_angle"] <- 0
      #####whenever the time separating 2 successive detections is greater than max_time, enter NA; because then the turn angle is not meaningful
      times_1 <- c(NA,traj$time);times_2 <- c(traj$time,NA);time_diff <- times_2-times_1;time_diff <- time_diff[1:(length(time_diff)-1)]
      time_diff <- round(time_diff/(1/FRAME_RATE))*(1/FRAME_RATE)
      traj[which(time_diff>max_time),"turn_angle"] <- NA
      # }else{
      #   traj <- read.table(get(paste("traj_file",suffix,sep="")),header=T,comment.char="%")
      #   traj <- traj[,which(!names(traj)%in%c("type","bout_id"))]
      # }
      ####################################################################    
      #######2. cut trajectory into bouts of activity vs. inactivity #####
      ####################################################################      
      
      #####define an object that does not contain NAs in the dist_per_frame column
      traj_noNA <- traj[which(!is.na(traj$dist_per_frame)),]
      
      #####the analysis is only possible if there is enough data in traj_noNA
      if (nrow(traj_noNA)>2*Lmin){
        #####find segments based on distance moved
        #####apply distance threshold
        traj_noNA[traj_noNA$dist_per_frame<= dist_threshold,"dist_per_frame"] <- 0
        
        #####use function cpt.meanvar to find breakpoints in the trajectory data
        if (length(unique(traj_noNA$dist_per_frame))>2){
          segments_dist <- changepoint::cpt.meanvar(traj_noNA$dist_per_frame, method = "PELT",minseglen=Lmin)#####best by far
          breakpoints <- segments_dist@cpts
          rm(list="segments_dist")
        }else{
          breakpoints <- nrow(traj_noNA)
        }
        #####use breakpoints to define start and end times for bouts, and add the information into thr traj file 
        breakpoints <- match(traj_noNA[breakpoints,"time"],traj$time)
        breakpoints[length(breakpoints)] <- max(nrow(traj),breakpoints[length(breakpoints)])
        if (min(breakpoints)==nrow(traj)) { #AW # if breakpoints has only one element, it is th minimum value
          bout_indices <- data.frame(start=1,end=nrow(traj)) #issue in original code: if no breakpoint, create single bout
        }else{
          #breakpoints <-  breakpoints[!is.na(breakpoints)]#AW
          breakpoints <- c(1,breakpoints)
          first_index <- breakpoints[1]
          bout_start_indices <- c(first_index,(breakpoints+1)[2:(length(breakpoints)-1)])
          bout_end_indices <- breakpoints[2:(length(breakpoints))]
          bout_indices <- data.frame(start=bout_start_indices,end=bout_end_indices)
        } #AW
        bout_count_index <- NULL
        for (i in 1:nrow(bout_indices)){
          bout_count_index <- c(  bout_count_index , rep(i,(bout_indices[i,"end"]-bout_indices[i,"start"]+1)))
        }
        bout_count_index <- bout_counter + bout_count_index
        traj$bout_index <- bout_count_index
        bout_counter <- max(bout_count_index,na.rm=T)
      }#if (nrow(traj_noNA)>2*Lmin)
      #assign(paste("traj",suffix,sep=""),traj)
      #rm(list=c("traj"))
      #}##for (suffix in c("_Pre","_Post"))
      ###Now combine the two traj objects into a single one for the meta-analysis of active vs. inactive
      #if (ncol(traj_Pre)==ncol(traj_Post)){###do it only if bouts were defined for both periods
      if (bout_counter>2){###do it only if more than 1 bout is defined
        #trajectory_table <- rbind(data.frame(period="before",traj_Pre),data.frame(period="after",traj_Post))
        #TEMP
        #trajectory_table <- traj
        ###Perform cluster analysis on the whole table to determine active/inactive using the mean and standard deviation of speed (dist_per_frame) and turn_angle as an input
        for_multi <- aggregate(na.action="na.pass",cbind(turn_angle,(dist_per_frame)^power)~bout_index,function(x)cbind(mean(x,na.rm=T),sd(x,na.rm=T)),data=traj)
        for_multi <- data.frame(for_multi$bout_index,for_multi$turn_angle,for_multi$V2);names(for_multi) <- c("bout_index","turn_angle","turn_angle_SD","Dist","Dist_SD")
        cluster_BOTH <- kmeans( na.omit( for_multi   [c("turn_angle","Dist","turn_angle_SD","Dist_SD")]),  centers=2) ## clustering on ave & sd of relative turn angle & average distance
        ####Distinguish between active and inactive bout using the speed (active corresponds to higher speed)
        types_BOTH <- c("inactive","active")[order(data.frame(cluster_BOTH["centers"])$centers.Dist)]
        traj$type <- types_BOTH[cluster_BOTH$cluster[traj$bout_index]]
        
        ####Use the results from the cluster analysis to define bouts and interbouts
        if (length(unique(traj$type))==1){
          if (unique(traj$type)=="active"){to_fill <- "bout1"}else{to_fill <- "interbout1"}
          traj[1:nrow(traj),"bout_id"] <- to_fill
          
        }else{
          ##find indices of changes in activity type, in order to pool successive bouts of the same type together
          changes <- data.frame(type=c(NA,traj$type),type.2=c(traj$type,NA))
          changes$change <- changes$type==changes$type.2
          breakpoints <- which(!changes$change[1:(length(changes$change)-1)])
          
          start_indices <- c(1,breakpoints)
          end_indices <- c(breakpoints-1,nrow(traj))
          
          if(traj[1,"type"]!="active"){
            interbout_indices <- 1+2*c(0:(ceiling(length(start_indices)/2)-1))
            bout_indices <- 2*c(1:floor(length(start_indices)/2))
          }else{
            bout_indices <- 1+2*c(0:(ceiling(length(start_indices)/2)-1))
            interbout_indices <- 2*c(1:floor(length(start_indices)/2)) 
          }
          
          bout_indices <- data.frame(start=start_indices[bout_indices],end=end_indices[bout_indices])
          interbout_indices <- data.frame(start=start_indices[interbout_indices],end=end_indices[interbout_indices])
          
          ###copy the new bout information into the traj object
          bout_index_list <- NULL
          bout_count_index <- NULL
          for (i in 1:nrow(bout_indices )){
            bout_index_list <- c(bout_index_list, seq(bout_indices[i,"start"],bout_indices[i,"end"]))
            bout_count_index <- c(  bout_count_index , rep(i,(bout_indices[i,"end"]-bout_indices[i,"start"]+1)))
          }
          interbout_index_list <- NULL
          interbout_count_index <- NULL
          for (i in 1:nrow(interbout_indices )){
            interbout_index_list <- c(interbout_index_list, seq(interbout_indices[i,"start"],interbout_indices[i,"end"]))
            interbout_count_index <- c(  interbout_count_index , rep(i,(interbout_indices[i,"end"]-interbout_indices[i,"start"]+1)))
          }
          
          traj[bout_index_list,"bout_id"] <- paste("bout",bout_count_index,sep="")
          traj[interbout_index_list,"bout_id"] <- paste("interbout",interbout_count_index,sep="")
        }
        ####Finally, find large gaps in trajectory and subdivide each bout that was defined across the gap
        ####separate bout_id into root and index
        temp <- gregexpr("[0-9]+", traj$bout_id)
        traj$idx <- as.numeric(unlist(regmatches(traj$bout_id,temp)))
        traj[paste("bout",traj$idx,sep="")==traj$bout_id,"root"] <- "bout"
        traj[paste("interbout",traj$idx,sep="")==traj$bout_id,"root"] <- "interbout"
        
        ###get times where more than 15 min gap
        times_1 <- c(NA,traj$time);times_2 <- c(traj$time,NA);time_diff <- times_2-times_1;time_diff <- time_diff[1:(length(time_diff)-1)]
        indices <- which(time_diff>(60*max_gap))
        for (index in indices){
          traj[(index:nrow(traj)),][traj[(index:nrow(traj)),"root"]==traj[index,"root"],"idx"] <-       traj[(index:nrow(traj)),][traj[(index:nrow(traj)),"root"]==traj[index,"root"],"idx"] + 1
        }
        
        ###copy new bout id into the table
        traj <- within(traj, bout_id <- paste(root,idx,sep=""))
        
        # ###Finally, divide the trajectory:table into pre- and post-treatment trajectories
        # traj_Pre <-trajectory_table[which(trajectory_table$period=="before"),]
        # traj_Post <-trajectory_table[which(trajectory_table$period=="after"),]
      }else{
        #   ###In case one of the periods did not have enough data for bouts to be defined, fill the columns with NA 
        #   traj_Pre  <- data.frame(period="before",traj_Pre); traj_Pre["bout_id"] <- NA; traj_Pre["type"] <- NA
        #   traj_Post <- data.frame(period="after",traj_Post); traj_Post["bout_id"] <- NA; traj_Post["type"] <- NA
        traj$bout_id <- NA; traj$type <- NA #AW
      }
      #for (suffix in c("_Pre","_Post")){
      # ###Add metadata and for each of the pre-treatment and pre-treatment trajectories
      # splitinfo           <- get(paste("splitinfo",suffix,sep=""))
      # traj                <- get(paste("traj",suffix,sep=""))
      # indices             <- unlist(lapply(traj$time,function(x)max(which(x>=splitinfo$time))))
      # traj["time_hours"]  <- splitinfo$time_hours[indices]
      # traj["time_of_day"] <- splitinfo$time_of_day[indices]
      # traj["colony"]      <- colony
      # traj["treatment"]   <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"treatment"]
      #traj                <- traj[c("time","frame","x","y","dist_per_frame","total_dist","turn_angle","type","bout_id")] #"colony","treatment","box","period","time_hours","time_of_day",
      traj                <- traj[c("time","x","y","zone","dist_per_frame","total_dist","turn_angle","type","bout_id")] #"colony","treatment","box","period","time_hours","time_of_day",   "frame",
      ##################################
      ####### Now analyse trajectory ###
      ##################################
      # ####prepare output
      # tab <- expand.grid(colony=colony,
      #                    tag=gsub("\\.txt","",gsub("ant_","",traj_file)),
      #                    time_hours=splitinfo$time_hours,
      #                    proportion_time_active=NA,
      #                    average_bout_speed_pixpersec=NA,
      #                    total_distance_travelled_pix=NA)
      ###Analyse separately for each time point
      #for (time_point in unique(tab$time_hours)){
      #subtraj <- traj[which(traj$time_hours==time_point),]
      
      # calculate summary statistics
      tab <- data.frame(
        antID_str=antID,
        nb_frames_outside=NA ,
        nb_frames_inside=NA ,
        prop_time_outside=NA ,
        proportion_time_active=NA ,
        average_bout_speed_pixpersec=NA ,
        total_distance_travelled_pix=NA 
      )
      
      if (!all(is.na(traj$type))){
        bouts <- traj[which(traj$type=="active"),]
        interbouts <- traj[which(traj$type=="inactive"),]
        
        bout_number <- length(unique(bouts$bout_id))
        interbout_number <- length(unique(interbouts$bout_id))
        if (!(bout_number==0&interbout_number==0)){###if there is no bout data, leave everything as NA
          if (bout_number==0){###if ant completely inactive
            #positions_summaries[ant_index,c("proportion_time_active","total_distance_travelled_pix")] <- 0
            tab$proportion_time_active <- 0;  tab$total_distance_travelled_pix <- 0
            # tab[which(tab$time_hours==time_point),c("proportion_time_active","total_distance_travelled_pix")] <- 0
          }else {###if ant shows at least some activity
            # tab[which(tab$time_hours==time_point),"total_distance_travelled_pix"] <- sum(bouts$total_dist,na.rm=T)
            #positions_summaries[ant_index,"total_distance_travelled_pix"] <- sum(bouts$total_dist,na.rm=T)
            tab$total_distance_travelled_pix <- sum(bouts$total_dist,na.rm=T)
            bout_speed <- aggregate(na.rm=T,na.action="na.pass",dist_per_frame~bout_id,FUN=mean,data=bouts)
            #positions_summaries[ant_index,"average_bout_speed_pixpersec"] <- mean(bout_speed$dist_per_frame,na.rm = T)*FRAME_RATE ####multiply by 2 because each frame lasts 0.5sec
            tab$average_bout_speed_pixpersec <- mean(bout_speed$dist_per_frame,na.rm = T)*FRAME_RATE ####multiply by 2 because each frame lasts 0.5sec
            #tab[which(tab$time_hours==time_point),"average_bout_speed_pixpersec"] <- mean(bout_speed$dist_per_frame,na.rm = T)*2 ####multiply by 2 because each frame lasts 0.5sec
            if (interbout_number==0){###if ant completely active
              #positions_summaries[ant_index,"proportion_time_active"] <- 1
              tab$proportion_time_active <- 1
              # tab[which(tab$time_hours==time_point),"proportion_time_active"] <- 1
            }else{###if ant shows both activity and inactivity
              ###calculate cumulated bout duration
              bout_starts <- aggregate(na.rm=T,na.action="na.pass",time~bout_id,FUN=min,data=bouts); names(bout_starts)[names(bout_starts)=="time"] <- "Starttime"
              bout_ends <- aggregate(na.rm=T,na.action="na.pass",time~bout_id,FUN=max,data=bouts); names(bout_ends)[names(bout_ends)=="time"] <- "Stoptime"
              bout_durations <- merge(bout_starts,bout_ends)
              bout_duration <-  sum(bout_durations$Stoptime-bout_durations$Starttime+(1/FRAME_RATE),na.rm=T)
              ###calculate cumulated interbout duration                    
              interbout_starts <- aggregate(na.rm=T,na.action="na.pass",time~bout_id,FUN=min,data=interbouts); names(interbout_starts)[names(interbout_starts)=="time"] <- "Starttime"
              interbout_ends <- aggregate(na.rm=T,na.action="na.pass",time~bout_id,FUN=max,data=interbouts); names(interbout_ends)[names(interbout_ends)=="time"] <- "Stoptime"
              interbout_durations <- merge(interbout_starts,interbout_ends)
              interbout_duration <-  sum(interbout_durations$Stoptime-interbout_durations$Starttime+(1/FRAME_RATE),na.rm=T)
              
              #positions_summaries[ant_index,"proportion_time_active"] <- bout_duration/(bout_duration+interbout_duration)
              tab$proportion_time_active <- bout_duration/(bout_duration+interbout_duration)
              #tab[which(tab$time_hours==time_point),"proportion_time_active"] <- bout_duration/(bout_duration+interbout_duration)
            }
          }
        }
      }
      #}
      
      #TAB_ALL <- rbind(TAB_ALL,tab)
      ###Finally, write results
      ###individual behaviour
      # behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
      # if (!"proportion_time_active"%in%names(behav)){
      #   behav[c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")] <- NA
      # }
      # behav[match(as.character(interaction(tab$colony,tab$tag,tab$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")]  <- tab[c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")]
      # options(digits=3)
      # write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
      # options(digits=16)
      # write.table(traj,file=get(paste("traj_file",suffix,sep="")), row.names=F, col.names=T,append=F,quote=F)
      #}
      # }else{#if (!grepl("bout_id",output))
      #   print("Ant already processed.")
      # }
      #}#else (from if (!traj_file %in% alive))
      #clean()
      #}#for (traj_file in traj_files)
      #}# for (traj_folder in traj_list)
      #to_keep <- to_keep_ori
      #--------------------------------------------------------------- 
      
      # OUTPUTS
      # positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
      # positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
      # positions_summaries[ant_index,"prop_time_outside"] <- positions_summaries[ant_index,"nb_frames_outside"] /(positions_summaries[ant_index,"nb_frames_outside"] + positions_summaries[ant_index,"nb_frames_inside"])
      tab$nb_frames_outside <- length(which(traj[,"zone"]%in%foraging_zone))
      tab$nb_frames_inside  <- length(which(traj[,"zone"]%in%nest_zone))
      tab$prop_time_outside <- tab$nb_frames_outside /(tab$nb_frames_outside + tab$nb_frames_inside)
      
      # save summary statistics as new columns in the summary dataframe
      #positions_summaries[ant_index, names(tab)] <- as.vector(unlist(tab))
      
      
      #may conflict with parallelization?
      #rm(list=ls()[which(!ls()%in%c(to_keep))]) #close experiment
      #}
    }else{ # LS: <= 3 trajectories
      # warning(paste("in replicate",colony, "ant with ID", antID, "has fewer than 3 trajectories in positions_list"))
      
      tab <- data.frame(
        antID_str=antID,
        nb_frames_outside=NA ,
        nb_frames_inside=NA ,
        prop_time_outside=NA ,
        proportion_time_active=NA ,
        average_bout_speed_pixpersec=NA ,
        total_distance_travelled_pix=NA
      )

      tab$antID_str <- antID
      tab$nb_frames_outside <- NA
      tab$nb_frames_inside <- NA
      tab$prop_time_outside <- NA
      tab$proportion_time_active <- NA
      tab$average_bout_speed_pixpersec <- NA
      tab$total_distance_travelled_pix <- NA
    
    }
    return(list(antID_str=antID, data=tab))
    
  } # function SpaceUse_ant
  
  SpaceUse_loop_end_time <- Sys.time()
  print(paste("spaceUse 3h chunk took ", round(as.numeric(difftime(SpaceUse_loop_end_time, SpaceUse_loop_start_time, units = "secs")),1), " sec to complete"))
  
  # stop the parallel backend
  stopImplicitCluster()
  
  # convert the results to a list by individual index
  results_list <- lapply(seq_along(results), function(i) {
    # extract the result corresponding to the current individual index
    result <- results[[i]]
    if (!is.null(result)) {
      # return the summary statistics as a data frame
      return(result$data)
    } else {
      # return an empty data frame
      return(data.frame())
    }
  })
 
  # combine the results into a single data frame
  summary_df <- do.call(rbind, results_list)
  positions_summaries <- merge(positions_summaries, summary_df,by="antID_str",all.x=T) 
  #positions_summaries[100:(length(positions_summaries$index)),]
  #match antID and tagID (live tracking gives tagID). 
  IDs <- e$identificationsAt(fmTimeCreate(offset=fmQueryGetDataInformations(e)$start))
  IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
  IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
  positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
  
  # SpaceUsage <- aggregate(cbind(nb_frames_outside,
  #                               nb_frames_inside,
  #                               prop_time_outside,
  #                               proportion_time_active,
  #                               average_bout_speed_pixpersec,
  #                               total_distance_travelled_pix
  #                               ) ~ antID + tag_hex_ID, FUN = sum,na.action=na.pass,positions_summaries ) #, na.rm=T
  
  # LS: the line below previously removed all the ants that did not move so they were treated as missing even though they were detected # fixed
  SpaceUsage <- aggregate(na.rm=T,na.action="na.pass",cbind(nb_frames_outside, nb_frames_inside, total_distance_travelled_pix) ~ antID, data = positions_summaries, sum)
  # LS: use merge instead of cbind here (safer)
  SpaceUsage <- merge(SpaceUsage, aggregate(na.rm=T,na.action="na.pass",cbind(proportion_time_active, average_bout_speed_pixpersec) ~ antID, data = positions_summaries, mean))
  
  #it does not make sense to calculate the prop_time_outside beforehand as, if you have a doubled id, it will be messed up. 
  SpaceUsage$prop_time_outside <- SpaceUsage$nb_frames_outside /(SpaceUsage$nb_frames_outside + SpaceUsage$nb_frames_inside)
  
  #order output
  SpaceUsage <- SpaceUsage[, c("antID","nb_frames_outside","nb_frames_inside","prop_time_outside","proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")]
  
  #round
  SpaceUsage[, sapply(SpaceUsage, is.numeric)] <- round(SpaceUsage[, sapply(SpaceUsage, is.numeric)], 3)
  
  # rm(list=ls()[which(!ls()%in%c("SpaceUsage","e","foraging_zone","nest_zone","AntID_list","hour_chunk_start","TimeWindow","loop_N"))]) #close experiment
  # gc()
  # mallinfo::malloc.trim(0L)
  # 
  
  # #add missing ants
  missing_ants <- subset(AntID_list, !(AntID_list %in% SpaceUsage$antID))
  missing_ants_table <- data.frame()
  
  for (MISSING in missing_ants) {
    for (id in length(e$ants[[MISSING]]$identifications)) {
      #print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
      missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING)) # , tag_hex_ID= e$ants[[MISSING]]$identifications[[id]]$tagValue
    }}
  
  if (nrow(missing_ants_table)>0) {
    #add empty cols
    missing_ants_table[setdiff(names(SpaceUsage),names(missing_ants_table))] <- NA
    SpaceUsage <- rbind(SpaceUsage, missing_ants_table)
  }
  SpaceUsage <- SpaceUsage[order(SpaceUsage$antID),]
  #SpaceUsage$tag_hex_ID <- NULL
  
  rm(list=ls()[which(!ls()%in%c("SpaceUsage"))]) #close experiment
  gc()
  mallinfo::malloc.trim(0L)
  
  #print("SpaceUse computed")
  ##RETURN OUTPUT
  return(SpaceUsage)
  
}
