
####################################################################################
#### THIS SCRIPT CONTAINS:
#### DATA MANIPULATION TO COVERT GROOMING INTERACTIONS FILE TO STROEYMEYT 2018 COMPATIBLE FORMAT 
#### the input files taken by this script could be the files produced by the grooming inference but her I provide a single file as produced by "plot_grooming_pre-post.R"
####################################################################################

#### LIBRARIES
gc()
mallinfo::malloc.trim(0L)

# starting params
# USER <- "2A13_Office" # Nath_office 
USER <- "2A13_Office" # Nath_office 

if (USER == "2A13_Office") {
  usr <- "cf19810"
} else {
  usr <- "bzniks"
}
### DIRECTORIES
WORKDIR <- paste("/media/",usr,"/DISK4/EXP1_base_analysis",sep="")
DATADIR <-  paste(WORKDIR,"Data/Adriano_extracted_groomings_final",sep="/")
METADATADIR <- paste(WORKDIR,"EXP_summary_data",sep="/")
INTDIR <- paste("/media/",usr,"/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment_grooming/intermediary_analysis_steps",sep="")


# PARAMETERS
N_DECIMALS <- 3 ## when assigning time in seconds, the number of decimals to preserve when rounding to match between interactions & collisions
TimeWind <- 3600 ## in seconds (3600 is an hour)

# Base files
metadata <- read.table(paste(METADATADIR, "/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
metadata_info         <- read.csv("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv", sep = ",")
inferred <- read.table(paste(WORKDIR,"/Data/inferred_groomings_ALL_FINAL.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
#check expected format
example <- read.table("/media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/intermediary_analysis_steps/full_interaction_lists/PostTreatment/observed/colony09SS_control.small_PostTreatment_interactions.txt",header=T,stringsAsFactors = F, sep="\t")
str(example)
str(inferred)
### assign time windows using information present in the file (check from grooming_script)

## define time windows
time_window_all        <- data.frame(
  PERIOD = "all"
  ,REPLICATE = metadata_info$REP_treat
  ,return_time = metadata_info$ReturnExposed_time #time_start
)
# 
time_window_all$return_time <- as.POSIXct(time_window_all$return_time, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
time_window_all$time_start <- time_window_all$return_time - 27*60*60 #(24h + 3h gap)
time_window_all$time_stop  <- time_window_all$return_time + 24*60*60
# 
Time_dictionary <- data.frame(time_hours = -36:35, time_of_day = rep(0:23, 3)) # time_hours= time since exposure (negatives pre-exposire, positive post exposure), time_of_day=time hour slot of the day
Time_dictionary <- Time_dictionary[which(Time_dictionary$time_hours <= 21 & Time_dictionary$time_hours >= -27), ]
Time_dictionary$period <- ifelse(Time_dictionary$time_hours<0, "pre", "post")

###### convert names etc as done in EXP1_base_analysis and Network script

# round digits
inferred$ant1.mean.x          <- round(as.numeric(inferred$ant1.mean.x),N_DECIMALS)
inferred$ant1.mean.y          <- round(as.numeric(inferred$ant1.mean.y),N_DECIMALS)
inferred$ant1.mean.angle      <- round(as.numeric(inferred$ant1.mean.angle),N_DECIMALS)
inferred$ant2.mean.x          <- round(as.numeric(inferred$ant2.mean.x),N_DECIMALS)
inferred$ant2.mean.y          <- round(as.numeric(inferred$ant2.mean.y),N_DECIMALS)
inferred$ant2.mean.angle      <- round(as.numeric(inferred$ant2.mean.angle),N_DECIMALS)

# Rename variables in a SCIENCE 2018 COMPATIBLE FORMAT
# EXPECTED FORMAT:: Tag1,Tag2,Startframe,Stopframe,Starttime,Stoptime,Box, Xcoor1,Ycoor1,Angle1,Xcoor2,Ycorr2,Angle2,Direction (types),Detections,
# DIRECTION → the capsule types matching
# DETECTION → the rate of frames in which both ants where present
names(inferred)[which(names(inferred) == "ant1")] <- "Tag1" # beware: these are antID but called Tags for compatibility
names(inferred)[which(names(inferred) == "ant2")] <- "Tag2"
names(inferred)[which(names(inferred) == "types")] <- "Direction"
names(inferred)[which(names(inferred) == "detections")] <- "Detections"
names(inferred)[which(names(inferred) == "startframe")] <- "Startframe"
names(inferred)[which(names(inferred) == "endframe")] <- "Stopframe"
names(inferred)[which(names(inferred) == "T_start_sec")] <- "Starttime"
names(inferred)[which(names(inferred) == "T_stop_sec")] <- "Stoptime"
names(inferred)[which(names(inferred) == "ant1.mean.x")] <- "Xcoor1"
names(inferred)[which(names(inferred) == "ant1.mean.y")] <- "Ycoor1"
names(inferred)[which(names(inferred) == "ant1.mean.angle")] <- "Angle1"
names(inferred)[which(names(inferred) == "ant2.mean.x")] <- "Xcoor2" # fixed typo in original naming Xcoor2
names(inferred)[which(names(inferred) == "ant2.mean.y")] <- "Ycoor2"
names(inferred)[which(names(inferred) == "ant2.mean.angle")] <- "Angle2"
names(inferred)[which(names(inferred) == "PERIOD")] <- "period"
names(inferred)[which(names(inferred) == "TREATMENT")] <- "treatment"

#inferred$Box <- e$spaces[[1]]$name # tracking system
# remove extras
inferred$frame_start    <- NULL
inferred$frame_stop     <- NULL
inferred$start          <- NULL
inferred$end            <- NULL
inferred$space          <- NULL
inferred$T_start_UNIX   <- NULL
inferred$T_stop_UNIX    <- NULL
inferred$ant1ID_str     <- NULL
inferred$ant2ID_str     <- NULL
inferred$REPLICATE      <- NULL
inferred$unique_interaction_id    <- NULL
inferred$predicted_Hit            <- NULL
inferred$ReturnExposed_time       <- NULL
inferred$ReturnExposed_time       <- NULL


# reassing var value
# Rename by name
inferred$treatment <- as.factor(inferred$treatment)
levels(inferred$treatment)[levels(inferred$treatment)=="BS"] <- "control.big"
levels(inferred$treatment)[levels(inferred$treatment)=="BP"] <- "pathogen.big"
levels(inferred$treatment)[levels(inferred$treatment)=="SS"] <- "control.small"
levels(inferred$treatment)[levels(inferred$treatment)=="SP"] <- "pathogen.small"

setdiff(names(example),names(inferred))
setdiff(names(inferred),names(example))


### save as individual files, completed with information
for (REP_TREAT in unique(inferred$REP_treat)) {
  
  #Interactions_REP_TREAT <- data.frame()
  
  # # # Extract the minimum "From" time and maximum "To" time for each PERIOD
  Period_windows <- data.frame(Period = c("pre","post")
                               , From = c(time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"time_start"], #skip 3h gap
                                          time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"]
                               ))
  Period_windows$To <-   Period_windows$From + 24*60*60
  
  
  
  ######################
  ### add  "time_hours"  "time_of_day" "colony" 
  #split by 3h blocks
  
  for (PERIOD in c("pre","post")) {
    
    #segment the time dictionary
    Time_dictionary_PERIOD <- Time_dictionary[which(Time_dictionary$period==PERIOD),]
    #conform naming to science2018
    # colony code
    REP_NUM           <- substring(REP_TREAT, 2, nchar(REP_TREAT))
    colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
    # colony_status
    status_char       <- substr(REP_TREAT, nchar(REP_TREAT), nchar(REP_TREAT))
    colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
    # treatment code
    period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
    size_char         <- substr(REP_TREAT, nchar(REP_TREAT)-1, nchar(REP_TREAT)-1)
    size_status       <- ifelse(size_char == "S", "small", ifelse(size_char == "B", "big", NA))
    treatment_code    <- paste(colony_status,size_status,sep=".")
    
    # Select the full PERIOD (24h)
    cat(paste("#######","period:", PERIOD, sep= " "))

      # DESTINATION FOLDER
      INTERACTIONS_FULL   <-  file.path(INTDIR,"full_interaction_lists",period_code,"observed", 
                                        paste(colony,treatment_code,period_code,"interactions.txt",sep="_"))
      
      #subset file by REP and PERIOD
      inferred_REP        <- inferred[which(inferred$REP_treat==REP_TREAT & inferred$period==PERIOD),]
      
      # Remove interactions involving dead ants
      dead_by_REP <- metadata[which(metadata$IsAlive==FALSE & metadata$REP_treat==REP_TREAT),"antID"]
      inferred_REP <- inferred_REP[!(inferred_REP$Tag1 %in% dead_by_REP | inferred_REP$Tag2 %in% dead_by_REP), ]
      
      # create time vars
      inferred_REP$time_hours   <- NA
      inferred_REP$time_of_day  <- NA
      # TIME_HOURS zero is the moment of exposed ants return
      for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours

        # TIME_OF_DAY
        TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
        #time windows
        From_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind),N_DECIMALS)
        To_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind),N_DECIMALS)
        
        #assing time labels
        inferred_REP[which(inferred_REP$Starttime >= From_TIME_HOURS & inferred_REP$Starttime <= To_TIME_HOURS),"time_hours"]   <- TIME_HOURS
        inferred_REP[which(inferred_REP$Starttime >= From_TIME_HOURS & inferred_REP$Starttime <= To_TIME_HOURS),"time_of_day"]  <- TIME_OF_DAY
                 
      }
      
      # Add metadata info
      inferred_REP <- cbind(data.frame(
        inferred_REP, colony=colony, #treatment=treatment_code,
        stringsAsFactors = F
      ))
      
      # extras not present in Science files: (added at the end of the file output)
      #"REP_treat","period","ant1.zones","ant2.zones","duration"
      inferred_REP <- inferred_REP[which(inferred_REP$time_hours!=-3), ##Remove Gap data
                                   c("Tag1","Tag2","pair","Act_Name","Rec_Name","Startframe","Stopframe","Starttime","Stoptime","Xcoor1","Ycoor1","Angle1","Xcoor2","Ycoor2","Angle2","Direction","Detections","time_hours","time_of_day","colony","treatment","REP_treat","period","ant1.zones","ant2.zones","duration","ReturnExposed_time_sec","start_time_sec","end_time_sec","preReturn_gap_time_sec")]

      ## Interactions save (saved INSIDE the Network_analysis folder)
      #if (file.exists(INTERACTIONS_FULL)) {
      #  write.table(Interactions, file = INTERACTIONS_FULL, append = T, col.names = F, row.names = F, quote = F, sep = ",")
      #} else {
      write.table(inferred_REP, file = INTERACTIONS_FULL, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
      #}
      
      ### split output into bins
      print(paste0("split files into 3-hours bins"))
      
      for (TIME_HOURS in unique(inferred_REP$time_hours)) {
        #labels for subsets
        TH <- paste0("TH",unique(inferred_REP[which(inferred_REP$time_hours==TIME_HOURS),"time_hours"]))
        TD <- paste0("TD",unique(inferred_REP[which(inferred_REP$time_hours==TIME_HOURS),"time_of_day"]))
        cat("\rTIME_HOURS", TH,"TIME_OF_DAY", TD,"PERIOD",PERIOD)
        
        INTERACTIONS_BINNED <-  file.path(INTDIR,"binned_interaction_lists",period_code,"observed", 
                                          paste(colony,treatment_code,period_code,TH,TD,"interactions.txt",sep="_"))
        #save object by TH and TD
        write.table(inferred_REP[which(inferred_REP$time_hours==TIME_HOURS),], file = INTERACTIONS_BINNED, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
      }

    
  } # PERIOD
}



