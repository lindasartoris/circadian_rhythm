#to clone the grooming capsules Def3 onto the newest myrmidon files (incl. metadata and re-orientation wherever necessary)
#remove Def2018 and add Def3

rm(list=ls())
gc()

########################################################################
################# ORIENT ALL FILES (AntsCreated) #######################
# Once DataPrep1_Clone-capsule-manual-to-manual_v082.R is run, check in fort-studio that capsules are present and normal and that Queen infos are overwritten (tag size, manual orientation, manual capsules).
# In auto_orientation_loop.R, load these 5 oriented and capsule provided files and orient all of the other files

"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(data.table)
library(circular)
library(R.utils)
library(stringr)

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

create_metadata <- function (dir_data){
  metadata <- data.frame(file_name=list.dirs(path = dir_data,recursive=F,full.names = F),tracking_system=unlist(lapply (   list.dirs(path = dir_data,recursive=F,full.names = F)    ,  function(x)unlist(strsplit (x,split="_"))[1])))
  return(metadata)
  
}

###source C++ movement direction program
sourceCpp("~/Documents/R_scripts/Get_Movement_Angle.cpp")

### directory of data and myrmidon files
dir_data <- '/media/lsartori/LS_1/circadian_rhythm_2022_experiment/tracking'
metadata <- create_metadata(dir_data)

#### ACCESS FILES
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(dir_data, n = 1)

### manually oriented ref file name
manually_oriented_files_CapsuleDef3 <- list.files(path=dir_data, pattern="manual_orientation_CapsuleDef3", recursive=T)
#ref_orient_caps_file <-  manually_oriented_files_CapsuleDef3

#FILES TO WHICH ASSIGN CAPSULE
ToOrient_file <-  list.files(path=dir_data, pattern="automatically_oriented_CapsuleDef2018_metadata_deaths_q_zones.myrmidon", recursive=T)



### Loop through all the directories in the dir_folder

# loop through the unique Tracking systems

to_keep_1 <- c(ls(),"to_keep_1","TS","CAPSULEDEF","data_list")
#unique(ref_orient_caps_file$TrackSys_name)

#to test the loop
#manually_oriented_files_text <- c("manually_oriented_files_CapsuleDef3")
#tracking_system <- "prideaux" 
#myrmidon_file <- as.character(data_list)

for ( manually_oriented_files_text in c("manually_oriented_files_CapsuleDef3")){
  manually_oriented_files <- get(manually_oriented_files_text)
  
  for (tracking_system in unique(metadata$tracking_system)){
    print(paste("START: Auto-orient + capsule for",tracking_system, sep =" "))
    #################################################################################################################################################################################################################
    ###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
    ###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
    ###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
    #################################################################################################################################################################################################################
    data_list         <- list (paste(dir_data,manually_oriented_files[which(grepl(tracking_system,manually_oriented_files))],sep="/")) ###here list all the myrmidon files containing oriented data
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
    ####Linda: adapting script from this point forward #### before this point script from "2_AUTO_ORIENT", beyond this point "CAPSULE_REASSIGNMENT_Luke"
      
    ###################################################################################
    ############### list the experiment files to ORIENT ###############################
    ToOrient_data_list         <- paste(dir_data,ToOrient_file[which(grepl(tracking_system,ToOrient_file))],sep="/")
    
    to_keep_2 <- c(ls(),"to_keep_2","ToOrient_myr_file")
    for (ToOrient_myr_file in ToOrient_data_list){
      
      ToOrient_exp_name_full <- unlist(strsplit(ToOrient_myr_file,split="/"))[length(unlist(strsplit(ToOrient_myr_file,split="/")))]
      ToOrient_exp_name <- unlist(strsplit(ToOrient_exp_name_full,split="[.]"))[1]
      
      print(paste("Assign Capsules to",ToOrient_exp_name, sep =" "))
      
      ###open the myrmidon file to orient nd the metadata myrmidon file
      ToOrient_data <- fmExperimentOpen(ToOrient_myr_file)
      ToOrient_ants <- ToOrient_data$ants
      
      
      #### DELETE ANY CAPSULE ALREADY PRESENT
      
      # delete individuals' capsule data IF PRESENT
      for (ant in ToOrient_data$ants){
        ToOrient_data$ants[[ant$ID]]$clearCapsules()
      }
      
      # delete the capsule shapes IF PRESENT
      if (length(ToOrient_data$antShapeTypeNames)>0) {
        for (caps in 1:length(ToOrient_data$antShapeTypeNames)){
          ToOrient_data$deleteAntShapeType(caps)
        }
      }
      
      
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
      
      # ###POSITIONS for all, cut trajs after
      # FROM <- fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$start) ###experiment start time
      # TO   <- fmTimeCreate(offset=fmQueryGetDataInformations(ToOrient_data)$start + 12*3600  ) ###experiment start time
      # max_gap <- fmHour(24*365)  ###this parameter is very important to - use a super large value to make sure you get only one trajectory per ant!!!
      # positions <- fmQueryComputeAntTrajectories(ToOrient_data,start = FROM,end = TO,maximumGap = max_gap,computeZones = TRUE)
      # 
      
      ###copy individual metadata, orient ant and add capsules
      to_keep_3 <- c(ls(),"to_keep_3","ant_INDEX")
      for (ant_INDEX in 1:length(ToOrient_data$ants)){   ####LOOP OVER ANTS
        #   ###get corresponding ant from metadata
        
        ###orient ant and define capsule: loop over identifications
        ###SKIP FOR QUEEN SO THAT YOU CAN ORIENT THE QUEEN AND CREATE CAPSULES MANUALLY IN FORTSTUDIO
        if (ToOrient_data$ants[[ant_INDEX]]$ID!=1){
          to_keep_4 <- c(ls(),"to_keep_4","identif")
          # for (identif in 1:length(ToOrient_data$ants[[ant_INDEX]]$identifications)){
          # }## identif
          ###finally, once ant has been oriented for all identifications, add its capsules
          ##assign capule numbers that match the order of the looped capsule names positions
          capsule_number <- 0
          for (capsule_name in unlist(ToOrient_data$antShapeTypeNames)) {
            capsule_number <- capsule_number +1
            #MAKE SURE THERE IS CAPSULE MATCHING, TO AVOID MIXING UP SHAPE INDEXES
            capsule_ratios <- capsule_list[[capsule_name]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
            capsule_coords <- ant_length_px*capsule_ratios
            
            ToOrient_data$ants[[ant_INDEX]]$addCapsule(capsule_number, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
            
          }  #ADD CAPSULES #IF THIS RESULTS IN A ERROR, CHECK PREVIOUS VERSION ON GIT OR COMPARE WITH DataPrep4_Clone-capule-queens-only_v082.R
          
          #} ## if no trajectory (eg. dead from)
        }###if not queen
        #rm(list=ls()[which(!ls()%in%to_keep_3)])
        #gc()
      }# LOOP ANTS
      
      #
      # SHORTHEN FILE NAMES!!!!
      #
      
      
      ToOrient_treatment <- unlist(strsplit(ToOrient_exp_name,split="_"))[1]
      ToOrient_colony <- unlist(strsplit(ToOrient_exp_name,split="_"))[3]
      output_file_name <- paste0(dir_data, "/", tracking_system, "_", ToOrient_colony, "/", ToOrient_treatment, "_", tracking_system, "_", ToOrient_colony,  "_automatically_oriented_CapsuleDef3_metadata_deaths_q_zones.myrmidon")
      ToOrient_data$save(output_file_name)
      print(paste("Saving ", tracking_system, " ", ToOrient_colony, " as ", output_file_name,sep =""))
      
      rm(list=(c("ToOrient_data"))) #remove experiment
      # }else{print(paste0(paste0(sub("\\..*", "", ToOrient_myr_file),"_AutoOriented.myrmidon")," already exists! Skip!"))} # ToOrient FILES LOOP #closes: for (ToOrient_myr_file in ToOrient_data_list){
      
      #rm (list = ls() [which(!ls()%in%to_keep_2)] )
      #gc()
      #} # IF EXISTS, SKIP PROCESSING
      
      #rm(list=ls()[which(!ls()%in%to_keep_1)])
      #gc()
    } #To orient LOOP
  } # CAPSULE DEF
}
  
  cat("LOOP ENDED!! \n Go to fort-studio to check things look all right
\n and enjoy the capsules!)
")
  
  