rm(list=ls())

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(circular)
library(R.utils)

###source C++ movement direction program
sourceCpp("/home/cf19810/Documents/scriptsR/Automated_Ant_Orientation_Nathalie/Get_Movement_Angle.cpp")

#################################################################################################################################################################################################################
###########STEP 1 - use manually oriented data to extract important information about AntPose and Capsules
###########          IN YOUR PARTICULAR SPECIES AND EXPERIMENTAL SETTINGS  
###########IMPORTANT: FOR THIS TO WORK YOU NEED TO HAVE DEFINED THE SAME CAPSULES WITH THE SAME NAMES ACROSS ALL YOUR MANUALLY ORIENTED DATA FILES
#################################################################################################################################################################################################################
data_list         <- list ("/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4.myrmidon") ###here list all the myrmidon files containing oriented data
#ALL COLONIES OF THE SAME TRACKNG BOX ORIENTED HAVE TO BE LOADED AT THIS STAGE 

oriented_metadata <- NULL
capsule_list <- list()
for (myrmidon_file in data_list){
  experiment_name <- unlist(strsplit(myrmidon_file,split="/"))[length(unlist(strsplit(myrmidon_file,split="/")))]
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
    
    
    ###extract offset betewen tag centre and ant centre
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
#mean_y_ant_coord <- mean(oriented_metadata$y_ant_coord) ###this is expected to be 0 or near 0 (check!) because the tag should be as likely to be to the right or to the left of the ant's bilateral symmetry line
mean_y_ant_coord <- 0 ##set it to zero manually because we don't expect there to be a consistent bias in deviation
###Furthermore, get the average worker length from the data
mean_worker_length_px <- mean(oriented_metadata$length_px)
###Finally, get information on each capsule
for (caps in 1:length(capsule_list)){
  capsule_list[[caps]] <- colMeans(capsule_list[[caps]][,which(grepl("ratio",names(capsule_list[[caps]])))])
}

###create blank myrmidon file
tracking_data <- fmExperimentCreate('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon') # no file created yet
tracking_data$save('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon') # file now exists

###create space
s <- tracking_data$createSpace('prideaux')
printf("Space '%s' has ID: %d\n",s$name,s$ID)
# outputs: Space 'nest' has ID: 1
tracking_data$save('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon') # file now exists


###add tracking data directory
tddURI_1  <- tracking_data$addTrackingDataDirectory(s$ID,"/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21.0000")
tddURI_2  <- tracking_data$addTrackingDataDirectory(s$ID,"/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21.0001")
tddURI_3  <- tracking_data$addTrackingDataDirectory(s$ID,"/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21.0002")

tracking_data$save('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon') # file now exists

###create ants
tag_statistics <- fmQueryComputeTagStatistics(tracking_data)
for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}
tracking_data$save('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon')

##Check: print identifications
ants <- tracking_data$ants
for (a in ants) {
  printf("Ant %s is identified by:\n", fmFormatAntID(a$ID))
  for (i in a$identifications){
    printf(" * %s\n", capture.output(i))
  }
}

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
###option 1: using all data - using start and end time from the experiment metadata
#from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) ###experiment start time
#to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$end  ) ###experiment end time

###option 2: using all data - using since Ever and Forwever functions
# from <- fmTimeSinceEver()
# to   <- fmTimeForever()

###option 3: if dataset is too large - too time consuming, use only 24 hours (for example)
from <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start) ###experiment start time
to   <- fmTimeCreate(offset=fmQueryGetDataInformations(tracking_data)$start + 24*3600  ) ###experiment start time

###option 4:..... whateever makes most sense 

###Then compute trajectories
max_gap <- fmHour(24*365)  ###this parameter is very important to - use a super large value to make sure you get only one trajectory per ant!!!
positions <- fmQueryComputeAntTrajectories(tracking_data,start = from,end = to,maximumGap = max_gap,computeZones = TRUE)

#hard-wire ant correspondence between trajectories_summary and trajectorues
positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list

##define a max temporal gap for which you are happy to calculate a movement angle; e.g. 0.5 s
max_time_gap <- 0.2 #0.5

##define a minimum distance moved, as you don't want to use noise or small shifts in position in this calculation; e.g. 30 pix (to think about)
min_dist_moved <- 15 #30

for (i in 1:length(ants)){
  ####to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
  traj <- positions$trajectories [[   positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"antID_str"]    ]]
  
  ###feed traj to c++ program
  traj <- cbind(traj,add_angles(traj,max_time_gap,min_dist_moved))
  
  ## get mean deviation angle between body and tag - the ant angle is equal to minus the Tag minus Movement angle output by C++ program
  AntAngle <- as.numeric(- mean(circular(na.omit(traj$Tag_minus_Movement_Angle),units="radians",zero=0)))
  ##now use trigonometry to calculate the pose, using AntAngle
  x_tag_coord <- mean_x_ant_coord*cos(AntAngle) - mean_y_ant_coord*sin(AntAngle)
  y_tag_coord <- mean_x_ant_coord*sin(AntAngle) + mean_y_ant_coord*cos(AntAngle)
  
  ##write this into ant metadata
  for (id in ants[[i]]$identifications){
    id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
  }
  
  ##also add this to trajectories_summary
  positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"ant_angle"] <- AntAngle
  positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"x_tag_coord"] <- x_tag_coord
  positions$trajectories_summary[which(positions$trajectories_summary$antID==ants[[i]]$ID),"y_tag_coord"] <- y_tag_coord
  
  ###finally, for each ant, add capsules using mean_ant_length and capsule_list
  for (caps in 1:length(capsule_list)){
    capsule_ratios <- capsule_list[[caps]]; names(capsule_ratios) <- gsub("_ratio","",names(capsule_ratios))
    capsule_coords <- mean_worker_length_px*capsule_ratios
    
    
    ants[[i]]$addCapsule(caps, fmCapsuleCreate(c1 = c(capsule_coords["c1_x"],capsule_coords["c1_y"]), c2 = c(capsule_coords["c2_x"],capsule_coords["c2_y"]), r1 = capsule_coords["r1"], r2 = capsule_coords["r2"] ) )
    
  }
  
}
tracking_data$save('/home/cf19810/Documents/Ants_behaviour_analysis/Data/R3SP_13-03-21_Capsule_Zones_defined_BODY-HEAD_20April22_4_auto_oriented_meanlength.myrmidon')
##LASTLY!! Go to fort-studio to check things look all right
### AND TO OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)


#HOW TO:
#assign ant size

#MAYBE IN A SIMILAR FASHION TO

##write this into ant metadata
# for (id in ants[[i]]$identifications){
#   id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)
# }

# #THE OBJECT TO ACCESS SHOULD BE
# oriented_data$measurementTypeNames
# #IT IS POSSIBLE TO EXTRACT THE INFORMATION, HOW TO REINTRODUCE IT?
# fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px
# 
# #python syntax FOR THE HEAD_TAIL_MEASUREMENT_TYPE_ID
# measurementTypeID=exp.HEAD_TAIL_MEASUREMENT_TYPE_ID
# 
# 
# 
# #possibly involeved functions
# measurementTypeID
# SetMeasurementTypeName
# setUserDefinedAntPose
# #id$setUserDefinedAntPose(c(x_tag_coord,y_tag_coord), AntAngle)

## check
# https://formicidae-tracker.github.io/myrmidon/latest/api/cpp/experiment.html?highlight=createmeasurementtype#_CPPv4N4fort8myrmidon10Experiment21CreateMeasurementTypeERKNSt6stringE
# https://formicidae-tracker.github.io/myrmidon/latest/api/cpp/ant_identification.html#classfort_1_1myrmidon_1_1Identification_1a532b4676510381aef6e174f25a595211
# https://formicidae-tracker.github.io/myrmidon/latest/api/cpp/typedef.html#_CPPv4N4fort8myrmidon26HEAD_TAIL_MEASUREMENT_TYPEE
# https://github.com/formicidae-tracker/myrmidon/blob/master/bindings/R/FortMyrmidon/src/ant.cpp


