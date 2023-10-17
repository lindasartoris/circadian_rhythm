###compatible with Myrmidon 0.8.3

###edit data_path #done
###create a copy of the myrmidon file that contains the zones and rename it "[...]_zones_Linda_2023" #done
###to avoid losing data don't overwrite but create new myrmidon files with ending "[...]_zones"


rm(list=ls())

library(FortMyrmidon)

data_path <-  '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/tracking'

file_with_zones    <- paste(data_path,"alleline_c05/DF_alleline_c05_automatically_oriented_CapsuleDef2018_metadata_deaths_q_zones_Linda_2023.myrmidon",sep="/")

files_without_zone <- list.files(path=data_path,pattern="automatically_oriented_CapsuleDef2018_metadata_deaths_q.myrmidon",full.names=T, recursive = T)
files_without_zone <- files_without_zone[which(!grepl(file_with_zones,files_without_zone))]
files_without_zone <- files_without_zone[which(!grepl("zones_Linda_2023",files_without_zone))] 

###extract zone information
zoned <- fmExperimentOpen(file_with_zones) 
zone_list <- list()

###get list of zones
for (space in zoned$spaces){
  zone_list <- c(zone_list,space$zones)
}

###get information for each zone
for (zone in zone_list){
  ID   <- zone$ID
  name <- zone$name
  start <- zone$definitions[[1]]$start
  end   <- zone$definitions[[1]]$end
  shapes <- zone$definitions[[1]]$shapes
}

###now create zones for each unzoned file 
for (unzoned_file in files_without_zone){
  if (!file.exists(gsub("\\.myrmidon","_zones.myrmidon",unzoned_file))){
    print(unzoned_file)
    print("Reading tracking data...")
    unzoned <- fmExperimentOpen(unzoned_file) 
    
    ###check that there are no zones, or delete existing zones 
    for (space in unzoned$spaces){
      if (length(space$zones)>0){
        print("Deleting existing zones...")
        for (zone in space$zones){
          space$deleteZone(zone$ID)
        }
      }
    }
    
    print("Duplicating desired zones...")
    for (space in unzoned$spaces){ ## do this for each space
      for (zone in zone_list){
        ###create zone
        space$createZone(zone$name)
        ###load new zone into object
        new_zone <- space$zones[[length(space$zones)]]
        ###add a new zone definition
        new_zone$addDefinition(zone$definitions[[1]]$shapes ,   zone$definitions[[1]]$start , zone$definitions[[1]]$end  )
      }
    }
    
    print("Saving file...")
    unzoned$save(gsub("\\.myrmidon","_zones.myrmidon",unzoned_file))
    
    rm(list=c("unzoned"))
    
  }else{
    print(paste(unzoned_file),": zoned file already exists; skipping.")
  }
}
