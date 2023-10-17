#be aware that myrmidon file names might have changed
#treatment (DF, DS, NF or NS) was added to the myrmidon file name [treatment]_[tracking system]_[colony]


#### prerequisites ####
rm(list=ls())

# load libraries
library(FortMyrmidon) # R bindings
library(R.utils)      # printf()
library(Rcpp)         # contains sourceCpp (used for ant orientation)
library(circular)     # used for ant orientation
library(data.table)   # used to save files fwrite(list(myVector), file = "myFile.csv")
library(stringr)
library(reader)



dir_data <- '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/tracking'
dir_data_foragerID <- '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/forager_IDs'



#target_file <- metadata[which(metadata$tracking_system==tracking_system),"file_name"]

#output_file_name <- paste(dir_data,"/",target_file,"/",target_file,"_automatically_oriented_",suffix,".myrmidon",sep="")

#foragerID_file_name <- paste(dir_data_foragerID,"/",target_colony,"/",target_colony,".myrmidon",sep="")

###functions
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
  metadata <- data.frame(file_name=list.dirs(path = dir_data,recursive=F,full.names = F),tracking_system=unlist(lapply (   list.dirs(path = dir_data,recursive=F,full.names = F)    ,  function(x)unlist(strsplit (x,split="_"))[1])), colony=unlist(lapply (   list.dirs(path = dir_data,recursive=F,full.names = F)    ,  function(x)unlist(strsplit (x,split="_"))[2])))
  return(metadata)
  
}

create_metadata_foragerID <- function (dir_data){
  metadata_foragerID <- data.frame(file_name=list.dirs(path = dir_data,recursive=F,full.names = F), colony=unlist(lapply (   list.dirs(path = dir_data,recursive=F,full.names = F)    ,  function(x)unlist(strsplit (x,split="_"))[3])))
  return(metadata_foragerID)
  
}


##### 5. Add meta data
##### 5.1 create myrmidon files for forager_ID tracking
##### create myrmidon files for the exposed foragers (forager_ID tracking files)
dir_data_foragerID <- '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/forager_IDs'

metadata_foragerID <- create_metadata_foragerID(dir_data_foragerID)

colony <- metadata_foragerID$colony
tag_size <- 0.76
for (target_colony  in    colony){
  
  ####define the output myridon file name
  output_file_name <- paste(dir_data_foragerID,"/forager_ID_",target_colony,"/foragerID_",target_colony,".myrmidon",sep="")
  
  ###create blank myrmidon file
  tracking_data <- fmExperimentCreate(output_file_name) # no file created yet
  tracking_data$save(output_file_name) # file now exists
  
  ###create space
  s <- tracking_data$createSpace("esterhase")
  printf("Space '%s' has ID: %d\n",s$name,s$ID)
  # outputs: Space 'nest' has ID: 1
  tracking_data$save(output_file_name) # file now exists
  
  
  ###add tracking data directory
  tracking_folder <- list.dirs.depth.n(paste(dir_data_foragerID,"/forager_ID_",target_colony,sep=""), n = 1)
 
  tddURI <- tracking_data$addTrackingDataDirectory(s$ID,tracking_folder,TRUE)
  tracking_data$save(output_file_name) # file now exists
  
  
  tracking_data$defaultTagSize <- tag_size
  ###create ants
  tag_statistics <- fmQueryComputeTagStatistics(tracking_data)
  for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
      a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
      identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
      print(identification)
  }
  tracking_data$save(output_file_name) 
}



#### 5.2.1 Create metadata keys
# add to auto_orient files!
### directory of data and list of myrmidon files we need 
dir_data <- '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/tracking'
dir_data_foragerID <- '/media/lsartori/Seagate Portable Drive/circadian_rhythm_2022_experiment/forager_IDs'

automatically_oriented_files_CapsuleDef2018 <- list.files(path=dir_data, pattern="automatically_oriented_CapsuleDef2018", recursive=T)

for ( automatically_oriented_files_text in c("automatically_oriented_files_CapsuleDef2018")){
  automatically_oriented_files <- get(automatically_oriented_files_text)
  
  for ( myrmidon_file in automatically_oriented_files){
    colony <- unlist(strsplit(myrmidon_file,split="_"))[3] #this would have to be adapted for new file names in case it is run again
    print(colony)
    target_file <- unlist(strsplit(myrmidon_file,split="/"))[1]
    print(target_file)
    output_file_name <- paste(dir_data,"/",target_file,"/",target_file,"_automatically_oriented_CapsuleDef2018_metadata.myrmidon",sep="")
    fort_data <- fmExperimentOpen(paste(dir_data,"/",myrmidon_file,sep=""))
    fort_data$setMetaDataKey(key = "IsQueen",     default_Value = FALSE)
    fort_data$setMetaDataKey(key = "Exposed",   default_Value = FALSE)
    fort_data$setMetaDataKey(key = "IsAlive",     default_Value = TRUE)
    fort_data$setMetaDataKey(key = "treatment",   default_Value = "NA") # treated ants will get sham or fungus if necessary as key
    fort_data$setMetaDataKey(key = "tag_reoriented", default_Value = FALSE)
    #fort_data$spaces[[1]]$createZone(name = "nest") # create zones to be defined manually in the fort files 
    #fort_data$spaces[[1]]$createZone(name = "arena") # not necessary in my case (do once manually and then clone)
    #fort_data$spaces[[1]]$createZone(name = "water_right") 
    #fort_data$spaces[[1]]$createZone(name = "sugar_left")
    fort_data$save(output_file_name)
    
    #### 5.2.2 Add treated forager ID myrmidon files  
    if (colony != "c38") { # because this colony does not have a myrmidon file (will be added manually)
    treatment_data <- fmExperimentOpen(paste(dir_data_foragerID,"/forager_ID_",colony,"/foragerID_", colony,".myrmidon",sep=""))     # create vector of the treated ants (separately recorded, separate myrmidon file)
    treated_ants <- treatment_data$ants
    tag_value_vector <- NULL
    tag_values <- NULL
    
      for (z in treated_ants) { #go through every ant of the separately recorded individuals
        tag_values <- z$identifications[[1]]$tagValue #extract tag value (= tag ID)
        tag_value_vector <- rbind(tag_value_vector, data.frame(tag_values))
      }
      ants <- fort_data$ants  # for each ant adjust the meta data if it is the queen or a treated worker
      
      for (x in ants) {
        if (x$identifications[[1]]$tagValue==0) { # Queen always has the 000 tag
          x$setValue("IsQueen", TRUE, time = fmTimeSinceEver())}
        if(is.element(x$identifications[[1]]$tagValue, as.matrix(tag_value_vector))) {
          x$setValue(key="Exposed", value = TRUE, time = fmTimeSinceEver())}
      }
      fort_data$save(output_file_name)
    } else {
      print(paste("no myrmidon file for", colony, sep=" "))
    }
  }
}


#find all the replicates that have a corrupted first streaming file
corrupted_first_stream <- list.files(path=dir_data, pattern="stream.0000.dis.mp4", recursive=T)
length(corrupted_first_stream)
#36 files
