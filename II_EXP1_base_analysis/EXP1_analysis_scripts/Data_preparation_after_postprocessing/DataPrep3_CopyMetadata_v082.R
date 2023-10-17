rm(list=ls())
gc()
mallinfo::malloc.trim(0L)


########################################################################
############# ASSIGN METADATA FROM METADATA-RICH FILES #################
# this script should be executed after the orientation of the files
# AntsCreated_AutoOriented files can inherit Metadata info  from unoriented files ("AntsCreated_DeathRecord_NoOrient.myrmidon") or from manually oriented files ("ManOriented.myrmidon")
# This loop takes approx 1 minute per file

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

######load necessary libraries
library(FortMyrmidon) ####R bindings

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


### directory of data and myrmidon files
# dir_data <- '/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA'
dir_data <- "/media/bzniks/DISK4/ADRIANO/EXPERIMENT_DATA"


#### ACCESS FILES
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(dir_data, n = 1)

#select REP folders
files_list <- files_list[grep("REP",files_list)]

#### OPEN REPLICATE
EXP_list <- NULL
# replicate folder
for (REP.n in 1:length(files_list)) {
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = ".myrmidon")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  for (variable in REP.files) {
    
    #get substring in variable until R9BS_
    REP_treat_name <- sub("\\_.*", "", variable)
    treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
    
    EXP_list <- rbind(EXP_list ,data.frame(REP_treat_name,path_name=paste(REP.folder,variable,sep="/")))
    
  }
}

#Select auto-oriented files
AntsCreated_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "AutoOriented")),]
AntsCreated_list <- AntsCreated_list[which(!grepl(AntsCreated_list$path_name,pattern = "_withMetaData")),]


# Select the metadata-rich files
Metadata_list <- EXP_list[which(grepl(EXP_list$path_name,pattern = "DeathRecord_NoOrient|ManOriented")),]
Metadata_list <- Metadata_list[which(!grepl(Metadata_list$path_name,pattern = "base")),]

#check that each metadata REP file has its associated AutoOriented file
Metadata_list$REP_treat_name %in% AntsCreated_list$REP_treat_name

### Loop through all the directories in the dir_folder
for (REP in Metadata_list$REP_treat_name){
  print(paste("Assign Metadata to",REP, sep =" "))
  
# temporary file/vars assignment
#AntsCreated_myr_file <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated.myrmidon"
#Metadata_myr_file <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated_DeathRecord_NoOrient.myrmidon"
# AntsCreated_myr_file <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated.myrmidon"
# Metadata_myr_file <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated_DeathRecord_NoOrient.myrmidon"

AntsCreated_myr_file <- AntsCreated_list[which(AntsCreated_list$REP_treat_name==REP),]$path_name
Metadata_myr_file <- Metadata_list[which(Metadata_list$REP_treat_name==REP),]$path_name

# if the _AutoOriented file file doesn't exist, then continue
if ( !file.exists(paste0(sub("\\..*", "", AntsCreated_myr_file),"_withMetaData.myrmidon"))) {
  

Metadata_exp <- fmExperimentOpen(Metadata_myr_file) 
Metadata_ants <- Metadata_exp$ants

### open antsCreated base file
AntsCreated  <- fmExperimentOpen(AntsCreated_myr_file)
creat_exp_name <- unlist(strsplit(AntsCreated_myr_file,split="/"))[length(unlist(strsplit(AntsCreated_myr_file,split="/")))]
AntsCreated_ants <- AntsCreated$ants
#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myr_file),"_withMetaData.myrmidon")) # file now exists



### CHECK VISUALLY IF RETAG HAPPENED OR FROM METADATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## THEN EXPLORE THE STRUCURE, MAY BE NEEDED LOOP OVER LIST ELEMENTS

# # check if retags happened:
# for (ant in Metadata_ants){
#  # individual  <- ant$ID
#   print( ant$ID)
#   print(ant$identifications)
# }


# ############# ASSIGN IDENTIFICATION ########################
for (ant in Metadata_ants){
  individual  <- ant$ID
  #print(ant$identifications)
  #if AntID matches
  if (AntsCreated_ants[[individual]]$identifications[[1]]$targetAntID==Metadata_ants[[individual]]$identifications[[1]]$targetAntID){
    #remove identification
    #add identification
    AntsCreated$addIdentification(antID= Metadata_ants[[individual]]$identifications[[1]]$targetAntID,
                                tagID= Metadata_ants[[individual]]$identifications[[1]]$tagValue,
                                start= Metadata_ants[[individual]]$identifications[[1]]$start,
                                end=   Metadata_ants[[individual]]$identifications[[1]]$end)
  }
}


# #CHECK copy
# for (ant in AntsCreated_ants){
#   print(ant$identifications)
#   ant$identifications
# }

############# ASSIGN METADATA KEYS AND VALUES #########################
list_keys <- list()
#assign metadata keys
for (KEY in   1:length(Metadata_exp$metaDataKeys)) {
  key <- names(Metadata_exp$metaDataKeys[KEY])
  defaultvalue <- unname(Metadata_exp$metaDataKeys[KEY][[1]])
  AntsCreated$setMetaDataKey(key,defaultvalue)
  #check
  #AntsCreated$metaDataKeys[KEY]
  list_keys <- c(list_keys,key)
}


for (ant in Metadata_ants){
  individual  <- ant$ID
  #extract metadata info per key
  for (METADATA_KEY in list_keys) {
    for (ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
      # AntsCreated$ants[[individual]]$deleteValue(key=METADATA_KEY,time=NA) #not possible to delete the first row of metadata as NA is not a valid time
      # if ("metadata value is equal in both files, skip assign") { #not possible as above
      if ( is.na(ant$getValues(METADATA_KEY)[ROW,"times"])){ # if there is an NA, set time as epoch start (fmTimeSinceEver doesn't work as it assigns currrent time)
        AntsCreated$ants[[individual]]$setValue(key=METADATA_KEY,value=ant$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeCreate(offset = 0)) #fmTimeSinceEver()
      }else{
        AntsCreated$ants[[individual]]$setValue(METADATA_KEY,value=ant$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeCreate(offset=ant$getValues(METADATA_KEY)[ROW,"times"]))
      }#assign value
     # }#skip if same
    }#ROW
  }#METADATA_KEY
}#ant

#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myr_file),"_withMetaData.myrmidon")) # file now exists

#check
# for (ant in AntsCreated_ants){
#   individual  <- ant$ID
#   #extract metadata info per key
#   for (METADATA_KEY in list_keys) {
#     cat(paste("ANT",individual,"\n","key:",METADATA_KEY,sep=" "))
#     
#     print(ant$getValues(METADATA_KEY))
#   }}

# cleaning: for value = base.value and time is 1970-01-01 01:00:00, delete event 
for (ant in AntsCreated_ants){
  individual  <- ant$ID
  #extract metadata info per key
  for (METADATA_KEY in list_keys) {
    for (ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
      #if there is an assigned time
      if (!is.na(ant$getValues(METADATA_KEY)[ROW,"times"])){
      #if value is equal to base.value 
        if (ant$getValues(METADATA_KEY)[ROW,"values"]==AntsCreated$metaDataKeys[[METADATA_KEY]]) {
          AntsCreated$ants[[individual]]$deleteValue(METADATA_KEY,time=fmTimeCreate(offset = 0))
        }}}}}


############# ASSIGN ZONES #########################
#assign zones
list_zones <- NULL
for (ZONE in   1:length(Metadata_exp$spaces[[1]]$zones)) {
  zone <- Metadata_exp$spaces[[1]]$zones[[ZONE]]$name
  #defaultvalue <- unname(Metadata_exp$metaDataKeys[ZONE][[1]])
  AntsCreated$spaces[[1]]$createZone(zone)
  #check
  #AntsCreated$spaces[[1]]$zones
  list_zones <- c(list_zones,zone)
}

#Assign SHAPE to ZONE (geometry to the nest zone)
#extract metadata info per key
for (ZONE.1 in   1:length(AntsCreated$spaces[[1]]$zones)) {
#for (ZONE_KEY in list_zones) {
  zone_definition <- Metadata_exp$spaces[[1]]$zones[[ZONE.1]]$definitions
  #assign shapes
  AntsCreated$spaces[[1]]$zones[[ZONE.1]]$addDefinition(shapes=zone_definition[[1]][["shapes"]],start= fmTimeSinceEver(),end=fmTimeForever())
  }#ZONE.1

#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myr_file),"_withMetaData.myrmidon"))
print(paste("Saving",creat_exp_name, "as _withMetaData.myrmidon",sep =" "))
}else{print(paste0(paste0(sub("\\..*", "", AntsCreated_myr_file),"_withMetaData.myrmidon")," already exists! Skip!"))} # withMetaData FILES LOOP
  rm(list=(c("Metadata_exp","AntsCreated"))) #remove experiments
} # IF EXISTS, SKIP PROCESSING

warning("ASSIGN IDENTIFICATION not working! Currently turned off")

print("NOTE: \nAs it is not possible to assign a metadata value change without time info (NA), value for \"Is queen\" results as a timed change starting on epoch start (1st Jan 1970)")

cat("LOOP ENDED!! \n Go to fort-studio to check things look all right
\n AND OVERWRITE INFO FOR QUEEN!!! (tag size, manual orientation, manual capsules)
\n \n Do check identifications to make sure there is a exact correspondance between ants in both Ants Created files and metadata provided files
  ")
