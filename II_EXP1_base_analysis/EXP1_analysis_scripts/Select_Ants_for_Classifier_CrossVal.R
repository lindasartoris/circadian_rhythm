
###########################################################################
## RANDOMLY SELECT NON-FOCAL INDIVIDUALS FOR CLASSIFIER CROSS-VALIDATION ##
###########################################################################

## The selected ants should be followed for: 
# 30 mins after exposure (between 2 and 32 mins post-isolation in 15 mins blocks as done with Vasudha)
# only for grooming

DATADIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data/"

Metadata_Exp1 <- read.table(paste(DATADIR,"Metadata_Exp1_2021_2022-10-12.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#select only non-exposed individuals
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$Exposed==FALSE),]
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$IsAlive==TRUE),]
#exclude tag falls
Metadata_Exp1 <- Metadata_Exp1[!grepl("TagFall", Metadata_Exp1$Comment),]

# RANDOMLY 1 ant per colony
# randomly choose only one row in each Replicate
Metadata_Exp1$Chosen <- 0
Metadata_Exp1[-tapply(-seq_along(Metadata_Exp1$REP_treat),Metadata_Exp1$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_Chosen <- Metadata_Exp1[which(Metadata_Exp1$Chosen==1),]

# exclude R3SP and R9SP
Metadata_Chosen <- Metadata_Chosen[!Metadata_Chosen$REP_treat %in% c("R3SP","R9SP"),]
# -  1 col x treatment - SP, BP, SS, BS 
# Metadata_Chosen[-tapply(-seq_along(Metadata_Chosen$treatment),Metadata_Chosen$treatment, sample, size=1),]$Chosen <- 2
# Metadata_Chosen <- Metadata_Chosen[which(Metadata_Chosen$Chosen==2),]

#check if the data for CrossValidation has already been collected
## TO BE RERUN ON THE UPDATED FILE Classifier_CrossVal_2022-08-26.txt WHEN THE FILES AT THE BACK ARE FIXED!
Pre_selected_REPS <- read.table(paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_2022-08-17.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
Pre_selected_REPS <- Pre_selected_REPS$REP_treat
Metadata_Chosen <- Metadata_Chosen[!Metadata_Chosen$REP_treat %in% Pre_selected_REPS, ]

#Save output
write.table(Metadata_Chosen,file=paste(DATADIR,"/Select_NON_FOCAL_Ants_for_Classifier_CrossVal_",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")




###########################################################################
### RANDOMLY SELECT FOCAL INDIVIDUALS FOR CLASSIFIER CROSS-VALIDATION #####
###########################################################################

## The selected ants should be followed for: 
# 30 mins after exposure (between 2 and 32 mins post-isolation in 15 mins blocks as done with Vasudha)
# only for grooming

DATADIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

Metadata_Exp1 <- read.table(paste(DATADIR,"/Metadata_Exp1_2021.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#select only exposed nurses
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$Exposed==TRUE),]
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$IsAlive==TRUE),]
#exclude tag falls
Metadata_Exp1 <- Metadata_Exp1[!grepl("TagFall", Metadata_Exp1$Comment),]

# RANDOMLY 1 ant per colony (only exposed)
# randomly choose only one row in each Replicate
Metadata_Exp1$Chosen <- 0
Metadata_Exp1[-tapply(-seq_along(Metadata_Exp1$REP_treat),Metadata_Exp1$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_Chosen <- Metadata_Exp1[which(Metadata_Exp1$Chosen==1),]

# exclude R3SP and R9SP
Metadata_Chosen <- Metadata_Chosen[!Metadata_Chosen$REP_treat %in% c("R3SP","R9SP"),]
# -  1 col x treatment - SP, BP, SS, BS 
# Metadata_Chosen[-tapply(-seq_along(Metadata_Chosen$treatment),Metadata_Chosen$treatment, sample, size=1),]$Chosen <- 2
# Metadata_Chosen <- Metadata_Chosen[which(Metadata_Chosen$Chosen==2),]

#check if the data for CrossValidation has already been collected
## TO BE RERUN ON THE UPDATED FIILE Classifier_CrossVal_2022-08-26.txt WHEN THE FILES AT THE BACK ARE FIXED!
Pre_selected_REPS <- read.table(paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_2022-08-17.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
Pre_selected_REPS <- Pre_selected_REPS$REP_treat
Metadata_Chosen <- Metadata_Chosen[!Metadata_Chosen$REP_treat %in% Pre_selected_REPS, ]

#Save output
write.table(Metadata_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")


