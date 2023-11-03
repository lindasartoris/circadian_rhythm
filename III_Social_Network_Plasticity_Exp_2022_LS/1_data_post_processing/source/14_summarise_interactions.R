####14_summarise_interactions.R#####

####Takes an interaction list as an input, and calculates:
###                     - the total duration of interactions of each ant to each of four categories of individuals (queen, nurses, untreated foragers and treated foragers)
###                     - the overall distribution of interactions according to task groups and ages

###Created by Nathalie Stroeymeyt
###Modified by Adriano Wanderlingh to work with FORT formicidae Tracking data. Mods tagged with the comment "AW". script wide mods cited here below.
###Modified by Nathalie Stroeymeyt to include number of events in addition to duration
###with adaptations by Linda Sartoris

#Script wide mods AW
# - replaced before/after with pre/post

# when ran for the grooming interactions it should only be run for the "observed" folders
####################################
to_keep_ori <- to_keep

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#### get input file list
input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
setwd(input_path)  
input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
input_folders        <- input_folders[which(input_folders!="")]

outputfolder1 <- paste(data_path,"/processed_data/individual_behaviour/random_vs_observed",sep="")
if (!file.exists(outputfolder1)){dir.create(outputfolder1,recursive = T)}

summary_dol <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_file","network_files","summary_interactions","summary_interactions_grooming","summary_pairs","all_interactions")
for (input_folder in input_folders){ ### LS: change this back so it iterates through all input_folders
  print(input_folder)
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  if (input_folder=="observed"&grepl("main",data_path)){
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
  }
  summary_interactions <- NULL
  summary_interactions_grooming <- NULL
  summary_pairs        <- NULL
  all_interactions     <- NULL
  for (network_file in network_files){
    # network_file <- network_files[1] # temp
    cat("\r",network_file)
    ####get file metadata
    root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))])  # LS: replace grepl("colony", ...) with grepl("interactions")
    components         <- unlist(strsplit(root_name,split="_"))
    # colony             <- components[grepl("colony",components)]
    colony             <- unlist(strsplit(root_name,split="_"))[1] # LS
    treatment          <- info[which(info$colony==colony),"treatment"] #AW: no need for as.numeric() 
    treatment_circadian<- unlist(strsplit(unlist(strsplit(root_name,split="_"))[2],split="\\."))[2] # LS
    colony_size        <- info[which(info$colony==colony),"colony_size"]
    
    if (!all(!grepl("PreTreatment",components))){period <- "pre"}else{period <- "post"} 
    time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
    time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
    
    # LS: add period_detail so it can be added in extra column later
    if (period=="pre"){
      if (components[4] %in% c("TH-24", "TH-21", "TH-18")){
        period_detail <- "pre1"
      }else{
        period_detail <- "pre2"
      }
    }else{
      period_detail <- "post"
    }
    
    # LS: add period_circadian
    if (period_detail=="pre1"|period_detail=="post"){
      period_circadian <- treatment_circadian
    }else{
      period_circadian <- ifelse(treatment_circadian=="day","night","day")
    }
    
    if (grepl("age",data_path)){
      colony_ages <- ages[which(ages$colony==colony),]
    }
    
    ####get appropriate task_group list, treated list and tag
    colony_treated     <- treated[which(treated$colony==colony),"tag"] #AW
    colony_task_group  <- task_groups[which(task_groups$colony==colony),]
    queenid            <- as.character(colony_task_group[which(colony_task_group$task_group=="queen"),"tag"]) #AW: call specific queen tag instead of fixed 665. it has to be a character to work with igraph
    #tagfile            <- tag_list[which(grepl(colony,tag_list))]
    tag <- read.tag(tag_list) #AW
    # if (length(tagfile)>1){
    #   tagfile <- tagfile[grepl(components[grepl("Treatment",components)],tagfile)]
    # }
    # tag                <- read.tag(tagfile)$tag
    #names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
    tag[which(tag$age==0),"age"]   <- NA ###unknown ages are coded as 0 in the tag file
    #tag <-tag[which(tag$final_status=="alive"),] ###remove dead ants from tag file
    ####read interactions
    interactions       <- read.table(network_file,header=T,stringsAsFactors = F)  
    alive <- tag$tag #AW
    #remove dead ants from interactions list #AW
    interactions <- subset(interactions, Tag1 %in% alive)
    interactions <- subset(interactions, Tag2 %in% alive)
    
    # if(any(interactions$Tag2==1)){print(paste(root_name, "has queen ID in Tag2"))} # LS: for "observed" there are never any queen's (1) in Tag2
    #### add a column containing interaction duration in min
    interactions["duration_min"] <- (interactions$Stoptime - interactions$Starttime + (1/FRAME_RATE)) /60 ###duration in minutes (one frame = 0.125 second) #AW
    interactions$N               <- 1 ###duration in minutes (one frame = 0.125 second) #AW
    #### add a column containing the status of tag 1 and the status of tag2
    foragers <- colony_task_group[which(colony_task_group$task_group=="forager"),"tag"]
    nurses   <- colony_task_group[which(colony_task_group$task_group=="nurse"),"tag"]
    
    #####1. calculate within/between caste interactions for pre period
    interactions[c("status_Tag1","status_Tag2")] <- NA
    
    interactions[which(interactions$Tag1%in%foragers),"status_Tag1"] <- "forager"
    interactions[which(interactions$Tag2%in%foragers),"status_Tag2"] <- "forager"
    
    interactions[which(interactions$Tag1%in%nurses),"status_Tag1"] <- "nurse"
    interactions[which(interactions$Tag2%in%nurses),"status_Tag2"] <- "nurse"
    
    interactions[which(interactions$Tag1==queenid),"status_Tag1"] <- "queen"
    interactions[which(interactions$Tag2==queenid),"status_Tag2"] <- "queen"
    
    if (grepl("grooming",input_path)){
      interactions$Actor <- gsub("ant_","",interactions$Act_Name)
      interactions$Receiver <- gsub("ant_","",interactions$Rec_Name)
      
      interactions[c("status_Actor","status_Receiver")] <- NA
      
      interactions[which(interactions$Actor%in%foragers),"status_Actor"] <- "forager"
      interactions[which(interactions$Receiver%in%foragers),"status_Receiver"] <- "forager"
      
      interactions[which(interactions$Actor%in%nurses),"status_Actor"] <- "nurse"
      interactions[which(interactions$Receiver%in%nurses),"status_Receiver"] <- "nurse"
      
      interactions[which(interactions$Actor==queenid),"status_Actor"] <- "queen"
      interactions[which(interactions$Receiver==queenid),"status_Receiver"] <- "queen"
    }
    
    if (period=="pre"){
      inter <- interactions
      inter <- inter[,sort(names(inter))]
      all_interactions <- rbind(all_interactions,data.frame(randy=input_folder,colony_size=colony_size,inter)) # LS: remove period=period because column already exists in "interactions" when we read in the file # although it will not appear twice in the files that are ultimately saved anyway
    }
    
    # #####2. continue calculations for pre vs post
    interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
    interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
    
    
    #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
    aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag1+status_Tag2,FUN=sum,data=interactions)
    names(aggregated1)          <- c("tag","partner_status","duration_min","number_contacts")
    aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag2+status_Tag1,FUN=sum,data=interactions)
    names(aggregated2)          <- c("tag","partner_status","duration_min","number_contacts")
    aggregated                  <- rbind(aggregated1,aggregated2)
    aggregated                  <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,number_contacts)~tag+partner_status,FUN=sum,data=aggregated)
    
    full_table                  <- expand.grid(tag=tag[which(tag$final_status=="alive"),"tag"],partner_status=c("queen","nurse","forager","treated"),stringsAsFactors = F)
    full_table                  <- merge(full_table,aggregated[c("tag","partner_status","duration_min","number_contacts")],all.x=T,all.y=T)
    full_table[is.na(full_table$duration_min),"duration_min"] <- 0
    full_table[is.na(full_table$number_contacts),"number_contacts"] <- 0    
    
    full_table                  <- merge(full_table,tag[c("tag","group")]); names(full_table)[names(full_table)=="group"] <- "status"
    # full_table                  <- merge(full_table,colony_task_group[c("tag","task_group")]); full_table[which(full_table$status=="treated"),"task_group"] <- "treated" # LS: there are no "treated" in status so this just adds the same column as before (="status" with three levels: queen, nurse, forager)
    # LS: fix line above by adding treated to task_group using "colony_treated"
    full_table                  <- merge(full_table,colony_task_group[c("tag","task_group")]); full_table[which(full_table$tag%in%colony_treated), "task_group"] <- "treated"
    
    if (!grepl("age",data_path)){
      full_table$age <- NA
    }else{
      full_table                <- merge(full_table,colony_ages,all.x=T,all.y=F)
    }
    full_table           <- full_table[c("tag","age","task_group","status","partner_status","duration_min","number_contacts")]
    summary_interactions <- rbind(summary_interactions,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,period_detail=period_detail,period_circadian=period_circadian,time_hours=time_hours,time_of_day=time_of_day,full_table,stringsAsFactors = F)) # LS: add period_detail & period_circadian
    
    if (grepl("grooming",input_path)){ # LS: did not check/change this part yet
      #####2. continue calculations for pre vs post
      interactions[which(interactions$Actor%in%colony_treated),"status_Actor"] <- "treated"
      interactions[which(interactions$Receiver%in%colony_treated),"status_Receiver"] <- "treated"
      
      
      #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
      rm(list=ls()[which(grepl("aggregate",ls()))])
      
      full_table                  <- expand.grid(tag=tag[which(tag$final_status=="alive"),"tag"],stringsAsFactors = F)
      
      try(aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Receiver,FUN=sum,data=interactions),silent=T)
      if(exists("aggregated1")){
        names(aggregated1)          <- c("tag","duration_grooming_received_min","number_contacts_received")
        full_table                  <- merge(full_table,aggregated1,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_received_min),"duration_grooming_received_min"] <- 0
        full_table[is.na(full_table$number_contacts_received),"number_contacts_received"] <- 0
      }else{
        full_table$duration_grooming_received_min <- 0 
        full_table$number_contacts_received <- 0
      }
      
      try(aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Actor,FUN=sum,data=interactions[which(interactions$status_Receiver=="treated"),]),silent=T)
      if(exists("aggregated2")){
        names(aggregated2)          <- c("tag","duration_grooming_given_to_treated_min","number_contacts_given_to_treated")
        full_table                  <- merge(full_table,aggregated2,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_given_to_treated_min),"duration_grooming_given_to_treated_min"] <- 0
        full_table[is.na(full_table$number_contacts_given_to_treated),"number_contacts_given_to_treated"] <- 0
      }else{
        full_table$duration_grooming_given_to_treated_min <- 0 
        full_table$number_contacts_given_to_treated <- 0 
      }
      
      
      try(aggregated3                 <- aggregate(na.rm=T,na.action="na.pass",c(duration_min,N)~Receiver,FUN=sum,data=interactions[which(interactions$ant1.zones==1),]),silent=T)
      if(exists("aggregated3")){
        names(aggregated3)          <- c("tag","duration_grooming_received_min_zone1","number_contacts_received_zone1")
        full_table                  <- merge(full_table,aggregated3,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_received_min_zone1),"duration_grooming_received_min_zone1"] <- 0
        full_table[is.na(full_table$number_contacts_received_zone1),"number_contacts_received_zone1"] <- 0
      }else{
        full_table$duration_grooming_received_min_zone1 <- 0 
        full_table$number_contacts_received_zone1 <- 0 
      }
      
      try(aggregated4                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Actor,FUN=sum,data=interactions[which(interactions$ant1.zones==1&interactions$status_Receiver=="treated"),]),silent=T)
      if(exists("aggregated4")){
        names(aggregated4)          <- c("tag","duration_grooming_given_to_treated_min_zone1","number_contacts_given_to_treated_zone1")
        full_table                  <- merge(full_table,aggregated4,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_given_to_treated_min_zone1),"duration_grooming_given_to_treated_min_zone1"] <- 0
        full_table[is.na(full_table$number_contacts_given_to_treated_zone1),"number_contacts_given_to_treated_zone1"] <- 0
      }else{
        full_table$duration_grooming_given_to_treated_min_zone1 <- 0 
        full_table$number_contacts_given_to_treated_zone1 <- 0 
      }
      
      try(aggregated5                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Receiver,FUN=sum,data=interactions[which(interactions$ant1.zones==2),]),silent=T)
      if(exists("aggregated5")){
        names(aggregated5)          <- c("tag","duration_grooming_received_min_zone2","number_contacts_received_zone2")
        full_table                  <- merge(full_table,aggregated5,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_received_min_zone2),"duration_grooming_received_min_zone2"] <- 0
        full_table[is.na(full_table$number_contacts_received_zone2),"number_contacts_received_zone2"] <- 0
      }else{
        full_table$duration_grooming_received_min_zone2 <- 0 
        full_table$number_contacts_received_zone2 <- 0 
      }
      
      try(aggregated6                 <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Actor,FUN=sum,data=interactions[which(interactions$ant1.zones==2&interactions$status_Receiver=="treated"),]),silent=T)
      if(exists("aggregated6")){
        names(aggregated6)          <- c("tag","duration_grooming_given_to_treated_min_zone2","number_contacts_given_to_treated_zone2")
        full_table                  <- merge(full_table,aggregated6,all.x=T,all.y=T)
        full_table[is.na(full_table$duration_grooming_given_to_treated_min_zone2),"duration_grooming_given_to_treated_min_zone2"] <- 0
        full_table[is.na(full_table$number_contacts_given_to_treated_zone2),"number_contacts_given_to_treated_zone2"] <- 0
      }else{
        full_table$duration_grooming_given_to_treated_min_zone2 <- 0 
        full_table$number_contacts_given_to_treated_zone2 <- 0 
      }
      
      full_table                  <- merge(full_table,tag[c("tag","group")]); names(full_table)[names(full_table)=="group"] <- "status"
      full_table                  <- merge(full_table,colony_task_group[c("tag","task_group")]); full_table[which(full_table$status=="treated"),"task_group"] <- "treated"
      if (!grepl("age",data_path)){
        full_table$age <- NA
      }else{
        full_table                <- merge(full_table,colony_ages,all.x=T,all.y=F)
      }
      summary_interactions_grooming <- rbind(summary_interactions_grooming,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_hours=time_hours,time_of_day=time_of_day,full_table,stringsAsFactors = F)) # LS: add period_detail?
      
    }
    
    
    
    clean()
  }
  #####if folder = observed, use the summary_interactions table to compute inter-caste contacts #####
  if (grepl("main",data_path)&input_folder=="observed"){
    summary_interactions <- summary_interactions [which(!summary_interactions$partner_status%in%c("treated","queen")),]
    summary_interactions <- summary_interactions [which(!summary_interactions$task_group%in%c("treated","queen")),]
    
    summary_interactions["within_vs_between"] <- summary_interactions$partner_status==summary_interactions$task_group
    summary_interactions <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,number_contacts)~.,FUN=sum,data=summary_interactions[c("colony","period","period_detail","period_circadian","time_hours","time_of_day","tag","task_group","status","partner_status","duration_min","number_contacts","within_vs_between")]) # LS: add period_detail and period_detail
    summary_interactions_between <- summary_interactions[which(!summary_interactions$within_vs_between),!names(summary_interactions)%in%c("partner_status","within_vs_between")];
    names(summary_interactions_between)[names(summary_interactions_between)=="duration_min"] <- "inter_caste_contact_duration"
    names(summary_interactions_between)[names(summary_interactions_between)=="number_contacts"] <- "inter_caste_contact_number"
    
    #####add information to individual behaviour file
    behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
    behav <- behav[,which(!names(behav)%in%c("inter_caste_contact_duration","inter_caste_contact_number",names(behav)[which(grepl("duration_grooming",names(behav)))]))]
    
    behav <- merge(behav,summary_interactions_between[c("colony","tag","time_hours","inter_caste_contact_duration","inter_caste_contact_number")],all.x=T,all.y=F,sort=F)
    
    if(grepl("grooming",input_path)){
      behav <- merge(behav,summary_interactions_grooming[c("colony","tag","time_hours",names(summary_interactions_grooming)[which(grepl("duration",names(summary_interactions_grooming)))])],all.x=T,all.y=F,sort=F)
      behav <- merge(behav,summary_interactions_grooming[c("colony","tag","time_hours",names(summary_interactions_grooming)[which(grepl("number",names(summary_interactions_grooming)))])],all.x=T,all.y=F,sort=F)
      
    }
    behav <- behav[order(behav$colony,behav$tag,behav$time_hours),]
    
    write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
  }
  
  ###Use all_interactions to obtain information about age-based DoL pre treatment
  ###first check that all_interactions only has Pre-treatment
  all_interactions <- all_interactions[which(all_interactions$period=="pre"),]
  ###summ all interactions for each pair of ants
  all_interactions <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~randy+colony+colony_size+period+period_detail+period_circadian+Tag1+Tag2+status_Tag1+status_Tag2+treatment,FUN=sum,data=all_interactions) # LS: add period_detail & period_circadian
  ###calculate intra_caste_over_inter_caste_WW_contact_duration
  all_interactions["same_caste"] <- all_interactions$status_Tag1==all_interactions$status_Tag2
  same_caste_interactions <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~randy+colony+colony_size+treatment+period+period_detail+period_circadian,FUN=sum,data=all_interactions[which((all_interactions$status_Tag1!="queen")&(all_interactions$status_Tag2!="queen")&(all_interactions$same_caste)),]) # LS: add period_detail & period_circadian
  names(same_caste_interactions)[names(same_caste_interactions)=="duration_min"] <- "duration_min_within"
  names(same_caste_interactions)[names(same_caste_interactions)=="N"] <- "number_contacts_within"
  inter_caste_interactions <- aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~randy+colony+colony_size+treatment+period+period_detail+period_circadian,FUN=sum,data=all_interactions[which((all_interactions$status_Tag1!="queen")&(all_interactions$status_Tag2!="queen")&(!all_interactions$same_caste)),]) # LS: add period_detail & period_circadian
  names(inter_caste_interactions)[names(inter_caste_interactions)=="duration_min"] <- "duration_min_between"
  names(inter_caste_interactions)[names(inter_caste_interactions)=="N"] <- "number_contacts_between"
  inter_intra_caste_interactions <- merge(same_caste_interactions,inter_caste_interactions,all.x=T,all.y=T)
  inter_intra_caste_interactions["intra_caste_over_inter_caste_WW_contact_duration"] <- inter_intra_caste_interactions$duration_min_within/inter_intra_caste_interactions$duration_min_between
  inter_intra_caste_interactions["intra_caste_over_inter_caste_WW_contact_number"] <- inter_intra_caste_interactions$number_contacts_within/inter_intra_caste_interactions$number_contacts_between
  dol <-   inter_intra_caste_interactions[!names(inter_intra_caste_interactions)%in%c("duration_min_within","duration_min_between","number_contacts_between","number_contacts_within")]
  
  ###calculate queen contact with nurses vs. workers
  ## skip queen in grooming interactions
  if (!grepl("grooming",input_path)) { #AW
    queen_interactions <- all_interactions[which(all_interactions$status_Tag1=="queen"|all_interactions$status_Tag2=="queen"),]
    queen_interactions[which(queen_interactions$status_Tag1=="queen"),"partner"] <- queen_interactions[which(queen_interactions$status_Tag1=="queen"),"Tag2"]
    queen_interactions[which(queen_interactions$status_Tag1=="queen"),"partner_status"] <- queen_interactions[which(queen_interactions$status_Tag1=="queen"),"status_Tag2"]
    queen_interactions[which(queen_interactions$status_Tag2=="queen"),"partner"] <- queen_interactions[which(queen_interactions$status_Tag2=="queen"),"Tag1"]
    queen_interactions[which(queen_interactions$status_Tag2=="queen"),"partner_status"] <- queen_interactions[which(queen_interactions$status_Tag2=="queen"),"status_Tag1"]
    
    interaction_with_nurses <-aggregate (na.rm=T,na.action="na.pass",cbind(duration_min,N)~randy+colony+colony_size+period+period_detail+period_circadian+treatment,FUN=sum,data=queen_interactions[which(queen_interactions$partner_status=="nurse"),]) # LS: add period_detail & period_circadian
    names(interaction_with_nurses)[names(interaction_with_nurses)=="duration_min"] <- "duration_min_with_nurses"
    names(interaction_with_nurses)[names(interaction_with_nurses)=="N"] <- "number_contacts_with_nurses"
    interaction_with_forager <-aggregate (na.rm=T,na.action="na.pass",cbind(duration_min,N)~randy+colony+colony_size+period+period_detail+period_circadian+treatment,FUN=sum,data=queen_interactions[which(queen_interactions$partner_status=="forager"),]) # LS: add period_detail & period_circadian
    names(interaction_with_forager)[names(interaction_with_forager)=="duration_min"] <- "duration_min_with_foragers"
    names(interaction_with_forager)[names(interaction_with_forager)=="N"] <- "number_contacts_with_foragers"
    queen_interac <- merge(interaction_with_nurses,interaction_with_forager,all.x=T,all.y=T)
    queen_interac["QNurse_over_QForager_contact_duration"] <- queen_interac$duration_min_with_nurses/queen_interac$duration_min_with_foragers
    queen_interac["QNurse_over_QForager_contact_number"] <- queen_interac$number_contacts_with_nurses/queen_interac$number_contacts_with_foragers
    dol <- merge(dol,queen_interac[c("randy","colony","period","period_detail","period_circadian","QNurse_over_QForager_contact_duration","QNurse_over_QForager_contact_number")])  # LS: add period_detail & period_circadian
  }
  
  ###if necessary: add age
  if (grepl("age",data_path)){
    partner_ages <- ages; names(partner_ages) <- c("colony","partner","partner_age")
    queen_interactions <- merge(queen_interactions,partner_ages,all.x=T,all.y=F)
    
    ages_Tag1                    <- ages; names(ages_Tag1) <- c("colony","Tag1","age_Tag1")
    ages_Tag2                    <- ages; names(ages_Tag2) <- c("colony","Tag2","age_Tag2")
    all_interactions             <- merge(all_interactions,ages_Tag1,all.x=T,all.y=F)
    all_interactions             <- merge(all_interactions,ages_Tag2,all.x=T,all.y=F)
    all_interactions["age_diff"] <- abs(all_interactions$age_Tag2-all_interactions$age_Tag1)
    
    ###write ordered pair of interacting ants
    all_interactions["ordered"]<- all_interactions[,"Tag1"]<all_interactions[,"Tag2"]
    all_interactions[which(all_interactions$ordered),"new_Tag1"] <-  all_interactions[which(all_interactions$ordered),"Tag1"]
    all_interactions[which(all_interactions$ordered),"new_Tag2"] <-  all_interactions[which(all_interactions$ordered),"Tag2"]
    all_interactions[which(!all_interactions$ordered),"new_Tag1"] <-  all_interactions[which(!all_interactions$ordered),"Tag2"]
    all_interactions[which(!all_interactions$ordered),"new_Tag2"] <-  all_interactions[which(!all_interactions$ordered),"Tag1"]
    all_interactions["pair"]   <- as.character(interaction(all_interactions$colony,all_interactions$new_Tag1,all_interactions$new_Tag2))
    ###to get slope: get each pair of possibly interacting ants
    for (colony in unique(all_interactions$colony)){
      colony_ages <- ages[which(ages$colony==colony),]
      colony_ages_Tag1                    <- colony_ages; names(colony_ages_Tag1) <- c("colony","Tag1","age_Tag1")
      colony_ages_Tag2                    <- colony_ages; names(colony_ages_Tag2) <- c("colony","Tag2","age_Tag2")
      tagfile            <- tag_list[which(grepl(colony,tag_list))][1]
      tag <- read.tag(tagfile)$tag
      names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
      
      ####list ants of known ages
      known_ages                <- colony_ages[which(!is.na(colony_ages$age)),"tag"]
      ###remove queen from list
      known_ages                <- known_ages[known_ages!=queenid]
      ###remove dead ants from list
      known_ages                <- known_ages[which(known_ages%in%tag[which(tag$final_status=="alive"),"tag"])]
      
      ### prepare complete pair table
      summ_pairs           <- data.frame(combinations(n=length(known_ages),r=2,v=known_ages))
      summ_pairs["ordered"]<- summ_pairs[,1]<summ_pairs[,2]
      summ_pairs[which(summ_pairs$ordered),"Tag1"] <-  summ_pairs[which(summ_pairs$ordered),1]
      summ_pairs[which(summ_pairs$ordered),"Tag2"] <-  summ_pairs[which(summ_pairs$ordered),2]
      summ_pairs[which(!summ_pairs$ordered),"Tag1"] <-  summ_pairs[which(!summ_pairs$ordered),2]
      summ_pairs[which(!summ_pairs$ordered),"Tag2"] <-  summ_pairs[which(!summ_pairs$ordered),1]
      summ_pairs              <- merge(summ_pairs,colony_ages_Tag1,all.x=T,all.y=F)
      summ_pairs              <- merge(summ_pairs,colony_ages_Tag2,all.x=T,all.y=F)
      summ_pairs["age_diff"]  <- abs(summ_pairs$age_Tag2-summ_pairs$age_Tag1)
      summ_pairs["pair"]   <- as.character(interaction(colony,summ_pairs$Tag1,summ_pairs$Tag2))
      summ_pairs <- merge(summ_pairs,unique(all_interactions[which(all_interactions$colony==colony),c("colony","randy","colony_size","period","treatment")])) # LS: add period_detail?
      ###merge it with interactions whose age diff is known 
      summ_pairs <- merge(summ_pairs,all_interactions[which(all_interactions$colony==colony),c("colony","randy","colony_size","period","treatment","pair","duration_min","N")],all.x=T,all.y=F) # LS: add period_detail?
      ###and fill intreactions which did not happen with 0
      summ_pairs[which(is.na(summ_pairs$duration_min)),"duration_min"] <- 0
      summ_pairs[which(is.na(summ_pairs$N)),"N"] <- 0
      model_WW <- lm(duration_min~age_diff,data=summ_pairs)
      dol[dol$colony==colony,"slope_WW_contact_duration_f_age_diff"] <- coef(model_WW)["age_diff"]
      model_WW <- lm(N~age_diff,data=summ_pairs)
      dol[dol$colony==colony,"slope_WW_contact_number_f_age_diff"] <- coef(model_WW)["age_diff"]
      
      ###do the same for queen interactions
      colony_ages <- colony_ages[which(!is.na(colony_ages$age)),]
      names(colony_ages) <-c("colony","partner","age")
      queen_W <- merge(colony_ages,unique(queen_interactions[which(queen_interactions$colony==colony),c("colony","randy","colony_size","period","treatment")])) # LS: add period_detail?
      ###merge it with interactions whose age diff is known 
      queen_W <- merge(queen_W,queen_interactions[which(queen_interactions$colony==colony),c("colony","randy","colony_size","period","treatment","partner","duration_min","N")],all.x=T,all.y=F) # LS: add period_detail?
      ###and fill intreactions which did not happen with 0
      queen_W[which(is.na(queen_W$duration_min)),"duration_min"] <- 0
      queen_W[which(is.na(queen_W$N)),"N"] <- 0
      model_QW <- lm(duration_min~age,data=queen_W)
      dol[dol$colony==colony,"slope_QW_contact_duration_f_W_age"] <- coef(model_QW)["age"]
      model_QW <- lm(N~age,data=queen_W)
      dol[dol$colony==colony,"slope_QW_contact_number_f_W_age"] <- coef(model_QW)["age"]
      
    }
  }
  summary_dol <- rbind(summary_dol,dol)  
  
}
if(!file.exists(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed",sep=""))){dir.create(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed",sep=""),recursive=T)}
if(!file.exists(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed/interactions.dat",sep=""))){write.table(summary_dol,file=paste(data_path,"/processed_data/collective_behaviour/random_vs_observed/interactions.dat",sep=""),col.names=T,row.names=F,append=F,quote=F)}
to_keep <- to_keep_ori
