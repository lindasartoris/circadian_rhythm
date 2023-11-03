####13_network_analysis.R#####

#### Takes an interaction list as an input, builds a network, and analyse its properties

###Created by Nathalie Stroeymeyt
###Modified by Adriano Wanderlingh to work with FORT formicidae Tracking data. Mods tagged with the comment "AW". script wide mods cited here below.
###Modified by Nathalie Stroeymeyt to include number of events in addition to duration
###with adaptations by Linda Sartoris


#Script wide mods AW
# - replaced before/after with pre/post

####################################
to_keep_ori <- to_keep
# library(outliers)
#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)
edge_weights <-"number"  ### "duration" or "number" # LS: run once for duration of interactions and once for number of interactions

###remove output file 
if (file.exists(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))){
  file.remove(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))
}

#### get input file list
if (!grepl("survival",data_path)){
  input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
  setwd(input_path)  
  input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
  input_folders        <- input_folders[which(input_folders!="")]
}else{
  input_path           <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
  setwd(input_path)  
  input_folders        <- list.dirs(recursive=T,path="PostTreatment",full.names=F)
  input_folders        <- input_folders[which(input_folders!="")]
}

queen_community_summary <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_files","options","option","summary_collective","summary_individual","outputfolder","network_file","queenid", "edge_weights")
for (input_folder in input_folders){
  print(input_folder)
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  if (input_folder=="observed"){
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
    if(grepl("main",data_path)){
      options <- c("all_workers","untreated_only")
    }else{
      options <- c("all_workers")
    }
  }else{
    options <- c("all_workers")
  }
  for (option in options){
    print(option)
    outputfolder <- paste(data_path,"/processed_data/network_properties_edge_weights_",edge_weights,sep="")
    
    summary_collective <- NULL
    summary_individual <- NULL
    for (network_file in network_files){
      # network_file <- "PostTreatment/observed/colony05BP_pathogen.big_PostTreatment_TH6_TD18_interactions.txt"

      # #### TEMPORARY EXCLUSION OF "PreTreatment/observed/colony08BP_pathogen.big_PreTreatment_TH-12_TD0_interactions.txt" as:
      # ## Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate
      #   if (any(grepl("colony08BP|colony06SP|colony02SP", network_file))) {
      #   warning(paste(input_folder, "has too few interactions to aggregate the cumulated duration of interaction with treated workers
      #                 \nor the queen does not belong to any community"))
      # }else{
      #  has too few interactions to aggregate the cumulated duration of interaction with treated workers: "PreTreatment/observed/colony08BP","PreTreatment/observed/colony06SP"
      #  the queen does not belong to any community: "PreTreatment/observed/colony08BP" , PostTreatment/observed/colony08BS
      
      cat("\r",network_file)
      ####get file metadata
      root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))]) # LS: replace grepl("colony", ...) with grepl("interactions")
      components         <- unlist(strsplit(root_name,split="_"))
      # colony             <- components[grepl("colony",components)]
      colony             <- unlist(strsplit(root_name,split="_"))[1] # LS
      treatment          <- info[which(info$colony==colony),"treatment"] #AW: no need for as.numeric() 
      treatment_circadian<- unlist(strsplit(unlist(strsplit(root_name,split="_"))[2],split="\\."))[2] # LS
      colony_size        <- info[which(info$colony==colony),"colony_size"]
      
      if (!all(!grepl("PreTreatment",components))){period <- "pre"}else{period <- "post"} 
      
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
      
      if(!grepl("survival",data_path)){
        time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
        time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
      }else{ # LS: would have to change but probably never needed bc I don't have survival experiment. Different depending on Day or Night, pre1 or pre2
        if (period=="post"){
          time_hours   <- 0
          time_of_day <- 12 
        }else{
          time_hours   <- -27 #AW: shifted time window by 3h
          time_of_day <- 9 #AW: shifted time window by 3h
        }
      }
      
      ####get appropriate task_group list, treated list and tag
      colony_treated     <- treated[which(treated$colony==colony),"tag"] #AW
      colony_task_group  <- task_groups[which(task_groups$colony==colony),]
      queenid            <- as.character(colony_task_group[which(colony_task_group$task_group=="queen"),"tag"]) #AW: call specific queen tag instead of fixed 665. it has to be a character to work with igraph
      #tagfile            <- tag_list[which(grepl(colony,tag_list))]
      # if (length(tagfile)>1){
      #   tagfile <- tagfile[grepl(unlist(strsplit(gsub("\\.txt","",root_name),"_"))[grepl("Treatment",unlist(strsplit(gsub("\\.txt","",root_name),"_")))],tagfile)]
      # }
      
      ####read interactions
      interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
      tag <- read.tag(tag_list)
      #tag                <- read.tag(tagfile)$tag # FLAG
      #names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),] #AW
      #tag[which(tag$age==0),"age"]   <- NA ###unknown ages are coded as 0 in the tag file #AW
      # if (!grepl("survival",data_path)){
      #   alive      <- tag[which(tag$final_status=="alive"),"tag"]
      # }else{
      #   alive      <- tag[which(as.numeric(tag$death)==0|as.numeric(tag$death)>=max(interactions$Stopframe,na.rm=T)),"tag"]
      # }
      #tag <-tag[which(tag$tag%in%alive),] ###remove dead ants from tag file
      alive <- tag$tag #AW
      #remove dead ants from interactions list  #AW
      interactions <- subset(interactions, Tag1 %in% alive)
      interactions <- subset(interactions, Tag2 %in% alive)

      ####if untreated only, reduce tag , colony_task_group, and interactions
      if (option=="untreated_only"){
        colony_task_group <- colony_task_group[which(!colony_task_group$tag%in%colony_treated),]
        tag               <- tag[which(!tag$tag%in%colony_treated),]
        interactions      <- interactions[which(!((interactions$Tag1%in%colony_treated)|(interactions$Tag2%in%colony_treated))),]
      }
      
      actors             <- data.frame(name=as.character(tag$tag))
      
      #### add a column containing interaction duration in min
      interactions["duration_min"] <- (interactions$Stoptime - interactions$Starttime + (1/FRAME_RATE)) /60 ###duration in minutes (one frame = 0.125 second) #AW
      interactions$N               <- 1 ###each interaction is given a count of 1
      
      #### add a column containing the status of tag 1 and the status of tag2
      interactions[c("status_Tag1","status_Tag2")] <- "untreated"
      interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
      interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
      
      #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
      if (input_folder=="observed"&option=="all_workers"){
        
        ### AW: To overcome the "Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate", create empty df if no rows to aggregate on
        aggregated1 <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag1+status_Tag2,FUN=sum,data=interactions[which(interactions$status_Tag2=="treated"),]) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(), number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag1+status_Tag2,FUN=sum,data=interactions[which(interactions$status_Tag2=="treated"),])
        names(aggregated1)          <- c("tag","partner_status","duration_min","number_contacts")
        aggregated2 <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,N)~Tag2+status_Tag1,FUN=sum,data=interactions[which(interactions$status_Tag1=="treated"),]) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(), number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag2+status_Tag1,FUN=sum,data=interactions[which(interactions$status_Tag1=="treated"),])
        names(aggregated2)          <- c("tag","partner_status","duration_min","number_contacts")
        aggregated                  <- rbind(aggregated1,aggregated2)
        aggregated <- tryCatch(
          { aggregate(na.rm=T,na.action="na.pass",cbind(duration_min,number_contacts)~tag+partner_status,FUN=sum,data=aggregated) },
          error = function(e) {data.frame(tag= integer(),partner_status=character(), duration_min=numeric(),number_contacts=numeric())} # Return an empty data frame in case of an error
        )
        #aggregated                  <- aggregate(na.rm=T,na.action="na.pass",duration_min~tag+partner_status,FUN=sum,data=aggregated)
        interactions_with_treated   <- merge(data.frame(tag=tag[which(tag$tag%in%alive),"tag"],stringsAsFactors = F),aggregated[c("tag","duration_min","number_contacts")],all.x=T)
        interactions_with_treated[is.na(interactions_with_treated$duration_min),"duration_min"] <- 0
        interactions_with_treated[is.na(interactions_with_treated$number_contacts),"number_contacts"] <- 0
        names(interactions_with_treated) <- c("tag","duration_of_contact_with_treated_min","number_of_contact_with_treated")
        interactions_with_treated["colony"] <- colony
        interactions_with_treated["time_hours"] <- time_hours
        ###write results
        if (grepl("main",data_path)){
          behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
          if (!"duration_of_contact_with_treated_min"%in%names(behav)){
            behav[c("duration_of_contact_with_treated_min")] <- NA
            behav[c("number_of_contact_with_treated")] <- NA
          }
          behav[match(as.character(interaction(interactions_with_treated$colony,interactions_with_treated$tag,interactions_with_treated$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("duration_of_contact_with_treated_min")]  <- interactions_with_treated$duration_of_contact_with_treated_min
          behav[match(as.character(interaction(interactions_with_treated$colony,interactions_with_treated$tag,interactions_with_treated$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("number_of_contact_with_treated")]        <- interactions_with_treated$number_of_contact_with_treated
          options(digits=3)
          write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
          options(digits=16)
        }
        ###interactions with treated: post-treatment
        if (period=="post"){
          outputfoldy <- paste(data_path,"/processed_data/individual_behaviour/post_treatment",sep="")
          if(!file.exists(outputfoldy)){dir.create(outputfoldy,recursive=T)}
          int_with_treated <- data.frame(colony_size=colony_size,treatment=treatment,period=period, period_detail=period_detail, period_circadian=period_circadian, time_of_day=time_of_day,# LS: add period_detail & period_circadian here to keep the same columns for pre and post
                                         interactions_with_treated,stringsAsFactors = F)
          if (!file.exists(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))){
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
          }else{
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=F,row.names=F,quote=F,append=T)
          }
        }
      }
      
      ###build NETWORK
      if (!grepl("survival",data_path)){
        net <- graph.data.frame(interactions[c("Tag1","Tag2")],directed=F,vertices=actors)
        ### add edge weights
        if (edge_weights=="number"){
          E(net)$weight <- interactions[,"N"]
        }else if (edge_weights=="duration"){
          E(net)$weight <- interactions[,"duration_min"]        
         }

        ###simplify graph (merge all edges involving the same pair of ants into a single one whose weight = sum of these weights)
        net <- simplify(net,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")
        
        ###get all degress without removing any node
        degrees_all         <- degree(net,mode="all")
        ##################remove unconnected nodes
        unconnected <- actors[degree(net)==0,]
        net <- net - as.character(unconnected)
        ##get actor list from net
        actors <- get.vertex.attribute(net,"name")
        
      ####Part 1: individual network properties ####
        ###prepare table
        tag["status"] <- "untreated"; tag[which(tag$tag%in%colony_treated),"status"] <- "treated"
        individual <- data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,tag=tag$tag,age=tag$age,status=tag$group,period=period,period_detail=period_detail,period_circadian=period_circadian,time_hours=time_hours,time_of_day=time_of_day, # LS: add period_detail here
                                 degree=NA,
                                 aggregated_distance_to_queen=NA,
                                 mean_aggregated_distance_to_treated=NA,
                                 same_community_as_queen=NA)
        ##degree
         individual[match(names(degrees_all),individual$tag),"degree"] <- degrees_all
        
        ## skip queen in grooming interactions
        if (!grepl("grooming",input_path)) {
          communities             <- cluster_louvain(net, weights = E(net)$weight)
          community_membership    <- communities$membership
          ##same community as queen
          queen_comm <- community_membership[which(V(net)$name==queenid)]
          community_membership <- community_membership==queen_comm
          individual[match(V(net)$name,individual$tag),"same_community_as_queen"] <- community_membership
          ##path length to queen
          if (queenid%in%actors){
            path_length_to_queen <- t(shortest.paths(net,v=actors,to=queenid,weights=1/E(net)$weight))
            individual[match(colnames(path_length_to_queen),individual$tag),"aggregated_distance_to_queen"] <- as.numeric(path_length_to_queen )
          }
        }
        
        ########Mean path length to treated; aggregated_network
        if(option!="untreated_only"){
          path_length_to_treated                             <- as.data.frame(as.matrix(shortest.paths(net,v=actors,to=as.character(colony_treated)[as.character(colony_treated)%in%V(net)$name],weights=1/E(net)$weight)))
          path_length_to_treated["mean_distance_to_treated"] <- NA
          path_length_to_treated$mean_distance_to_treated    <- as.numeric(rowMeans(path_length_to_treated,na.rm=T))
          individual[match(rownames(path_length_to_treated),individual$tag),"mean_aggregated_distance_to_treated"] <- path_length_to_treated[,"mean_distance_to_treated"]
        }
        ###Add data to main data table
        summary_individual <- rbind(summary_individual,individual)
        
        
        ####Part 2: collective network properties ####
        ##################remove outliers (e.g. ants that have only 1 or 2 connections)
        # outlier_p <- 0
        # outlier_removed <- c()
        # while(outlier_p<0.05){
        #   outlier_p <- grubbs.test(sort(degree(net))[1:min(30,length(degree(net)))])$p.value ###do the test on lowest 30 values first as test does not work as well when there are too many nodes
        #   if (outlier_p<0.05){
        #     outlier_removed <- c(outlier_removed,actors[which.min(degree(net))])
        #     net <- net - as.character(actors[which.min(degree(net))])
        #   }
        # }
        outlier_removed <- actors[which(strength(net)<(5/60))]
        net <- net - as.character(outlier_removed)
        
        ##################update actor list
        actors <- get.vertex.attribute(net,"name")
        
        ##Assortativity  - Age
        ####if age experiment, get colony ages
        if (grepl("age",data_path)){
          colony_ages <- ages [which(ages$colony==colony),]
          ####set queen age to NA as this would bias the result (the queen is the oldest individual and interacts mostly with the young nurses)
          colony_ages[which(colony_ages$tag==queenid),"age"] <- NA
          ####order the age according to the order of the network's vertices
          ordered_ages        <- colony_ages[match(V(net)$name,as.character(colony_ages$tag)),"age"]
          #### calculate age assortativity
          age_assortativity <- assortativity(net-V(net)$name[is.na(ordered_ages)],types1=ordered_ages[!is.na(ordered_ages)],directed=F)
        }else{
          age_assortativity <- NA
        }
        ## Assortativity - Task
        ordered_task_groups <- colony_task_group[match(actors,as.character(colony_task_group$tag)),"task_group"]
        ordered_task_groups <- as.numeric(as.factor(ordered_task_groups))
        task_assortativity  <- assortativity_nominal(net,types=ordered_task_groups,directed=F)
        ##Clustering
        clustering <- mean(transitivity(net,type="barrat",weights=E(net)$weight,isolates = c("NaN")),na.rm=T)
        ##Degree mean and max
        degrees         <- degree(net,mode="all")
        degree_mean     <- mean(degrees,na.rm=T)
        degree_maximum  <- max(degrees,na.rm=T)
        ##Density
        density  <- igraph::edge_density(net)
        ##Diameter
        diameter <- igraph::diameter(net,directed=F,unconnected=TRUE,weights=(1/E(net)$weight)) ###here use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        ##Efficiency
        net_dist                    <- shortest.paths(net, weights=1/E(net)$weight, mode="all") ##again use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        net_dist[net_dist==0]       <- NA ##remove distances to self
        efficiency                  <- 1/net_dist ##transform each distance into an efficiency
        efficiency <- (1/((vcount(net)*(vcount(net)-1))))*(sum(efficiency,na.rm=TRUE))
        ## Modularity
        communities             <- cluster_louvain(net, weights = E(net)$weight)
        community_membership    <- communities$membership
        modularity              <- modularity(net,community_membership,weights=E(net)$weight)
        
        ###Add to data
        summary_collective <- rbind(summary_collective,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,period_detail=period_detail,period_circadian=period_circadian,time_hours=time_hours,time_of_day=time_of_day, # LS: add period_detail & period_circadian
                                                                  age_assortativity=age_assortativity,
                                                                  task_assortativity=task_assortativity,
                                                                  clustering=clustering,
                                                                  degree_mean=degree_mean,
                                                                  degree_maximum=degree_maximum,
                                                                  density=density,
                                                                  diameter=diameter,
                                                                  efficiency=efficiency,
                                                                  modularity=modularity,
                                                                  nb_unconnected=length(unconnected),
                                                                  nb_outliers_removed=length(outlier_removed),
                                                                  stringsAsFactors = F))
        
       }
      clean()
    #}# TEMPORARY EXCLUSION 
    }
    
    #print progress AW
    print(" End of network_files processing >> writing")
    
    #####write #####
    if (!grepl("survival",data_path)){
      if (input_folder=="observed"){
        ####Main experiment: write pre_vs_post_treatment data
        if (!grepl("age",data_path)){
          outputfolder2 <- paste(outputfolder,"pre_vs_post_treatment",option,sep="/")
          if(!file.exists(outputfolder2)){dir.create(outputfolder2,recursive=T)}
          write.table(summary_collective[,names(summary_collective)!="randy"],file=paste(outputfolder2,"/colony_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[,names(summary_individual)!="randy"],file=paste(outputfolder2,"/individual_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        }
        ####All workers: write pre_treatment data into random_vs_observed folder
        if (option=="all_workers"){
          outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
          if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)

          outputfolder4 <- paste(outputfolder,"post_treatment",sep="/")
          if(!file.exists(outputfolder4)){dir.create(outputfolder4,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="post"),],file=paste(outputfolder4,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) 
          write.table(summary_individual[which(summary_individual$period=="post"),],file=paste(outputfolder4,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) 
        }
        ######Main experiment, All workers: add pre_treatment node_properties information to pre_treatment behaviour file
        if (!grepl("age",data_path)&option=="all_workers"){
          pre_treatment_behav_file <- paste(data_path,"/processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat",sep="")
          pre_treatment_behav      <- read.table(pre_treatment_behav_file,header=T,stringsAsFactors = F)
          names(summary_individual)[which(names(summary_individual)=="aggregated_distance_to_queen")] <- paste("aggregated_distance_to_queen_edge_weights_",edge_weights,sep="")
          pre_treatment_behav      <- merge(pre_treatment_behav,summary_individual[which(summary_individual$period=="pre"),c("colony","tag","time_hours","degree",paste("aggregated_distance_to_queen_edge_weights_",edge_weights,sep=""))],all.x=T,all.y=T) 
          pre_treatment_behav      <- pre_treatment_behav[order(pre_treatment_behav$colony,pre_treatment_behav$tag,pre_treatment_behav$time_hours),]
          write.table(pre_treatment_behav, file=pre_treatment_behav_file,col.names=T,row.names=F,quote=F,append=F) # LS: this file adds columns to "network_position_vs_time_outside.dat" which already has the information of period_detail in it
        }
      }else{
        outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
        if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
        write.table(summary_collective[which(summary_collective$period=="pre"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) #
        write.table(summary_individual[which(summary_individual$period=="pre"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F) #
        
      }
    }
    ###Get characteristics of queen community vs. other communities (worker age, prop. of foragers)
    if (!grepl("survival",data_path)&option=="all_workers"){ 
      #print progress AW
      print(" End of writing >> Get characteristics of queen community vs. other communities")
      
      summary_individual_pre <- read.table(paste(outputfolder,"/random_vs_observed/node_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F) 
      summary_individual_pre <- summary_individual_pre[which(summary_individual_pre$period=="pre"),] 
      
      ####if necessary: add age
      if (grepl("age",data_path)){
        summary_individual_pre <- merge(summary_individual_pre[,which(names(summary_individual_pre)!="age")],ages,all.x=T,all.y=F)
      }else{
        summary_individual_pre$age <- NA
      }
      ####add task_group
      summary_individual_pre <- merge(summary_individual_pre,task_groups,all.x=T,all.y=F)
      ###remove queen
      summary_individual_pre <- summary_individual_pre[which(summary_individual_pre$tag!=queenid),]
      
      ###1. calculate mean proportion of foragers depending on with vs without queen
      summary_individual_pre["forager"] <- 0
      summary_individual_pre[which(summary_individual_pre$task_group=="forager"),"forager"] <- 1
      ## skip queen in grooming interactions
      if (!grepl("grooming",input_path)) { #AW
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+period_detail+period_circadian+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_pre) # LS: period_detail & period_circadian
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+period_detail+period_circadian+same_community_as_queen,FUN=mean,data=prop_foragers) # LS: period_detail & period_circadian
      names(prop_foragers)[names(prop_foragers)=="same_community_as_queen"] <- "in_queen_comm";names(prop_foragers)[names(prop_foragers)=="forager"] <- "proportion_of_foragers"
      }else{
        prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+period_detail+period_circadian+time_hours,FUN=mean,data=summary_individual_pre) # LS: period_detail & period_circadian
        prop_foragers$in_queen_comm <- NA
        names(prop_foragers)[names(prop_foragers)=="forager"] <- "proportion_of_foragers"
      }
      ###2. calculate mean age of workers depending on with vs without queen
      if (grepl("age",data_path)){
        mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+period_detail+period_circadian+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_pre) # LS: period_detail & period_circadian
        mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+period_detail+period_circadian+same_community_as_queen,FUN=mean,data=mean_age) # LS: period_detail & period_circadian
        names(mean_age)[names(mean_age)=="same_community_as_queen"] <- "in_queen_comm"
        prop_foragers <- merge(prop_foragers,mean_age,all.x=T)
      }else{
        prop_foragers$age <- NA
      }
      
      
      
      
      prop_foragers[which(prop_foragers$in_queen_comm=="FALSE"),"in_queen_comm"] <- "not_with_queen"
      prop_foragers[which(prop_foragers$in_queen_comm=="TRUE"),"in_queen_comm"] <- "with_queen"
      if (grepl("random",input_folder)){
        prop_foragers["randy"] <- "random"
      }
      queen_community_summary <- rbind(queen_community_summary,prop_foragers)
    }
  } 
}


if (!grepl("survival",data_path)){
  queen_community_summary <- aggregate(na.rm=T,na.action="na.pass",cbind(proportion_of_foragers,age)~.,FUN=mean,data=queen_community_summary)
  queen_community_summary <- queen_community_summary[order(queen_community_summary$randy,queen_community_summary$colony),]
  queen_community_summary$treatment <- queen_community_summary$randy
  if (!file.exists(paste(outputfolder,"/random_vs_observed",sep=""))){dir.create(paste(data_path,"/processed_data/network_properties/random_vs_observed",sep=""),recursive=T)}
  write.table(queen_community_summary,file=paste(outputfolder,"/random_vs_observed/queen_community.dat",sep=""),append=F,quote=F,row.names=F,col.names=T)
}
to_keep <- to_keep_ori
