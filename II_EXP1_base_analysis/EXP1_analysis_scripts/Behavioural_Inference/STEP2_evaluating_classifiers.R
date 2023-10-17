rm(list=ls())

####libraries
library(lme4)
library(car)
library(multcomp)

####user and folder
USER <- "Nathalie"
if (USER=="Adriano") {
  SAVEOUTPUT <- "/home/cf19810/Documents"
}
if (USER=="Nathalie"){
  SAVEOUTPUT <- "/media/bzniks/FiveTB"
}

####parameters
measure         <- "Fbeta" ##""Fbeta" "CSI"                               ####which measure to do the selection on
priority        <- "trends" ##~trends-overall
beta            <- 1                                                    ####For F score: which beta to use (CURRENT CHOICE IS BETWEEN 0.5 AND 1)
second_measure  <- c("CSI","Fbeta")[which(c("CSI","Fbeta")!=measure)]     ####which second measure to use to break ties

###read quality score files
output_names <- list.files(path= file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL"),pattern="quality_scores",full.names = T)
output_names <-output_names[which(!grepl("CHOSEN",output_names))]
outcomes <- NULL
for (output_name in output_names){
  outcomes <- rbind(outcomes,read.table(output_name,header=T,stringsAsFactors = F))
}

####keep only the lines with the correct beta
outcomes <- outcomes[which(outcomes$beta==beta),]

if (priority == "overall"){
 outcomes_to_keep <- outcomes[which.max(outcomes[,paste(measure,"_test",sep="")]),]
}else{
  outcomes_scores <- outcomes[which(names(outcomes) %in% names(outcomes)[which(grepl("test",names(outcomes))| grepl("training",names(outcomes)))])]
  outcomes_param  <- outcomes[which(!names(outcomes) %in% c(names(outcomes)[which(grepl("test",names(outcomes))| grepl("pre_classifier",names(outcomes))| grepl("training",names(outcomes)))],"proportion_truegrooming_detected","Loop_ID"))]
  
  outcomes_to_keep <- outcomes
  
  
  for (parameter in c("CAPSULE_FILE","MAX_INTERACTION_GAP",names(outcomes_param)[which(!names(outcomes_param)%in% c("CAPSULE_FILE","MAX_INTERACTION_GAP"))])){
    print(paste("Evaluating best values for parameter",parameter))
    if (length(unique(outcomes_param[,parameter]))>1){
      for (what in c("training","test")){ ####Best parameter values for training and for test, on average for all attempts
        ###Build data frame for stats
        dat <- data.frame(
          score= outcomes_scores[,paste(measure,"_",what,sep="")],
          parameter = as.factor(outcomes_param[,parameter]),
          other_params = interaction(outcomes_param[names(outcomes_param)!=parameter])
        )
        ###Fit stats model
        model <- lmer (score~ parameter + (1|other_params),data=dat)
        
        ###List parameter values in decreasing beta_coefficient order (higher coefficient in model = better)
        coefs         <- data.frame(summary(model)$coefficients)
        best_to_worse <-   c( levels(dat$parameter)[1], gsub("parameter","",rownames(coefs)[2:nrow(coefs)]   )  )[ order(-c(0,coefs$Estimate[2:nrow(coefs)]))]
        
        ###perform Tukey post-hoc comparisons to find out which parameter values are significantly different from which, and order from best to lest good parameter value
        Tukey_cld <- cld(summary(glht(model, linfct = mcp(parameter = "Tukey"))))$mcletters$Letters[best_to_worse]
        print(Tukey_cld)
        ###keep candidates: those that share the first letter in Tukey_cld
        assign ( paste("list_",what,sep="") ,  names(Tukey_cld[which(grepl(Tukey_cld[1],Tukey_cld))]))
        rm(list=c("Tukey_cld","model","coefs","best_to_worse"))
      }
      ###then keep the values that are best for both training and test
      if (!all(!list_training%in%list_test)){
        assign(paste(parameter,"_list",sep=""),     list_training[which(list_training %in% list_test)])
      }else{
        assign(paste(parameter,"_list",sep=""),     list_training)
      }
      
      ###reduce outcomes table
      if (parameter %in% c("CAPSULE_FILE","MAX_INTERACTION_GAP")){
        outcomes_scores <- outcomes_scores[which(outcomes_param[,parameter]%in%get(paste(parameter,"_list",sep=""))),]
        outcomes_param  <- outcomes_param[which(outcomes_param[,parameter]%in%get(paste(parameter,"_list",sep=""))),]
      }
      rm(list=c("list_test","list_training"))
    }else{
      assign(paste(parameter,"_list",sep=""),unique(outcomes_param[,parameter]))
    }
    outcomes_to_keep <- outcomes_to_keep[which(as.character(outcomes_to_keep[,parameter])  %in% as.character( get(paste(parameter,"_list",sep=""))) ), ]
  }
  
  print(outcomes_to_keep)
  
   ###among those, keep the one that has the best generalisation value
  outcomes_to_keep <- outcomes_to_keep[which(outcomes_to_keep[,paste(measure,"_test",sep="")]==max(outcomes_to_keep[,paste(measure,"_test",sep="")],na.rm=T)),]
  ###selected method: keep the ones with highest score_training - for priority measure
  outcomes_to_keep <- outcomes_to_keep[which(outcomes_to_keep[,paste(measure,"_training",sep="")]==max(outcomes_to_keep[,paste(measure,"_training",sep="")],na.rm=T)),]
  print(outcomes_to_keep)
  
  
  ###selected method: keep the ones with highest score_training - for second emasure, in case there are ex-aequos
  outcomes_to_keep <- outcomes_to_keep[which(outcomes_to_keep[,paste(second_measure,"_training",sep="")]==max(outcomes_to_keep[,paste(second_measure,"_training",sep="")],na.rm=T)),]
  ###among those, keep the one that has the best generalisation value
  outcomes_to_keep <- outcomes_to_keep[which(outcomes_to_keep[,paste(second_measure,"_test",sep="")]==max(outcomes_to_keep[,paste(second_measure,"_test",sep="")],na.rm=T)),]
  print(outcomes_to_keep)
  
  
  
}

if (measure=="Fbeta"){
  measure <- paste(measure,beta,sep="_")
}
write.table(outcomes_to_keep,file=file.path(SAVEOUTPUT, "MachineLearning_outcomes_FINAL",paste("quality_scores_",measure,"_priority",priority,"_CHOSEN.txt",sep="")),col.names=T,row.names=F,quote=F,append=F)
