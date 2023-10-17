##########################################################################################
############## BEH CLASSIFICATION ANALYSIS #################################################
##########################################################################################

#script dependant on BEH_MAIN_behaviours_analysis_fort081.R

#For previous versions of this script, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral




# FORECAST VERIFICATION
# https://www.swpc.noaa.gov/sites/default/files/images/u30/Forecast%20Verification%20Glossary.pdf
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity

fit_classifiers <- function (summary_AUTO){
  
  #remove ant names/rep/int/etc in pca (keep only vars)
  numeric_variable_list  <- names(summary_AUTO) [ which (!names(summary_AUTO) %in% c( "REPLICATE", "PERIOD","BEH","ROW","Act_Name","Rec_Name","ant1","ant2","INT","ACT","REC","pair","frame_start","frame_stop","disagreement","Hit","unique_interaction_id") )]
  metadata_variable_list <- names(summary_AUTO) [ which (names(summary_AUTO) %in% c( "REPLICATE", "PERIOD","BEH","ROW","Act_Name","Rec_Name","ant1","ant2","INT","ACT","REC","pair","frame_start","frame_stop","disagreement","Hit","unique_interaction_id") )]
  
  #######################################################
  ######### MISSING CASES ###############################
  #######################################################
  ###remove variables with too many NAs and keep only complete cases
  missing_hits   <- round((sapply(summary_AUTO[which(summary_AUTO$Hit==1),], function(x) sum(is.na(x))/length(x))),3)*100
  missing_misses <- round((sapply(summary_AUTO[which(summary_AUTO$Hit==0),], function(x) sum(is.na(x))/length(x))),3)*100
  to_remove      <- names(summary_AUTO)[which(missing_hits>5|missing_misses>15)] ##remove all variables that have more than 5% missing cases in minority class and more than 15% missing cases in majority case
  summary_AUTO.1 <- summary_AUTO[, which(!names(summary_AUTO)%in%to_remove)]
  
  # Many NAs are present only in the False positives (majority class), therefore dropping them is going to affect only very few cases of the True Postives
  ###NATH_FLAG: this removes 1 true possitives out of 37
  summary_AUTO_NAOmit <- summary_AUTO.1[complete.cases(summary_AUTO.1), ]
  
  ###############################################################################
  ###### NORMALISE VARS BEFORE CLASSIFICATION ###################################
  ###############################################################################
  #empty base
  summary_AUTO_NAOmit_transf <- data.frame()[1:nrow(summary_AUTO_NAOmit), ]
  #initialise list to contain the BestNormalize Objects to be reused for test data and to roll out to all colonies
  BN_list <- list()
  
  #transform variables if they are not normal
  for (variable in names(summary_AUTO_NAOmit)){
    if (variable %in%numeric_variable_list) {
      #find the best transformation
      #This function currently estimates the Yeo-Johnson, the Box Cox  (if the data is positive), the log_10(x+a), the square-root (x+a), the arcsinh and the ordered quantile normalization
      BNobject <- bestNormalize(summary_AUTO_NAOmit[,variable],standardize = FALSE,allow_orderNorm=F) #with standardize = FALSE also no_tranformation is attempted 
      
      ##fill in new object containing the transformed data
      ### DO NOT SCALE as would be impossible to transfer to different values 
      ###Though (scale and center variable reccomended for data prep for LDA in Metabolites 2014, 4, 433-452; doi:10.3390/metabo4020433)
      summary_AUTO_NAOmit_transf$var <- as.numeric(BNobject$x.t); 
      names(summary_AUTO_NAOmit_transf)[names(summary_AUTO_NAOmit_transf) == 'var'] <- paste(variable,class(BNobject$chosen_transform)[1],sep=".")
      
      ##save BN object to be reused later
      BNobject <- list(BNobject);names(BNobject) <- paste(variable,class(BNobject[[1]]$chosen_transform)[1],sep=".")
      BN_list <- c(BN_list,BNobject)
      
    }else{
      # print(paste("non numeric attribute. Pasting [",variable,"] in the new dataset",sep = " "))
      summary_AUTO_NAOmit_transf[variable]<- summary_AUTO_NAOmit[,variable]}
  }
  
  # transform Hit into a factor
  summary_AUTO_NAOmit_transf$Hit <- as.factor(summary_AUTO_NAOmit_transf$Hit)
  
  # create smaller dataset containing only variables of interest and Hit
  summary_LDA_vars_hit <- cbind(summary_AUTO_NAOmit_transf[,which(!names(summary_AUTO_NAOmit_transf)%in%metadata_variable_list)],Hit=summary_AUTO_NAOmit_transf$Hit) 
  summary_LDA_vars     <- summary_LDA_vars_hit[which(names(summary_LDA_vars_hit)!="Hit")]
  
  ##### Apply RELIEF ######
  estReliefF <- attrEval("Hit", summary_LDA_vars_hit, estimator="ReliefFequalK", # with ReliefFexpRank the output is very similar
                         ReliefIterations=30)
  
  estReliefF <- as.data.frame(estReliefF)
  estReliefF$Variable <- rownames(estReliefF); rownames(estReliefF) <- NULL
  #missing cases by var
  MissByVar <- round((sapply(summary_LDA_vars, function(x) sum(is.na(x)))/nrow(summary_LDA_vars)) * 100, 2)
  MissByVar <- as.data.frame(MissByVar)
  MissByVar$Variable <- rownames(MissByVar); rownames(MissByVar) <- NULL
  #add information of missing cases to the Relief feature score
  estReliefF <- plyr::join(x=estReliefF, y=MissByVar, type = "full", match = "all",by="Variable")
  #sort by score
  estReliefF <- estReliefF[order(estReliefF$estReliefF,decreasing = TRUE),]
  
  ##add statistical significance using Multisurf 
  neighbor.idx.observed <- find.neighbors(summary_LDA_vars, summary_LDA_vars_hit$Hit, k = 0, method = "multisurf")
  results.list <- stir(summary_LDA_vars, neighbor.idx.observed, k = k, metric = "manhattan", method = "multisurf")
  t_sorted_multisurf <- results.list$STIR_T
  t_sorted_multisurf$Variable <- rownames(t_sorted_multisurf)
  
  #add information of missing cases to the Relief feature score
  estReliefStir <- plyr::join(x=estReliefF, y=t_sorted_multisurf, type = "full", match = "all",by="Variable")
  
  #select only variables with p<0.01
  RELIEF_selected <- t_sorted_multisurf[which(t_sorted_multisurf$t.pval.adj<0.001),"Variable"]
  summary_LDA_vars_RELIEF <- summary_LDA_vars[, match(c(RELIEF_selected), names(summary_LDA_vars))] 
  BN_list  <- BN_list[match(RELIEF_selected,names(BN_list))]
  
  #check vars correlation
  results.cor_RELIEF <- cor(summary_LDA_vars_RELIEF) 
  #corrplot(results.cor_RELIEF, type="upper", tl.cex= .4,title = "\nresults.cor",)
  
  #add the Hit classes
  summary_LDA_vars_RELIEF_hit <- cbind(summary_LDA_vars_RELIEF,Hit=summary_AUTO_NAOmit_transf$Hit)
  
  ##### corelearn can also be used to try other models using CoreModel
  
  
  
  ###########################
  ###### DATA SAMPLING ######
  ###########################
  
  #TESTING SMOTE AS IT PROVED TO BE THE HISHEST SCORING FOR F1-SCORE FOR A DATASET (ABALONE) OF VERY SIMILAR CARACTHERISTICS (IMBALANCE, N CASES, N PARAMS)
  # IDEA DERIVED FROM TAHIR, Muhammad Atif; KITTLER, Josef; YAN, Fei. Inverse random under sampling for class imbalance problem and its application to multi-label classification. Pattern Recognition, 2012, 45.10: 3738-3750.
  
  #Now we test several sampling algorithms to balance the dataset.
  suppressMessages( train_sbc   <- SBC(summary_LDA_vars_RELIEF_hit,   "Hit") )#Under-Sampling Based on Clustering (SBC)
  suppressMessages( train_ros   <- ROS(summary_LDA_vars_RELIEF_hit,   "Hit") )#random over-sampling algorithm (ROS)
  suppressMessages( train_rus   <- RUS(summary_LDA_vars_RELIEF_hit,   "Hit") )#random under-sampling algorithm (RUS)
  suppressMessages( train_smote <- RSBID::SMOTE(summary_LDA_vars_RELIEF_hit, "Hit") )#synthetic minority over-sampling technique (SMOTE)
  
  try(suppressMessages( train_borderline_smote_type1 <-  BLSMOTE(X=summary_LDA_vars_RELIEF_hit[which(names(summary_LDA_vars_RELIEF_hit)!="Hit")],
                                                                 target=summary_LDA_vars_RELIEF_hit$Hit,
                                                                 K=10,C=10,dupSize=0,method =c("type1"))$data ),silent=T)
  try(suppressMessages( train_borderline_smote_type2 <-  BLSMOTE(X=summary_LDA_vars_RELIEF_hit[which(names(summary_LDA_vars_RELIEF_hit)!="Hit")],
                                                                 target=summary_LDA_vars_RELIEF_hit$Hit,
                                                                 K=10,C=10,dupSize=0,method =c("type2"))$data ),silent=T)
  if(exists("train_borderline_smote_type1")){
    names(train_borderline_smote_type1)[which(names(train_borderline_smote_type1)=="class")] <- "Hit";train_borderline_smote_type1$Hit <- as.factor(train_borderline_smote_type1$Hit); train_borderline_smote_type1 <- train_borderline_smote_type1[complete.cases(train_borderline_smote_type1),]
  }
  if(exists("train_borderline_smote_type2")){
    names(train_borderline_smote_type2)[which(names(train_borderline_smote_type2)=="class")] <- "Hit";train_borderline_smote_type2$Hit <- as.factor(train_borderline_smote_type2$Hit); train_borderline_smote_type2 <- train_borderline_smote_type2[complete.cases(train_borderline_smote_type2),]
  }
  
  ##############################
  ### LDA / QDA ################
  ##############################
  
  
  
  # ### Quadratic Discriminant Analysis
  # #train on the whole and class-balanced datasets
  # fit_QDA_all <- qda(Hit ~ ., data = summary_LDA_vars_RELIEF_hit)
  # fit_QDA_sbc <- qda(Hit ~ ., data = train_sbc)
  # fit_QDA_ros <- qda(Hit ~ ., data=train_ros)
  # fit_QDA_rus <- qda(Hit ~ ., data=train_rus)
  # fit_QDA_smote <- qda(Hit ~ ., data=train_smote)
  # fit_QDA_borderline_smote_type1 <- qda(Hit ~ ., data=train_borderline_smote_type1)
  # fit_QDA_borderline_smote_type2 <- qda(Hit ~ ., data=train_borderline_smote_type2)
  
  ##############################
  ### RANDOM FOREST ############
  ##############################
  
  #Random forests is an ensemble learning method for classification that operates by constructing a multitude of decision
  #trees at training time. For classification tasks, the output of the random forest is the class selected by most trees.
  # to plot RF: https://rpubs.com/markloessi/498787
  
  #### hyperparam tuning: Random Search
  # you can use the trainControl function to specify a number of parameters (including sampling parameters) in your model. 
  #The object that is outputted from trainControl will be provided as an argument for RF train
  control <- trainControl(method="repeatedcv", number=10, repeats=10, search="random") # Repeated K Fold Cross-Validation
  
  #train on the class-balanced dataset
  #train on the whole and class-balanced datasets
  fit_RF_all                    <- randomForest::randomForest(Hit ~ ., data = summary_LDA_vars_RELIEF_hit   , tuneLength=20, trControl=control)
  fit_RF_sbc                    <- randomForest::randomForest(Hit ~ ., data = train_sbc                     , tuneLength=20, trControl=control)
  fit_RF_ros                    <- randomForest::randomForest(Hit ~ ., data = train_ros                     , tuneLength=20, trControl=control)
  fit_RF_rus                    <- randomForest::randomForest(Hit ~ ., data = train_rus                     , tuneLength=20, trControl=control)
  fit_RF_smote                  <- randomForest::randomForest(Hit ~ ., data = train_smote                   , tuneLength=20, trControl=control)
  if(exists("train_borderline_smote_type1")){
    fit_RF_borderline_smote_type1 <- randomForest::randomForest(Hit ~ ., data = train_borderline_smote_type1  , tuneLength=20, trControl=control)
  }else{
    fit_RF_borderline_smote_type1 <- NULL
  }
  if(exists("train_borderline_smote_type2")){
    
    fit_RF_borderline_smote_type2 <- randomForest::randomForest(Hit ~ ., data = train_borderline_smote_type2  , tuneLength=20, trControl=control)
  }else{
    fit_RF_borderline_smote_type2 <- NULL
  }
    #tuneLength : It allows system to tune algorithm automatically. It indicates the number of different values to try for each tunning parameter.
  #trControl : looks for the random search parameters established by caret::trainControl
  #plot(fit_RF_sbc)
  
  
  ##############################
  ### SUPPORT VECTOR MACHINE####
  ##############################
  fit_SVM_all                    <- svm(Hit ~ ., data = summary_LDA_vars_RELIEF_hit  , scale = FALSE, kernel = "radial", cost = 5)
  fit_SVM_sbc                    <- svm(Hit ~ ., data = train_sbc                    , scale = FALSE, kernel = "radial", cost = 5)
  fit_SVM_ros                    <- svm(Hit ~ ., data = train_ros                    , scale = FALSE, kernel = "radial", cost = 5)
  fit_SVM_rus                    <- svm(Hit ~ ., data = train_rus                    , scale = FALSE, kernel = "radial", cost = 5)
  fit_SVM_smote                  <- svm(Hit ~ ., data = train_smote                  , scale = FALSE, kernel = "radial", cost = 5)
  if(exists("train_borderline_smote_type1")){
    fit_SVM_borderline_smote_type1 <- svm(Hit ~ ., data = train_borderline_smote_type1 , scale = FALSE, kernel = "radial", cost = 5)
  }else{
    fit_SVM_borderline_smote_type1 <- NULL
    }
  if(exists("train_borderline_smote_type2")){
    fit_SVM_borderline_smote_type2 <- svm(Hit ~ ., data = train_borderline_smote_type2 , scale = FALSE, kernel = "radial", cost = 5)
  }else{
    fit_SVM_borderline_smote_type2 <- NULL
  }
  
  
  
  # ############ RULE EXTRACTION FOR RF ###################
  # # Stable and Interpretable RUle Set for RandomForests,
  # # Should output an easily interpretable ruleset for classification of a RandomForest algorithm
  # # compared to the RandomForest output, the SIRUS rule selection is more stable with respect to data perturbation.
  # # (Benard et al. 2021)
  # summary_LDA_vars_RELIEF_hit$Hit  <- as.numeric(as.character(summary_LDA_vars_RELIEF_hit$Hit))
  # train_sbc$Hit                    <- as.numeric(as.character(train_sbc$Hit))
  # train_ros$Hit                    <- as.numeric(as.character(train_ros$Hit))
  # train_rus$Hit                    <- as.numeric(as.character(train_rus$Hit))
  # train_smote$Hit                  <- as.numeric(as.character(train_smote$Hit))
  # train_borderline_smote_type1$Hit                  <- as.numeric(as.character(train_borderline_smote_type1$Hit))
  # train_borderline_smote_type2$Hit                  <- as.numeric(as.character(train_borderline_smote_type2$Hit))
  # 
  # # train SIRUS on the whole and class-balanced dataset
  # 
  # TUNE_SIRUS <- FALSE
  # if (TUNE_SIRUS) {
  #   # Estimate the optimal hyperparameter p0 used to select rules in sirus.fit using cross-validation (Benard et al. 2020, 2021).
  #   # p0 is the only tunable hyperparameter.
  #   cross_val_p0_all <- sirus.cv(summary_LDA_vars_RELIEF_hit[which(names(summary_LDA_vars_RELIEF_hit)!="Hit")], summary_LDA_vars_RELIEF_hit$Hit, type = "classif") #takes a while  
  #   cross_val_p0_sbc <- sirus.cv(train_sbc[which(names(train_sbc)!="Hit")], train_sbc$Hit, type = "classif") #takes a while  
  #   cross_val_p0_ros <- sirus.cv(train_ros[which(names(train_ros)!="Hit")], train_ros$Hit, type = "classif") #takes a while  
  #   cross_val_p0_rus <- sirus.cv(train_rus[which(names(train_rus)!="Hit")], train_rus$Hit, type = "classif") #takes a while  
  #   cross_val_p0_smote <- sirus.cv(train_smote[which(names(train_smote)!="Hit")], train_smote$Hit, type = "classif") #takes a while  
  #   cross_val_p0_borderline_smote_type1 <- sirus.cv(train_borderline_smote_type1[which(names(train_borderline_smote_type1)!="Hit")], train_borderline_smote_type1$Hit, type = "classif") #takes a while  
  #   cross_val_p0_borderline_smote_type2 <- sirus.cv(train_borderline_smote_type2[which(names(train_borderline_smote_type2)!="Hit")], train_borderline_smote_type2$Hit, type = "classif") #takes a while  
  #   
  #   fit_RFsirus_all   <- sirus.fit(summary_LDA_vars_RELIEF_hit[which(names(summary_LDA_vars_RELIEF_hit)!="Hit")], summary_LDA_vars_RELIEF_hit$Hit, type = "classif",p0= cross_val_p0_all$p0.stab,verbose=F)
  #   fit_RFsirus_sbc   <- sirus.fit(train_sbc[which(names(train_sbc)!="Hit")], train_sbc$Hit, type = "classif",p0= cross_val_p0_sbc$p0.stab,verbose=F)
  #   fit_RFsirus_ros   <- sirus.fit(train_ros[which(names(train_ros)!="Hit")], train_ros$Hit, type = "classif",p0= cross_val_p0_ros$p0.stab,verbose=F) 
  #   fit_RFsirus_rus   <- sirus.fit(train_rus[which(names(train_rus)!="Hit")], train_rus$Hit, type = "classif",p0= cross_val_p0_rus$p0.stab,verbose=F) 
  #   fit_RFsirus_smote <- sirus.fit(train_smote[which(names(train_smote)!="Hit")], train_smote$Hit, type = "classif",p0= cross_val_p0_smote$p0.stab,verbose=F)
  #   fit_RFsirus_borderline_smote_type1 <- sirus.fit(train_borderline_smote_type1[which(names(train_borderline_smote_type1)!="Hit")], train_borderline_smote_type1$Hit, type = "classif",p0= cross_val_p0_borderline_smote_type1$p0.stab,verbose=F)
  #   fit_RFsirus_borderline_smote_type2 <- sirus.fit(train_borderline_smote_type2[which(names(train_borderline_smote_type2)!="Hit")], train_borderline_smote_type2$Hit, type = "classif",p0= cross_val_p0_borderline_smote_type2$p0.stab,verbose=F)
  #   
  # }else{
  #   fit_RFsirus_all   <- sirus.fit(summary_LDA_vars_RELIEF_hit[which(names(summary_LDA_vars_RELIEF_hit)!="Hit")], summary_LDA_vars_RELIEF_hit$Hit, type = "classif",verbose=F)
  #   fit_RFsirus_sbc   <- sirus.fit(train_sbc[which(names(train_sbc)!="Hit")], train_sbc$Hit, type = "classif",verbose=F)
  #   fit_RFsirus_ros   <- sirus.fit(train_ros[which(names(train_ros)!="Hit")], train_ros$Hit, type = "classif",verbose=F) 
  #   fit_RFsirus_rus   <- sirus.fit(train_rus[which(names(train_rus)!="Hit")], train_rus$Hit, type = "classif",verbose=F) 
  #   fit_RFsirus_smote <- sirus.fit(train_smote[which(names(train_smote)!="Hit")], train_smote$Hit, type = "classif",verbose=F)
  #   fit_RFsirus_borderline_smote_type1 <- sirus.fit(train_borderline_smote_type1[which(names(train_borderline_smote_type1)!="Hit")], train_borderline_smote_type1$Hit, type = "classif",verbose=F)
  #   fit_RFsirus_borderline_smote_type2 <- sirus.fit(train_borderline_smote_type2[which(names(train_borderline_smote_type2)!="Hit")], train_borderline_smote_type2$Hit, type = "classif",verbose=F)
  #   
  # }
  # 
  
  
  #print the rule set (saved at the bottom of MAIN)
  # SirusRules <- list(SIRUS_rules_all =sirus.print(fit_RFsirus_all, digits = 5),
  #                    SIRUS_rules_sbc =sirus.print(fit_RFsirus_sbc, digits = 5),
  #                    SIRUS_rules_ros =sirus.print(fit_RFsirus_ros, digits = 5),
  #                    SIRUS_rules_rus =sirus.print(fit_RFsirus_rus, digits = 5),
  #                    SIRUS_rules_smote =sirus.print(fit_RFsirus_smote, digits = 5),
  #                    SIRUS_rules_borderline_smote_type1 =sirus.print(fit_RFsirus_borderline_smote_type1, digits = 5),
  #                    SIRUS_rules_borderline_smote_type2 =sirus.print(fit_RFsirus_borderline_smote_type2, digits = 5)
  #                    )
  
  ###Support Vector Machine
  
  
  return (
    list(
      selected_variables_BNobject_list     = BN_list
      ,fits                                = list(
        # fit_QDA_all   = fit_QDA_all
        # ,fit_QDA_sbc   = fit_QDA_sbc
        # ,fit_QDA_ros   = fit_QDA_ros
        # ,fit_QDA_rus   = fit_QDA_rus
        # ,fit_QDA_smote = fit_QDA_smote
        # ,fit_QDA_borderline_smote_type1 = fit_QDA_borderline_smote_type1
        # ,fit_QDA_borderline_smote_type2 = fit_QDA_borderline_smote_type2
        
        fit_RF_all   = fit_RF_all
        ,fit_RF_sbc   = fit_RF_sbc
        ,fit_RF_ros   = fit_RF_ros
        ,fit_RF_rus   = fit_RF_rus
        ,fit_RF_smote = fit_RF_smote
        ,fit_RF_borderline_smote_type1 = fit_RF_borderline_smote_type1
        ,fit_RF_borderline_smote_type2 = fit_RF_borderline_smote_type2
        
        ,fit_SVM_all   = fit_SVM_all
        ,fit_SVM_sbc   = fit_SVM_sbc
        ,fit_SVM_ros   = fit_SVM_ros
        ,fit_SVM_rus   = fit_SVM_rus
        ,fit_SVM_smote = fit_SVM_smote
        ,fit_SVM_borderline_smote_type1 = fit_SVM_borderline_smote_type1
        ,fit_SVM_borderline_smote_type2 = fit_SVM_borderline_smote_type2
        
        # ,fit_RFsirus_all = fit_RFsirus_all
        # ,fit_RFsirus_sbc = fit_RFsirus_sbc
        # ,fit_RFsirus_ros = fit_RFsirus_ros
        # ,fit_RFsirus_rus = fit_RFsirus_rus
        # ,fit_RFsirus_smote = fit_RFsirus_smote
        # ,fit_RFsirus_borderline_smote_type1 = fit_RFsirus_borderline_smote_type1
        # ,fit_RFsirus_borderline_smote_type2 = fit_RFsirus_borderline_smote_type2
      )
      # ,SirusRules = SirusRules
      ,SirusRules = NULL
    )
  )
}

predict_class <- function (summary_AUTO,BN_list,classifier){
  ###I. to predict values from summary_AUTO, we need to:
  ###1/ select the appropriate variables (those should be evident from the names of BN_list)
  ###2/ perform the appropriate transformations (using BN_list)
  
  summary_LDA_vars <- data.frame()[1:nrow(summary_AUTO), ]
  for ( selected_variable in names(BN_list)){
    summary_LDA_vars[selected_variable] <-                   predict(BN_list[[selected_variable]], newdata = summary_AUTO[,unlist(strsplit(selected_variable,split="\\."))[1]])
  }  
  ###in very rare cases there may be a +Inf or -Inf produced in transformation. Replace those with NA
  summary_LDA_vars[summary_LDA_vars==-Inf] <- NA
  summary_LDA_vars[summary_LDA_vars==Inf] <- NA
  
  if (grepl("QDA",names(classifier))){
    predicted_class                                   <- as.numeric(as.character(predict(classifier[[names(classifier)]], summary_LDA_vars, type="response")$class ))
  }
  if (grepl("_RF_",names(classifier))){
    predicted_class                                   <- as.numeric(as.character(predict(classifier[[names(classifier)]], summary_LDA_vars, type="response")       ))
  }
  if (grepl("_RFsirus_",names(classifier))){
    predicted_class                                   <- rep(NA, nrow(summary_LDA_vars))
    predicted_class[complete.cases(summary_LDA_vars)] <-  as.numeric(sirus.predict(classifier[[names(classifier)]], summary_LDA_vars[complete.cases(summary_LDA_vars),])>0.5) ###assign class '1' whenever probability > 0.5
  }
  if (grepl("_SVM_",names(classifier))){
    predicted_class                                   <- rep(NA, nrow(summary_LDA_vars))
    predicted_class[complete.cases(summary_LDA_vars)] <-  as.numeric(as.character(predict(classifier[[names(classifier)]], summary_LDA_vars[complete.cases(summary_LDA_vars),])))
  }
  
  return(predicted_class)
}

quality_scores <- function(true_false_positive_negatives,beta){
  TP <- as.numeric(true_false_positive_negatives["true_positives"])
  FP <- as.numeric(true_false_positive_negatives["false_positives"])
  FN <- as.numeric(true_false_positive_negatives["false_negatives"])
  
  #################################################
  ### Algorithm quality scores ####################
  #################################################
  
  ### critical success index (CSI) https://www.swpc.noaa.gov/sites/default/files/images/u30/Forecast%20Verification%20Glossary.pdf
  # # CSI is a verification measure of categorical forecast performance
  # # calculated as CSI= (hits)/(hits + false alarms + misses)
  # # The CSI is not affected by the number of non-event forecasts that verify (correct rejections).
  CSI <- TP / (TP + FP + FN)
  
  ### Fbeta score https://en.wikipedia.org/wiki/F-score
  # F-score or F-measure is a measure of a test's accuracy. It is calculated from the precision and recall of the test, where the precision is the number of true positive results divided by the number of all positive results, including those not identified correctly, and the recall is the number of true positive results divided by the number of all samples that should have been identified as positive. Precision is also known as positive predictive value, and recall is also known as sensitivity in diagnostic binary classification. 
  Fbeta <- ( (1 + beta^2)*TP   )/(  (1+beta^2)*TP + beta^2 * FN + FP       )
  
  ### precision = positive predictive value https://en.wikipedia.org/wiki/Likelihood_ratios_in_diagnostic_testing#positive_likelihood_ratio
  # = proportion of detected groomings that are true groomings
  precision <- TP / (TP + FP)
  
  ### True positive rate (TPR) =recall = sensitivity
  # = proportion of true groomings that are detected
  sensitivity <- TP / (TP + FN)
  
  return(c(CSI=CSI,Fbeta=Fbeta,precision=precision,sensitivity=sensitivity))
}
