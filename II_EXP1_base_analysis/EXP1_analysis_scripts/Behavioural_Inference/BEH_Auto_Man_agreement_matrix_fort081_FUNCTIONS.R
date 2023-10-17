make_binary_interaction_matrix <- function(full_ant_list,focal_list,IF_frames,summary_object){
  # allpairs <- data.frame(combinations(length(ant_list), 2, v = sort(as.numeric(ant_list)) , repeats.allowed = FALSE))
    allpairs <- expand.grid(full_ant_list,focal_list)###list all pairs
    
    allpairs <- allpairs[which(allpairs[,1]!=allpairs[,2]),]##remove self pairs
    
    allpairs$pair <- apply(allpairs[,],1,function(x){paste(sort(x),collapse = "_") }) ###sort pair in increasing order
    allpairs <- allpairs[which(!duplicated(allpairs$pair)),]##remove duplicates
  
    ids_pairs <- subset(allpairs,select = "pair")
  
  # # initialize interaction matrix each rows represent a binary array, one for each ids pairs, with 1s on the interactions and 0s elsewhere
  int_mat <- matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames)[1] ))
  rownames(int_mat) <- ids_pairs$pair
  colnames(int_mat) <- c(IF_frames$frame_num)
  if (nrow(summary_object)>0){
    for (i in 1:nrow(summary_object))
    {
      PAIR <- summary_object$pair[i]
      int_mat[PAIR,] <- int_mat[PAIR,] + c(rep(0,( summary_object$frame_start[i]-1)),
                                           rep(1,(summary_object$frame_stop[i]) - summary_object$frame_start[i] + 1),
                                           rep(0,(length(IF_frames$frame_num) - summary_object$frame_stop[i])))
    }
  }
  return(int_mat)
}

auto_manual_agreement <- function(summary_AUTO , summary_MANUAL, list_IF_Frames, list_replicate_full_ant_list, list_replicate_focal_list ){
  
  ##########################################################################################
  ############## AUTO-MAN AGREEMENT MATRIX #################################################
  ##########################################################################################
  
  #script dependant on BEH_MAIN_behaviours_analysis_fort081.R
  
  #For previous versions of this script, explore: 
  # https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral
  # print(paste("AUTO-MAN AGREEMENT MATRIX"))
  
  ###initialise number of TN, TP, FN, FP
  true_negatives   <- 0
  true_positives   <- 0
  false_negatives  <- 0
  false_positives  <- 0
  
  ###initialise disagreement in summary_AUTO
  summary_AUTO$disagreement <- NA
  
  for (REPLICATE in unique(summary_AUTO$REPLICATE)){
    for (PERIOD in unique(summary_AUTO[which(summary_AUTO$REPLICATE==REPLICATE),"PERIOD"])){
      summary_AUTO_REP_PER <- summary_AUTO    [ which(summary_AUTO   $REPLICATE==REPLICATE & summary_AUTO   $PERIOD==PERIOD),]   
      summary_MAN_REP_PER  <- summary_MANUAL  [ which(summary_MANUAL $REPLICATE==REPLICATE & summary_MANUAL $PERIOD==PERIOD),]
      IF_frames            <- list_IF_Frames[[paste(c("IF_frames",REPLICATE,PERIOD),collapse="_")]]
      full_ant_list        <- list_replicate_full_ant_list[[paste(c("replicate_list",REPLICATE),collapse="_")]]
      focal_list           <- list_replicate_focal_list [[paste(c("replicate_list",REPLICATE),collapse="_")]]
      if (is.null(focal_list)){
        focal_list <- full_ant_list
      }
      
      # ###NATH_FLAG:INSTEAD USE COMBINATIONS, and extract your ant list from the objects you have
      # ant_list <- sort(na.omit(unique(c(summary_AUTO_REP_PER$ant1,summary_AUTO_REP_PER$ant2,summary_MAN_REP_PER$ant1,summary_MAN_REP_PER$ant2))))
      # if(focal){
      #   ant_list <- sort(na.omit(unique(c(ant_list,gsub(paste(REPLICATE,"_",sep=""),"",focal_list[which(grepl(REPLICATE,focal_list))])))))
      # }
      
      
      ###create interaction matrices
      int_mat_manual <- make_binary_interaction_matrix (full_ant_list,focal_list,IF_frames,summary_MAN_REP_PER )
      int_mat_auto   <- make_binary_interaction_matrix (full_ant_list,focal_list,IF_frames,summary_AUTO_REP_PER )
      
      #####NATH_FLAG: use these to calculate the number of true positives, number of true negatives, number of false positives and number of false negatives
      true_negatives   <- true_negatives  + sum(as.numeric(int_mat_manual==0) * as.numeric(int_mat_auto ==0))
      true_positives   <- true_positives  + sum(as.numeric(int_mat_manual==1) * as.numeric(int_mat_auto ==1))
      false_negatives  <- false_negatives + sum(as.numeric(int_mat_manual==1) * as.numeric(int_mat_auto ==0))
      false_positives  <- false_positives + sum(as.numeric(int_mat_manual==0) * as.numeric(int_mat_auto ==1))
      
      #####then for each automatic interaction, determine if it can be considered as a true positive (hit) or a false positive (miss)
      #Calculate the % disagreement per each interaction
      int_mat_err <- int_mat_manual - int_mat_auto
      for (i in which(summary_AUTO   $REPLICATE==REPLICATE & summary_AUTO   $PERIOD==PERIOD)) ###loop over the indices in summary auto object that correspond to replicate and period of interest
      {
        PAIR <- summary_AUTO$pair[i]
        Col_Indices <- summary_AUTO$frame_start[i] :  summary_AUTO$frame_stop[i]
        Row_Indices <- which(PAIR == rownames(int_mat_err))
        Overlap <- int_mat_err[ Row_Indices, Col_Indices]
        summary_AUTO$disagreement[i] <-   sum(Overlap)/length(Overlap)
      } 
      ## visualise the agrement:
      ## An interaction with mean =0 has total agreement (0)
      ## An interaction with mean <0 is a false positive (-1)
      ## An interaction with mean >0 is false negative (+1) ###NATH_FLAG cannot happen as here we are only looking at automatic interactions i.e.e 'positives'
      
      ## explore
      #plot(disagreement ~ duration_sec, summary_AUTO_REP_PER); abline(h=0, lty=2)
      
      
    }
  }
  ## APPLY Disagreement THRESHOLDS
  summary_AUTO$Hit <- NA
  summary_AUTO$Hit [which(abs(summary_AUTO$disagreement) <=  DISAGREEMENT_THRESH)] <- 1 ## 
  summary_AUTO$Hit [which(abs(summary_AUTO$disagreement) >   DISAGREEMENT_THRESH)] <- 0
  
  # print(paste("AUTO-MAN AGREEMENT MATRIX","PERFORMED"))
  
  true_false_positive_negatives <- data.frame(true_negatives,true_positives,false_negatives,false_positives)
  return(list(summary_AUTO=summary_AUTO,true_false_positive_negatives=true_false_positive_negatives))
}

