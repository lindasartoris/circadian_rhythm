library(dplyr)
library(caret)
library(MLeval)
library(randomForest)
library(datasets)
library(data.table)
library(ggplot2)
library(viridis)

#if (USER=="Supercomputer1") {
  WORKDIR <- "/home/cf19810/Documents"
  DATADIR <- paste(WORKDIR,"MachineLearning_outcomes_12-11-22",sep="/")
  #SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
  #EXPDATADIR <- "/media/bzniks/DISK4/ADRIANO/EXPERIMENT_DATA" #"/home/cf19810/Documents/Ants_behaviour_analysis/Data"
#}

###############################################################################
#### READ CHOSEN METHOD #######################################################
###############################################################################
chosen <- read.table(paste(DATADIR,"/quality_scores_CHOSEN.txt",sep=""),header=T,stringsAsFactors = F)
all_methods <- read.table(paste(DATADIR,"/quality_scores.txt",sep=""),header=T,stringsAsFactors = F)
  
###############################################################################
###### EXTRACT CHOSEN PARAMETERS FROM CHOSEN ##################################
subDir                      <- paste0("Loop_ID_",chosen[,"Loop_ID"]) 
  
##Load_classifier
rf_model        <- readRDS (file.path(DATADIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep="")))
rf_model_LIST       <- list(readRDS (file.path(DATADIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep=""))))
names(rf_model_LIST) <- chosen[,"classifier"]

#plot model using randomforest package tools
###Error plot
rf_model.legend <- if (is.null(rf_model$test$err.rate)) {colnames(rf_model$err.rate)} else {colnames(rf_model$test$err.rate)}
plot(rf_model, log="y", main = "log of error rate", )
legend("top", cex =1, legend=rf_model.legend, lty=c(1,2,3), col=c(1,2,3), horiz=T)

#####error plot in ggplot
# Get OOB data from plot and coerce to data.table
oobData <- as.data.table(plot(rf_model))
# Define trees as 1:ntree
oobData[, trees := .I]
# Cast to long format
oobData2 <- melt(oobData, id.vars = "trees")
setnames(oobData2, "value", "error")
# Plot using ggplot
ggplot(data = oobData2, aes(x = trees, y = error, color = variable)) + 
  theme_bw() +
  geom_line() +
  #scale_x_discrete(name = "Class", labels = c("OOB", "Non-grooming", "Grooming"))
  scale_color_manual(name="class",
                   labels=c("OOB", "Non-grooming", "Grooming"),
                   values=c("black","brown1","forestgreen"))

### Variables importance plot
#measure of the Mean decrease Gini
imp <- varImpPlot(rf_model) # let's save the varImp object

# this part just creates the data.frame for the plot part
imp <- as.data.frame(imp)
imp$varnames <- rownames(imp) # row names to column
rownames(imp) <- NULL
#assign category
imp$var_categ  <- ifelse(grepl('ACT', imp$varnames), 'ACT',
                       ifelse(grepl('REC', imp$varnames), 'REC', "BOTH"))


# this is the plot part, be sure to use reorder with the correct measure name
ggplot(imp, aes(x=reorder(varnames, MeanDecreaseGini), y=MeanDecreaseGini, color=as.factor(var_categ))) + 
  geom_point() +
  geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseGini)) +
  scale_color_discrete(name="Individual measured") +
  theme_bw() +
  theme(legend.position="top") +
  ylab("Mean Decrease Gini") +
  xlab("Variable Name") +
  coord_flip()

varUsed(rf_model, by.tree=FALSE, count=TRUE)


##### precision vs sensitivity
ggplot(all_methods$precision_test,all_methods$sensitivity_test)

  ggplot(all_methods, aes(x=precision_test, y=sensitivity_test)) + 
    geom_point(aes(colour = classifier)) +
    scale_colour_viridis_d(option = "plasma") + #inferno
    geom_point(data=chosen,
               aes(x=precision_test, y=sensitivity_test), 
               color='green3',
               size=3) +
    theme_bw() +
    #stat_ellipse(aes(color = classifier), geom="polygon",level=0.95,alpha=0) +
    geom_text(data=chosen,
              aes(x=precision_test, y=sensitivity_test,label=paste(chosen$precision_test,chosen$sensitivity_test,sep=", ")),
              hjust=-0.3, vjust=-0.3, color="green3") +
    labs(caption = paste0("N=",nrow(all_methods)))
  
  
  #uncouple classif from balancer to give shapes to balancer

#-------------------------------------

#Assuming you already have a vector of probabilities (called probs) computed with your model and the true class labels are in your data frame as df$label (0 and 1) this code should work:
# https://stats.stackexchange.com/questions/10501/calculating-aupr-in-r
fg <- probs[df$label == 1]
bg <- probs[df$label == 0]

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

#PS: The only disconcerting thing is you use scores.class0 = fg when fg is computed for label 1 and not 0.

#-------------------------------------



