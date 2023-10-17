library(ggplot2)
library(gridExtra)
ANNOTATIONDIR <- "/media/cf19810/DISK3/Ants_behaviour_analysis/Data"

candidate_groomings <- read.table(paste(ANNOTATIONDIR,"/candidate_groomings_Vasudha2021_2022-11-23.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#num NAs
na_count <-sapply(candidate_groomings, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count$name<-rownames(na_count)

#select vars with nas
NA_grooms <- candidate_groomings[,colSums(is.na(candidate_groomings)) > 0]

### FIX THIS scatterplots and make them decent mean by lenght or smthing:

#https://stackoverflow.com/questions/35067600/plotting-time-series-values-using-bar-chart

par(mfrow=c(4,3))
for (x in names(NA_grooms)) {
   title("My 'Title' in a strange place", line = -1, outer = TRUE)
    barplot(candidate_groomings$duration_sec, NA_grooms[[x]],ylab="value",xlab =" ")
    title(cex.main=0.8,paste(x, "\n NAs",round(na_count[x,"na_count"]/nrow(candidate_groomings),1),"%"))
}


# for (x in names(NA_grooms)) {
#     ggplot(candidate_groomings, aes(x = duration_sec, y = NA_grooms[[x]])) +
#     geom_point() +
#     geom_rug()
# }


#-----------------

#the current annotation file FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv underwent cross-validation by Adriano
annotations_all <- read.csv(paste(ANNOTATIONDIR,"/Grooming_Classifier_CrossVal_ANNOTATIONS.csv",sep = ""), sep = ",")

###define new time columns called T_start_UNIX and T_stop_UNIX and delete T_start and T_stop
annotations_all$T_start_UNIX  <- as.POSIXct(annotations_all$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations_all$T_stop_UNIX   <- as.POSIXct(annotations_all$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )

annotations_all$duration      <- as.numeric(annotations_all$T_stop_UNIX - annotations_all$T_start_UNIX)

#compare original annotations with candidate grooming events
#look only at events under 50 sec
# par(mfrow=c(1,2))
# hist(annotations_all$duration[annotations_all$duration<50],breaks = 10,
#      ylim=c(0,1200),
#      main="annotations")
# 
# # hist(candidate_groomings$duration_sec[candidate_groomings$duration_sec<50],
# #      ylim=c(0,1200),
# #      main="candidate_groomings")

p1 <- ggplot(annotations_all[annotations_all$duration<50,], 
            aes(duration)) +
            geom_histogram(binwidth = 10) +
            scale_y_continuous(limits = c(0,1200))+
            theme_classic()+
            ggtitle("annotations_all")


p2 <- ggplot(candidate_groomings[candidate_groomings$duration_sec<50,], 
            aes(duration_sec, fill = as.factor(predicted_Hit))) +
            geom_histogram(binwidth = 10) +
            scale_y_continuous(limits = c(0,1200))+
            theme_classic() +
            theme(legend.position = c(0.8, 0.8))+
            ggtitle("candidate_groomings")

            
grid.arrange(p1, p2, nrow = 1)


#------------------------------------------------------------------------
library(data.table)

#SCORES
quality_scores <- read.table(paste(ANNOTATIONDIR,"/quality_scores_ClassifierValidation_Vasudha2021_Adriano2022.txt",sep=""),header=T,stringsAsFactors = F, sep=" ")

#substr(quality_scores$replicate,1,nchar(quality_scores$replicate)-2)
quality_scores$treatment <- str_sub(quality_scores$replicate,start=-2)

#specify reference values
quality_scores_limits <- data.frame(value=c(0.77,0.71),score=c("precision","sensitivity"))

ggplot(quality_scores[,c(1,4:6)]) +
  geom_segment(aes(x = 1, xend = 2, y = precision, yend = sensitivity)) + #, colour = treatment
  theme_classic() +
  geom_point(aes(x = 1,  y = precision)) +
  geom_point(aes(x = 2, y = sensitivity)) +
  scale_x_discrete(
    breaks = c("1", "2"),
    labels = c("precision", "sensitivity"),
    limits = c(1, 2)) +
  facet_wrap(~treatment) +
  #theme(legend.position="none") +
  labs(y = "value", x = "score") +
 #geom_hline(yintercept=0.77, linetype="dashed", color = "red")
  geom_hline(aes(yintercept = value, color = score), data = quality_scores_limits) 




# geom_point(data=chosen,
#            aes(x=precision_test, y=sensitivity_test), 
#            color='green3',
#            size=3) +
