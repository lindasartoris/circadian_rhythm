####################################################################################
#### THIS SCRIPT CONTAINS:
####
####################################################################################
library(FortMyrmidon)
library(ggplot2)
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)
library(viridis)
library(stringr)
library(dplyr)

##PLOTTING
STYLE_NOVIR <- list(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                    theme_bw(),
                    scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

# Define custom color palette
custom_palette <- c("treated nurse" = "#333333", "untreated forager" = "#FDE725", "untreated nurse" = "#1F9E89")

# 
# WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
# DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"
# 
# #SpaceUse <- read.table(paste(DATADIR,"/AntTasks_SpaceUse_july2022.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2023-02-27.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# 
# SpaceUse$REP_treat <- SpaceUse$colony
# SpaceUse$delta_time_inside <- NULL
# SpaceUse$delta_time_outside <-SpaceUse$prop_inside_24hPRE -  SpaceUse$prop_inside_24hPOST
# 
# 
# ## Colony sizes info
# SpaceUse$size_treat <- str_sub(SpaceUse$colony,-2,-1)
# Mean_ants_exp <- aggregate(colony_size ~ size_treat, FUN=mean, na.action=na.omit, SpaceUse)
# SD_ants_exp <- aggregate(colony_size ~ size_treat, FUN=sd, na.action=na.omit, SpaceUse); colnames(SD_ants_exp) [match("colony_size",colnames(SD_ants_exp))] <- "SD_received"
# Colony.size    <-  plyr::join(x=Mean_ants_exp, y=SD_ants_exp, type = "full", match = "all")
# data.frame(size_treat=Colony.size$size_treat, Colony_size=sprintf("%.2f \U00B1 %.2f",Colony.size$colony_size,Colony.size$SD_received))
# 
# 
# #select some SpaceUse cols
# SpaceUse <- SpaceUse[,c("REP_treat","antID","delta_time_outside")] #remove the ant_task as it is old comapred to the one in metadata
# 
# ### Get info from metadata
# SpaceMeta <- list(SpaceUse,metadata)
# SpaceMeta <- Reduce(function(x, y) merge(x, y, all=TRUE), SpaceMeta)
# #clean
# SpaceMeta <- SpaceMeta[which(!is.na(SpaceMeta$antID)),]
# 
# N_exp <- as.data.frame(table(SpaceMeta$Exposed))
# 
# #LEVELS RENAMING
# # Rename by name
# SpaceMeta$size_treat <- as.factor(SpaceMeta$size_treat)
# levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="BS"] <- "Big Sham"
# levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="BP"] <- "Big Pathogen"
# levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="SS"] <- "Small Sham"
# levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="SP"] <- "Small Pathogen"
# 
# #remove queens and dead ants
# SpaceMeta <- SpaceMeta[which(SpaceMeta$AntTask!= "queen"),]
# SpaceMeta <- SpaceMeta[which(SpaceMeta$IsAlive== TRUE),]
# 
# 
# # Rename levels
# N_exp$Var1 <- as.factor(N_exp$Var1)
# levels(N_exp$Var1)[levels(N_exp$Var1)==TRUE] <- "treated"
# levels(N_exp$Var1)[levels(N_exp$Var1)==FALSE] <- "untreated"
# 
# 
# #CHECK IF THERE ARE ANY EXPOSED NURSES LABELED AS FORAGERS
# SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),] #only 3 individuals, no issues in code!
# #change task for the sake of plotting
# SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),"AntTask"] <- "nurse"
# 
# #aggregg for boxplot - ants means, by colony
# #mean by colony, AntTask, exposed, size_treat
# Mean_SpaceMeta <- aggregate(delta_time_outside ~ REP_treat + AntTask + Exposed + size_treat, FUN=mean, na.rm=T, na.action=na.pass, SpaceMeta)
# SD_SpaceMeta <- aggregate(delta_time_outside ~ REP_treat + AntTask + Exposed + size_treat, FUN=sd, na.rm=T, na.action=na.pass, SpaceMeta); colnames(SD_SpaceMeta) [match("delta_time_outside",colnames(SD_SpaceMeta))] <- "SD_delta_time_outside"
# 
# #merge dfs
# Mean_SpaceMeta <- plyr::join (x = Mean_SpaceMeta , y=SD_SpaceMeta, type = "right", match = "all")
# # Rename levels
# Mean_SpaceMeta$Exposed <- as.factor(Mean_SpaceMeta$Exposed)
# levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==TRUE] <- "treated"
# levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==FALSE] <- "untreated"
# 
# Mean_SpaceMeta$Ant_status <- paste(Mean_SpaceMeta$Exposed, Mean_SpaceMeta$AntTask)
# 
# ### for barplot - colony means
# Mean_SpaceMeta_Rep <- aggregate(delta_time_outside ~ Ant_status + size_treat, FUN=mean, na.rm=T, na.action=na.pass, Mean_SpaceMeta)
# SE_SpaceMeta_Rep   <- aggregate(delta_time_outside ~ Ant_status + size_treat, FUN=std.error, na.rm=T, na.action=na.pass, Mean_SpaceMeta); colnames(SE_SpaceMeta_Rep) [match("delta_time_outside",colnames(SE_SpaceMeta_Rep))] <- "SE_delta_time_outside"
# #merge dfs
# Mean_SpaceMeta_Rep <- plyr::join (x = Mean_SpaceMeta_Rep , y=SE_SpaceMeta_Rep, type = "right", match = "all")
# 
# 
# 
# # #### BOXPLOT OF USED SPACE FOR DELTA TIME INSIDE
# # ggplot(Mean_SpaceMeta, aes(x=size_treat, y=delta_time_outside,color=Ant_status))+
# #   # geom_errorbar( aes(x=size_treat,ymin=delta_time_outside-SD_delta_time_outside, ymax=delta_time_outside+SD_delta_time_outside),position=position_dodge2(width=0.8, preserve = "single"))+
# #   # geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
# #  # geom_jitter(aes(fill = Ant_status))+
# #   geom_boxplot(position = position_dodge(width = 0.8, preserve = "single"))+
# #   #geom_point(size=0.8,position=position_dodge2(width = 0.8, preserve = "single"))+
# #   #facet_wrap(~Exposed) +
# #   STYLE_NOVIR +
# #   labs(title= "Space Use, ant means by colony",
# #        subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))
# 
# 
# 
# 
# 
# 
# # Define custom color palette
# custom_palette <- c("treated nurse" = "#333333", "untreated forager" = "#FDE725", "untreated nurse" = "#1F9E89")
# 
# # Create the boxplot
# ggplot(Mean_SpaceMeta, aes(x=size_treat, y=delta_time_outside, color=Ant_status)) +
#   geom_boxplot(position = position_dodge(width = 0.8, preserve = "single")) +
#   scale_color_manual(values = custom_palette) + # Apply custom color palette
#   STYLE_NOVIR +
#   labs(x="")
#   #labs(title= "Space Use, ant means by colony",
#   #     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))
# 
# 
# 
# #### BARPLOT OF USED SPACE FOR DELTA TIME INSIDE
# ggplot(Mean_SpaceMeta_Rep, aes(x=size_treat, y=delta_time_outside,fill=Ant_status))+
#   geom_errorbar( aes(x=size_treat,ymin=delta_time_outside-SE_delta_time_outside, ymax=delta_time_outside+SE_delta_time_outside),position=position_dodge2(width=0.8, preserve = "single"))+
#   geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
#   #facet_wrap(~Exposed) +
#   scale_fill_manual(values = custom_palette) + # Apply custom color palette
#   STYLE_NOVIR +
#   labs(x="")
#   #+
#   #labs(title= "Space Use, colony means w/ SE",
#   #     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))
# 
# 
# 
# 
# 
# #####################################################
# 
# 

individual_data <- read.table("/media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",header=T,stringsAsFactors = F, sep=" ")
metadata <- read.table("/media/cf19810/DISK4/EXP1_base_analysis/EXP_summary_data/Metadata_Exp1_2021_2023-02-27.txt",header=T,stringsAsFactors = F, sep=",")

metadata$tag <- metadata$antID
individual_data <- merge(individual_data,metadata[,c("REP_treat","tag","AntTask")],all.x = T)
individual_data <- individual_data[which(!is.na(individual_data$AntTask) & individual_data$AntTask != "queen"), ]

individual_data$Ant_status <- paste(individual_data$status, individual_data$AntTask  ,sep=" ")

#aggregate per time_of_day
individual_data <- aggregate(prop_time_outside ~ REP_treat + Ant_status + tag + treatment + period + time_hour, FUN=mean, na.rm=T, na.action=na.pass, individual_data)

library(dplyr)
library(tidyr)


# Rename and order by name
individual_data$treatment <- as.factor(individual_data$treatment)
levels(individual_data$treatment)[levels(individual_data$treatment)=="control.big"] <- "Big Sham"
levels(individual_data$treatment)[levels(individual_data$treatment)=="pathogen.big"] <- "Big Pathogen"
levels(individual_data$treatment)[levels(individual_data$treatment)=="control.small"] <- "Small Sham"
levels(individual_data$treatment)[levels(individual_data$treatment)=="pathogen.small"] <- "Small Pathogen"
individual_data$treatment <- factor(individual_data$treatment, levels = sort(levels(individual_data$treatment)))

individual_data$Ant_status <- factor(individual_data$Ant_status, levels = c("untreated nurse" ,"untreated forager","treated nurse" ))


#difference in proportion time outside
individual_data_diff  <- individual_data %>%
  pivot_wider(names_from = period, values_from = prop_time_outside) %>%
  mutate(diff_prop_time_outside = post - pre) %>%
  select(REP_treat, Ant_status, treatment, diff_prop_time_outside,tag)

#aggregg for boxplot - ants means, by colony
#mean by colony, AntTask, exposed, size_treat
Mean_data_diff <- aggregate(diff_prop_time_outside ~ REP_treat + Ant_status  + treatment, FUN=mean, na.rm=T, na.action=na.pass, individual_data_diff)
SD_data_diff <- aggregate(diff_prop_time_outside ~ REP_treat + Ant_status  + treatment, FUN=sd, na.rm=T, na.action=na.pass, individual_data_diff); colnames(SD_data_diff) [match("delta_time_outside",colnames(SD_data_diff))] <- "SD_delta_time_outside"
#merge dfs
Mean_data_diff <- plyr::join (x = Mean_data_diff , y=SD_data_diff, type = "right", match = "all")

### for barplot - colony means
Mean_data_diff_Rep <- aggregate(diff_prop_time_outside ~ Ant_status + treatment, FUN=mean, na.rm=T, na.action=na.pass, Mean_data_diff)
SE_data_diff_Rep   <- aggregate(diff_prop_time_outside ~ Ant_status + treatment, FUN=std.error, na.rm=T, na.action=na.pass, Mean_data_diff); colnames(SE_data_diff_Rep) [match("diff_prop_time_outside",colnames(SE_data_diff_Rep))] <- "SE_delta_time_outside"
#merge dfs
Mean_data_diff_Rep <- plyr::join (x = Mean_data_diff_Rep , y=SE_data_diff_Rep, type = "right", match = "all")


# Create the boxplot
ggplot(Mean_data_diff, aes(x=treatment, y=diff_prop_time_outside, color=Ant_status)) +
  geom_boxplot(position = position_dodge(width = 0.8, preserve = "single")) +
  scale_color_manual(values = custom_palette) + # Apply custom color palette
  STYLE_NOVIR +
  labs(x="")+
  theme(legend.position = "bottom")
#labs(title= "Space Use, ant means by colony",
#     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))




#### BARPLOT OF USED SPACE FOR DELTA TIME INSIDE
ggplot(Mean_data_diff_Rep, aes(x = treatment, y = diff_prop_time_outside, fill = Ant_status)) +
  geom_errorbar(aes(x = treatment, ymin = diff_prop_time_outside - SE_delta_time_outside, ymax = diff_prop_time_outside + SE_delta_time_outside),
                position = position_dodge2(width = 0.8, preserve = "single")) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(values = custom_palette) +
  STYLE_NOVIR +
  labs(x = "") +
  theme(legend.position = "bottom")
#+
#labs(title= "Space Use, colony means w/ SE",
#     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))















library(dplyr)


min_sample <- individual_data
# min_sample$time_hours <- NULL
# min_sample$time_of_day <- NULL
# Create a minimal reproducible sample
min_sample <- min_sample %>%
  group_by(Ant_status, treatment,period) %>%
  slice_sample(n = 1) %>%
  ungroup()
dput(as.data.frame(min_sample[,c("Ant_status", "treatment", "period","prop_time_outside", "REP_treat")]))





min_sample <- iterations_data %>%
  group_by(i,x0,k,L,IMPUTATION) %>% #"x0","k","L","IMPUTATION"
  slice_sample(n = 1) %>%
  ungroup()
min_sample$rel_conc <- round(min_sample$rel_conc,1)
min_sample$rel_conc_imputed <- round(min_sample$rel_conc_imputed,1)
min_sample$x0 <- round(min_sample$x0,1)


dput(as.data.frame(min_sample[,c("i","Ant_status", "Treatment","rel_conc","rel_conc_imputed","LogFun_Category","propNA","x0","k","L","IMPUTATION")]))