
# # To Do
### FIND A WAY TO COPY RELEVANT DATA INFO (E.G. TREATMENT, ETC,) WHEN EXPAND GRIDDING!
# SHOULD JOIN BE USED INSTEAD OF MERGE? ANY THER ALTERNATIVE?

### FIX delta calcs  TO BE OPERATED DIRECTLY ON "inferred_ByAnt"!!!!!!!!!!!

### GIVE CLEARER NAMES TO THINGS:
# - MEANS OF EACH ANT IN REP: USED FOR DATA ANALYSIS AND PLOTTING (now called "ByAnt")
# - MEANS OF EACH REP: USED FOR PLOTTING (NOW CALLED "FOR PLOTS", SHOULD BE xANTxREP AND xREP ?)


####################################################################################
#### THIS SCRIPT CONTAINS:
#### DATA MANIPULATION TO PLOT INFERRED GROOMING
#### extra: SOME BASE PLOTS USED IN THE IUSSI SAN DIEGO 2022 PRESENTATION (MOVE THEM?)
####################################################################################

#### LIBRARIES
gc()
mallinfo::malloc.trim(0L)
# plot aesthetics
library(showtext)
font_add_google("Montserrat")
showtext_auto() # must be called to indicate that showtext is going to be automatically invoked to draw text whenever a plot is created.

library(FortMyrmidon)
library("ggplot2")
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)
library("viridis")
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# stats
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) # post-hoc comparisons
library(e1071) # calc skewness and other stuff
library(lawstat) # for levene test (homogeneity of variance)
library(lmPerm) # for permutations
library(lmerTest)
library(GLMMadaptive)
library(censReg) # random effect in tobit models in R using the censReg package?
library(tidyverse)
library(plm) #  pdata.frame()
library(mgcv) # fit Generalised additive mixed models
library(ggeffects)
library(itsadug)
library(performance)
library(remotes)
# install_version("MuMIn", "1.46.0")
library(MuMIn)
library(multcomp)
library(multcompView)


#### FUNCTIONS
### source function scripts
print("Loading functions and libraries...")
source(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/FUNCTIONS_Analysis_and_Styling.R", sep = "/"))

# list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n - 1)
    c(res, add)
  } else {
    res
  }
}

### DIRECTORIES
WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <- paste(WORKDIR, "Data/Adriano_extracted_groomings_final", sep = "/")
METADATADIR <- paste(WORKDIR, "EXP_summary_data", sep = "/")

# Base files
metadata <- read.table(paste(METADATADIR, "/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")

### TURN stats report on/off
REPORT <- TRUE
Report_stats <- list()

### create files for time bins
TIME_BINS <- FALSE

### Analyse either only treated or all
only_treated <- TRUE

### EXTRAS
Presentation_plots <- FALSE # script for single interactions heatmaps

# plot saving folder
if (only_treated == T) {
  save_dir_plots <- paste0(DATADIR, "/Grooming_plots/")
} else {
  save_dir_plots <- paste0(DATADIR, "/Grooming_plots_ALL_ANTS/")
}



tokeep <- ls()

###########################################################################################################
###### ALIGNMENT OF FILES STOP/START TIMES
###########################################################################################################

if (!file.exists(paste0(WORKDIR, "/Data/inferred_groomings_ALL_FINAL.txt"))) {
  ## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED!
  # window_shift <- 60*20 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition
  # window_shift_UNIX <- as.POSIXct(window_shift,origin = "1970-01-01",tz = "GMT")

  ### read return times
  returnTimes_Vasudha <- read.csv(paste(WORKDIR, "Data", "Exposed_nurses_R3SP_R9SP_RETURN_EXP_TIME_ZULU.csv", sep = "/"), sep = ",")
  returnTimes_Adriano <- read.csv(paste(WORKDIR, "Data", "Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU.csv", sep = "/"), sep = ",")
  returnTimes <- rbind(
    returnTimes_Vasudha[, c("REP_treat", "ReturnExposed_time")],
    returnTimes_Adriano[, c("REP_treat", "ReturnExposed_time")],
    structure(list(REP_treat = c("R9SS", "R9BS"), ReturnExposed_time = c(
      "2021-04-26T11:22:20.913120342Z",
      "2021-04-26T11:17:05.260434183Z"
    )), row.names = 1:2, class = "data.frame")
  )
  # write a summed up file with return times and upload it on the EXP_summary_folder
  write.table(returnTimes, file = paste(WORKDIR, "/Data/ALL_COLONIES_RETURN_EXP_TIME_ZULU.txt", sep = ""), append = F, col.names = T, row.names = F, quote = T, sep = ",")

  # Convert Zulu time to POSIXct
  returnTimes$ReturnExposed_time_sec <- as.POSIXct(returnTimes$ReturnExposed_time, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "GMT")
  # Convert POSIXct to seconds since epoch start
  returnTimes$ReturnExposed_time_sec <- as.numeric(returnTimes$ReturnExposed_time_sec)
  # #calc start and end time
  returnTimes$start_time_sec <- returnTimes$ReturnExposed_time_sec - 60 * 60 * 27
  returnTimes$end_time_sec <- returnTimes$ReturnExposed_time_sec + 60 * 60 * 24
  returnTimes$preReturn_gap_time_sec <- returnTimes$ReturnExposed_time_sec - 60 * 60 * 3

  ### DIRECTORIES
  # List all text files in the working directory
  filenames <- list.files(DATADIR, pattern = "\\whole_experiment.txt$")
  filenames <- filenames[grep("nongroom", filenames, invert = TRUE)]
  filesdir <- paste(DATADIR, filenames, sep = "/")
  # Read every text file with header, skipping the 1st row.
  result <- lapply(filesdir, function(x) read.table(x, header = T, stringsAsFactors = F, sep = ","))
  inferred <- do.call(rbind, result)

  inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX, format = "%Y-%m-%d %H:%M:%OS", origin = "1970-01-01", tz = "GMT")
  inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX, format = "%Y-%m-%d %H:%M:%OS", origin = "1970-01-01", tz = "GMT")

  # base file info
  # TREATMENT
  inferred$TREATMENT <- substr(inferred$REPLICATE, (nchar(inferred$REPLICATE) + 1) - 2, nchar(inferred$REPLICATE))
  # get rep-treat to match with exp filenames
  inferred$REP_treat <- sub(".*\\/", "", inferred$REPLICATE)

  # #inferred$return_time <- NA
  # inferred$end_time <- NA
  # inferred$exposed <- "no"
  # inferred$dead <- "no"
  # #inferred$time_since_start <- NA

  # ### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
  # #list subdirectories in parent folder EXPERIMENT_DATA
  # files_list <- list.dirs.depth.n("/media/cf19810/DISK3/ADRIANO/EXPERIMENT_DATA", n = 1)
  # #select REP folders
  # files_list <- files_list[grep("REP",files_list)]
  #
  # Reps_N_exposed <- data.frame(REP_treat= unique(inferred$REP_treat), N_ants = NA)

  # Add returnTimes to inferred
  # Convert data.frames to data.tables for faster processing
  returnTimes_dt <- as.data.table(returnTimes)
  inferred_dt <- as.data.table(inferred)
  # Find common column names
  common_col_names <- intersect(names(returnTimes_dt), names(inferred_dt))
  # Perform fast merge using data.table
  setkeyv(returnTimes_dt, common_col_names) # Set keys for returnTimes_dt
  setkeyv(inferred_dt, common_col_names) # Set keys for inferred_dt
  # Merge data.tables
  inferred_dt <- inferred_dt[returnTimes_dt, nomatch = 0]

  # # Add metadata to inferred
  # metadata_dt <- as.data.table(metadata)
  # # Find common column names
  # DIFFERENT STRUCTURE, WON'T WORK AS IN INFERRED THRE ARE TWO ANTS
  # common_col_names <- intersect(names(metadata_dt), names(inferred_dt))
  # # Perform fast merge using data.table
  # setkeyv(metadata_dt, common_col_names) # Set keys for metadata_dt
  # setkeyv(inferred_dt, common_col_names) # Set keys for inferred_dt
  # # Merge data.tables
  # inferred_dt <- inferred_dt[metadata_dt, nomatch = 0]

  # Optionally, convert the result back to data.frame
  inferred <- as.data.frame(inferred_dt)

  # Create the new column "period" using ifelse()
  inferred$PERIOD <- with(inferred, ifelse(
    T_start_sec >= ReturnExposed_time_sec, "post",
    ifelse(T_start_sec >= preReturn_gap_time_sec & T_start_sec < ReturnExposed_time_sec, "GAP", "pre")
  ))

  # quick visual check
  ggplot(inferred[which(inferred$REP_treat == "R9BS"), ], aes(x = T_start_sec, fill = PERIOD)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 1000) +
    geom_vline(aes(xintercept = preReturn_gap_time_sec), linetype = "dashed", color = "blue") +
    geom_vline(aes(xintercept = ReturnExposed_time_sec), linetype = "dashed", color = "red") +
    labs(
      title = "Simplified Visualization of Period Assignment",
      x = "T_start_sec",
      y = "Count",
      fill = "PERIOD"
    ) +
    theme_minimal()

  ### save the output
  write.table(inferred, file = paste(WORKDIR, "/Data/inferred_groomings_ALL_FINAL.txt", sep = ""), append = F, col.names = T, row.names = F, quote = T, sep = ",")
  # write.table(Reps_N_exposed,file=paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
}

#######################################################
# start from here

####################################################################################
########### PREPARE DATA FOR ANALYSES AND PLOTTING #################################
####################################################################################

###### AGGREGATE ALL VALUES FOR PRE.POST  #####
# LOAD
inferred <- read.table(paste(WORKDIR, "/Data/inferred_groomings_ALL_FINAL.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
Reps_N_exposed <- read.table(paste(METADATADIR, "/N_ants_exposed_xREP.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
metadata <- metadata[, c("REP_treat", "antID", "Exposed", "IsAlive", "AntTask")]
metadata <- metadata[which(metadata$AntTask != "queen"), ]

# remove EXP_GAP data
inferred <- inferred[which(inferred$PERIOD != "GAP"), ]
# remove data without zone
inferred <- inferred[which(inferred$ant1.zones != 0), ]
inferred <- inferred[which(inferred$ant2.zones != 0), ]
# ## add a common time to all rows
# # Negative before exposure, positive after
inferred$time_stop_since_treat <- (inferred$T_stop_sec - inferred$ReturnExposed_time_sec)
# add count column
inferred$Count_byAnt <- 1

# N of reps
table(str_sub(unique(inferred$REP_treat), -2, -1))

# rename columns for merging
metadata$antID <- paste0("ant_", metadata$antID)

#######################################################################
# dput(head(metadata))
# dput(head(inferred[,c("REP_treat","Act_Name","Rec_Name","frame_start","frame_stop","duration" )]))

# Rename columns in metadata
colnames(metadata)[colnames(metadata) == "Exposed"] <- "Exposed_Act"
colnames(metadata)[colnames(metadata) == "IsAlive"] <- "IsAlive_Act"
colnames(metadata)[colnames(metadata) == "AntTask"] <- "AntTask_Act"

# Merge based on Act_Name
Meta <- merge(inferred, metadata, by.x = c("REP_treat", "Act_Name"), by.y = c("REP_treat", "antID"), all.x = TRUE)

# Rename columns in metadata for Rec_Name
colnames(metadata)[colnames(metadata) == "Exposed_Act"] <- "Exposed_Rec"
colnames(metadata)[colnames(metadata) == "IsAlive_Act"] <- "IsAlive_Rec"
colnames(metadata)[colnames(metadata) == "AntTask_Act"] <- "AntTask_Rec"

# Merge based on Rec_Name
Meta <- merge(Meta, metadata, by.x = c("REP_treat", "Rec_Name"), by.y = c("REP_treat", "antID"), all.x = TRUE)

# Rename columns in metadata to original
colnames(metadata)[colnames(metadata) == "Exposed_Rec"] <- "Exposed"
colnames(metadata)[colnames(metadata) == "IsAlive_Rec"] <- "IsAlive"
colnames(metadata)[colnames(metadata) == "AntTask_Rec"] <- "AntTask"

# ### Get info from metadata
# Meta <- list(inferred,metadata)
# Meta <- Reduce(function(x, y) merge(x, y, all=TRUE), Meta)

# there are mismatched ants that appear in the grooming data but not in the metadata.
#  this is ~1% of ants that appear in the grooming but not in the metadata, remove
nrow(unique(Meta[which(is.na(Meta$AntTask_Act)), c("REP_treat", "Act_Name")]))
nrow(unique(Meta[which(is.na(Meta$AntTask_Rec)), c("REP_treat", "Rec_Name")]))

nrow(unique(Meta[which(!is.na(Meta$AntTask_Rec)), c("REP_treat", "Rec_Name")]))
nrow(unique(Meta[which(!is.na(Meta$AntTask_Act)), c("REP_treat", "Act_Name")]))

# remove mismatched ants
Meta <- Meta[which(!is.na(Meta$AntTask_Rec)), ]
Meta <- Meta[which(!is.na(Meta$AntTask_Act)), ]

# ensure removal of dead (should be 0)
# Meta <- Meta[which(Meta$IsAlive_Rec==TRUE),]

## expand grid on list of exposed per REP that have received no grooming in the specific timespan
# Meta_all_RecPer <- metadata %>%
#   group_by(REP_treat) %>%
#   tidyr::expand (antID,PERIOD= c("pre","post")) #specify as some colonies (REP_treat) don't have observations for pre to expand on
#


# expand grid and then filter b the subject ants
# Create a copy of metadata with period "pre"
metadata_pre <- metadata
metadata_pre$PERIOD <- "pre"
# Create a copy of metadata with period "post"
metadata_post <- metadata
metadata_post$PERIOD <- "post"
# Combine the two data frames using rbind
Meta_all_RecPer <- rbind(metadata_pre, metadata_post)

# select present receivers
Meta_all_RecPer <- Meta_all_RecPer[which(Meta_all_RecPer$Exposed == only_treated), ]
Meta <- Meta[which(Meta$Exposed_Rec == only_treated), ]

# remove extra vars not needed for merge
Meta_all_RecPer <- Meta_all_RecPer[, c("REP_treat", "antID", "PERIOD")]
colnames(Meta_all_RecPer)[colnames(Meta_all_RecPer) == "antID"] <- "Rec_Name"

# #select exposed only or non exposed
# if (only_treated==TRUE) {
#   print("Only treated nurses")
# }else{
#   print("Only untreated ants")
# }



### merge
Meta_all_combs <- list(Meta, Meta_all_RecPer)
Meta_all_combs <- Reduce(function(x, y) merge(x, y, all = TRUE), Meta_all_combs) # , by = c("REP_treat","PERIOD","Rec_Name")

## check if any data has NAs
# Meta_all_combs[is.na(Meta_all_combs$PERIOD),]

# #add 0 counts and durations
Meta_all_combs[is.na(Meta_all_combs$Act_Name), "Count_byAnt"] <- 0
Meta_all_combs[is.na(Meta_all_combs$Act_Name), "duration"] <- 0

# OVERWRITE
inferred <- Meta_all_combs

# Relevel Exposed
# for the moment, to avoid errors in subsettings etc, keep the new exposed var aside
inferred$Exposed <- as.factor(inferred$Exposed_Rec)
levels(inferred$Exposed)[levels(inferred$Exposed) == "TRUE"] <- "treated"
levels(inferred$Exposed)[levels(inferred$Exposed) == "FALSE"] <- "untreated"
# Create new Status category
inferred$Ant_status <- paste(inferred$Exposed, inferred$AntTask_Rec)
# remove the 6 mismatched non-treated receivers -empty grid - (why there?) created when pasting
inferred <- inferred[which(!inferred$Ant_status == "NA NA"), ]

# Rename by name
inferred$TREATMENT <- as.factor(inferred$TREATMENT)
levels(inferred$TREATMENT)[levels(inferred$TREATMENT) == "BS"] <- "Big Sham"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT) == "BP"] <- "Big Pathogen"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT) == "SS"] <- "Small Sham"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT) == "SP"] <- "Small Pathogen"


warning("make more summary tables! maybe with report(). use same functions as in immune_genes tables code")
as.data.frame(table(inferred$Ant_status))



###### AGGREGATE VALUES FOR PRE.POST: FULL DATASET | FOR STATISTICS #####
### AGGREGATE TO THE LEVEL OF ANT

## count the number of observations by ant and get mean
inferred_count_summary <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + REP_treat + Rec_Name + Ant_status, FUN = sum, na.action = na.pass, inferred)
## calculate mean durations for REP
inferred_dur_summary <- aggregate(duration ~ PERIOD + TREATMENT + REP_treat + Rec_Name + Ant_status, FUN = mean, na.rm = T, na.action = na.pass, inferred)
# sum by ant (SUM DUR MAY BE MORE INFORMATIVE AS THE GROOMING DETECTION MAY RESULT FRAGMENTED)
inferred_SUM <- aggregate(duration ~ PERIOD + TREATMENT + REP_treat + Rec_Name + Ant_status, FUN = sum, na.rm = T, na.action = na.pass, inferred)
names(inferred_SUM)[names(inferred_SUM) == "duration"] <- "SUM_duration"

# MERGE - REP_treatments
# THIS DATA WILL BE USED FOR STATS!
# MAKE SURE THAT THERE ARE TWO VALUES PER ANT (pre and post treatment)
inferred_ByAnt <- list(inferred_count_summary, inferred_dur_summary, inferred_SUM)
inferred_ByAnt <- Reduce(function(x, y) merge(x, y, all = TRUE), inferred_ByAnt)
# inferred_count_summ1    <-  plyr::join(x=inferred_count_summary, y=inferred_dur_summary, type = "full", match = "all")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### FIX THIS TO BE OPERATED DIRECTLY ON "inferred_ByAnt"!!!!!!!!!!!

## CALCULATING THE DIFFERENCES POST MINUS AFTER REQUIRES SUBTRACTING MEANS, THEREFORE SE HAS TO BE CALCULATED ACCORDINGLY
# https://rpubs.com/brouwern/SEdiff2means

## diff after relative to before
# StdError_diff_between_means

# Formula for POOLED standard deviation (to calculate the SE of the difference between means)
## Note the formulas squares SD to get variance
var.pooled <- function(N1, N2, SD1, SD2) {
  (N1 * SD1^2 + N2 * SD2^2) / (N1 + N2)
}
# Standard error of difference
## Note that this uses sample size, NOT degrees of freedom (N)
SE.diff <- function(var.pool, n1, n2) {
  sqrt(var.pool * (1 / n1 + 1 / n2))
}


## SAME BUT FOR DATA ANALYSIS :

# #use the functions above to calculate the pooled SD and the N samples!
# Counts_by_Behaviour_AllCombos1_DELTA <- Counts_by_Behaviour_AllCombos1 %>%
#   group_by(REP_treat) %>%
#   dplyr::summarise(Count_byAnt = Count_byAnt[match("post", PERIOD)] - Count_byAnt[match("pre", PERIOD)],
#                    duration = duration[match("post", PERIOD)] - duration[match("pre", PERIOD)],
#                    SUM_duration = SUM_duration[match("post", PERIOD)] - SUM_duration[match("pre", PERIOD)],
#                    TREATMENT = TREATMENT
#                    # SE_Count_byAnt = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_Count_byAnt[match("post", PERIOD)], SD_Count_byAnt[match("pre", PERIOD)])*(1/N_Count_REP[match("post", PERIOD)] + 1/N_Count_REP[match("pre", PERIOD)])),
#                    # SE_duration = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_duration[match("post", PERIOD)], SD_duration[match("pre", PERIOD)])*(1/N_Count_REP[match("post", PERIOD)] + 1/N_Count_REP[match("pre", PERIOD)])),
#                    # SE_SUM_duration = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_SUM_duration[match("post", PERIOD)], SD_SUM_duration[match("pre", PERIOD)])*(1/N_Count_REP[match("post", PERIOD)] + 1/N_Count_REP[match("pre", PERIOD)]))
#                    # SD_Count_byAnt[match("post", PERIOD)] - SD_Count_byAnt[match("pre", PERIOD)],
#                    # SE_duration = SE_duration[match("post", PERIOD)] - SE_duration[match("pre", PERIOD)],
#   )
#
# Counts_by_Behaviour_AllCombos1_DELTA <- dplyr::distinct(Counts_by_Behaviour_AllCombos1_DELTA)
# Counts_by_Behaviour_AllCombos1_DELTA <- as.data.frame(Counts_by_Behaviour_AllCombos1_DELTA)



###### AGGREGATE VALUES FOR PRE.POST: FULL DATASET | FOR PLOTS #####
### AGGREGATE TO THE LEVEL OF NEST (REP_treat)

# calculate mean by REP
inferred_count_summ1 <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + REP_treat + Ant_status, FUN = mean, na.action = na.pass, inferred_ByAnt)

# alternative expand grid for paired columns (REP_treat and TREATMENT)
combos1 <- unique(inferred[, c("REP_treat", "TREATMENT")])
combos1 <- combos1[!is.na(combos1$TREATMENT), ]
combos2 <- unique(inferred[, "PERIOD", drop = FALSE])
combos3 <- merge(combos1, combos2)
Counts_by_Behaviour_AllCombos1 <- merge(inferred_count_summ1, combos3, all = TRUE)

## replace the NAs with 0 counts
Counts_by_Behaviour_AllCombos1[is.na(Counts_by_Behaviour_AllCombos1)] <- 0

# add N events and SD to calculate the Std Err of the mean differences
infer_Nevents <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + Ant_status, FUN = length, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
colnames(infer_Nevents)[match("Count_byAnt", colnames(infer_Nevents))] <- "N_Count_REP"
infer_SD <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + Ant_status, FUN = sd, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
colnames(infer_SD)[match("Count_byAnt", colnames(infer_SD))] <- "SD_Count_byAnt"
colnames(infer_SD)[match("duration", colnames(infer_SD))] <- "SD_duration"
colnames(infer_SD)[match("SUM_duration", colnames(infer_SD))] <- "SD_SUM_duration"
## finally, get the mean & S.E. for each behav before/after  for barplots
infer_MEAN <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + Ant_status, FUN = mean, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
infer_SE <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + Ant_status, FUN = std.error, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
# rename cols
colnms <- c("Count_byAnt", "duration", "SUM_duration")
colnames(infer_MEAN)[match(colnms, colnames(infer_MEAN))] <- paste0("Mean_", colnms)
colnames(infer_SE)[match(colnms, colnames(infer_SE))] <- paste0("SE_", colnms)

# MERGE everything
infer_full <- list(infer_Nevents, infer_SD, infer_MEAN, infer_SE)
infer_full <- Reduce(function(x, y) merge(x, y, all = TRUE), infer_full)
infer_full[is.na(infer_full)] <- 0
# reorder levels to fix strange behaviour in plotting
infer_full$TREATMENT <- factor(infer_full$TREATMENT, levels = c("Big Pathogen", "Big Sham", "Small Pathogen", "Small Sham"))
# infer_full$TREATMENT <- factor(infer_full$TREATMENT, levels = c("Big Pathogen","Big Sham","Small Pathogen","Small Sham"))
# reorder levels of period
infer_full$PERIOD <- factor(infer_full$PERIOD, levels = c("pre", "post"))


### CALCULATING THE DIFFERENCES POST MINUS AFTER (DELTA)
# use the functions above to calculate the pooled SD and the N samples!
infer_full_DELTA <- infer_full %>%
  group_by(TREATMENT) %>%
  dplyr::summarise(
    Mean_Count_byAnt = Mean_Count_byAnt[match("post", PERIOD)] - Mean_Count_byAnt[match("pre", PERIOD)],
    Mean_duration = Mean_duration[match("post", PERIOD)] - Mean_duration[match("pre", PERIOD)],
    Mean_SUM_duration = Mean_SUM_duration[match("post", PERIOD)] - Mean_SUM_duration[match("pre", PERIOD)],
    SE_Count_byAnt = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_Count_byAnt[match("post", PERIOD)], SD_Count_byAnt[match("pre", PERIOD)]) * (1 / N_Count_REP[match("post", PERIOD)] + 1 / N_Count_REP[match("pre", PERIOD)])),
    SE_duration = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_duration[match("post", PERIOD)], SD_duration[match("pre", PERIOD)]) * (1 / N_Count_REP[match("post", PERIOD)] + 1 / N_Count_REP[match("pre", PERIOD)])),
    SE_SUM_duration = sqrt(var.pooled(N_Count_REP[match("post", PERIOD)], N_Count_REP[match("pre", PERIOD)], SD_SUM_duration[match("post", PERIOD)], SD_SUM_duration[match("pre", PERIOD)]) * (1 / N_Count_REP[match("post", PERIOD)] + 1 / N_Count_REP[match("pre", PERIOD)]))
    # SD_Count_byAnt[match("post", PERIOD)] - SD_Count_byAnt[match("pre", PERIOD)],
    # SE_duration = SE_duration[match("post", PERIOD)] - SE_duration[match("pre", PERIOD)],
  )


###### AGGREGATE VALUES FOR PRE.POST: TIME BINS | FOR STATS #####
### AGGREGATE TO THE LEVEL OF ANT
warning("THERE IS NO ANT_STATUS IMPLEMENTED HERE")


if (TIME_BINS) {

  # time bins for plotting
  time.break <- c("h4", "hour", "10min", "3h") # ,"10min"     "h24",

  for (TIME in time.break) {
    # clean all
    rm(list = setdiff(ls(), c("inferred", "infer_full_DELTA", "Counts_by_Behaviour_AllCombos1_DELTA", "Reps_N_exposed", "infer_full", "WORKDIR", "DATADIR", "TIME", "time.break", "tokeep", "inferred_ByAnt", "time_of_day", tokeep)))
    # create local copy
    inferred_bin <- inferred

    # add time.breaks to dataframe
    if (TIME == "hour") {
      inferred_bin$timespan <- inferred_bin$time_stop_since_treat / 3600
    } else if (TIME == "3h") {
      inferred_bin$timespan <- inferred_bin$time_stop_since_treat / 10800
    } else if (TIME == "10min") {
      inferred_bin$timespan <- inferred_bin$time_stop_since_treat / 600
    } else if (TIME == "h4") {
      inferred_bin$timespan <- inferred_bin$time_stop_since_treat / 14400 # create interval of 4h
    }

    inferred_bin$timespan <- floor(inferred_bin$timespan)

    #######################################
    # #expand grid of all possible combs of Rec_Name and period (pre-post) within group
    # Meta_all_RecTime <- inferred_bin %>%
    #   group_by(REP_treat) %>%
    #   tidyr::expand (Rec_Name,timespan)
    # Meta_all_RecTime <- Meta_all_RecTime[which(!is.na(Meta_all_RecTime$timespan)),]
    #
    # ### merge
    # Meta_all_combsH <- list(inferred_bin,Meta_all_RecTime)
    # Meta_all_combsH <- Reduce(function(x, y) merge(x, y, all=TRUE), Meta_all_combsH)
    # #clean
    # Meta_all_combsH <- Meta_all_combsH[which(!is.na(Meta_all_combsH$Rec_Name)),]
    # Meta_all_combsH <- Meta_all_combsH[which(!is.na(Meta_all_combsH$timespan)),]
    # #add 0 counts and durations
    # Meta_all_combsH[is.na(Meta_all_combsH$Act_Name),"Count_byAnt"] <- 0
    # Meta_all_combsH[is.na(Meta_all_combsH$Act_Name),"duration"]    <- 0
    # #Re-add missing treatments
    # Meta_all_combsH$TREATMENT <- substr(Meta_all_combsH$REP_treat,(nchar(Meta_all_combsH$REP_treat)+1)-2,nchar(Meta_all_combsH$REP_treat))
    # #Re-add missing PERIODs
    # #since these are times in seconds, it is ok to  have strict ">" instead of ">="
    # Meta_all_combsH$PERIOD <- ifelse(Meta_all_combsH$timespan >= 0, "post", "pre")
    #


    # Rename by name
    inferred_bin$TREATMENT <- as.factor(inferred_bin$TREATMENT)
    levels(inferred_bin$TREATMENT)[levels(inferred_bin$TREATMENT) == "BS"] <- "Big Sham"
    levels(inferred_bin$TREATMENT)[levels(inferred_bin$TREATMENT) == "BP"] <- "Big Pathogen"
    levels(inferred_bin$TREATMENT)[levels(inferred_bin$TREATMENT) == "SS"] <- "Small Sham"
    levels(inferred_bin$TREATMENT)[levels(inferred_bin$TREATMENT) == "SP"] <- "Small Pathogen"
    ################################################

    ## count the number of observations by ant and get mean
    inferred_count_bin_summary <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + timespan + REP_treat + Rec_Name, FUN = sum, na.action = na.pass, inferred_bin)
    ## calculate mean durations for REP
    inferred_dur_bin_summary <- aggregate(duration ~ PERIOD + TREATMENT + timespan + REP_treat + Rec_Name, FUN = mean, na.rm = T, na.action = na.pass, inferred_bin)
    # sum by ant (SUM DUR MAY BE MORE INFORMATIVE AS THE GROOMING DETECTION MAY RESULT FAGMENTED)
    inferred_bin_SUM <- aggregate(duration ~ PERIOD + TREATMENT + timespan + REP_treat + Rec_Name, FUN = sum, na.rm = T, na.action = na.pass, inferred_bin)
    names(inferred_bin_SUM)[names(inferred_bin_SUM) == "duration"] <- "SUM_duration"

    # MERGE - REP_treatments
    inferred_bin_ByAnt <- list(inferred_count_bin_summary, inferred_dur_bin_summary, inferred_bin_SUM)
    inferred_bin_ByAnt <- Reduce(function(x, y) merge(x, y, all = TRUE), inferred_bin_ByAnt)
    # inferred_count_bin_summ1    <-  plyr::join(x=inferred_count_bin_summary, y=inferred_dur_bin_summary, type = "full", match = "all")

    ### expand grid
    expected <- metadata %>%
      group_by(REP_treat) %>%
      tidyr::expand(Rec_Name, timespan = unique(inferred_bin$timespan)) # specify as some colonies (REP_treat) don't have observations for pre to expand on

    inferred_bin_ByAnt <- merge(expected, inferred_bin_ByAnt[which(names(inferred_bin_ByAnt) != "TREATMENT")], all.x = T)

    inferred_bin_ByAnt[which(is.na(inferred_bin_ByAnt$duration)), "duration"] <- 0
    inferred_bin_ByAnt[which(is.na(inferred_bin_ByAnt$SUM_duration)), "SUM_duration"] <- 0
    inferred_bin_ByAnt[which(is.na(inferred_bin_ByAnt$Count_byAnt)), "Count_byAnt"] <- 0
    inferred_bin_ByAnt$PERIOD <- ifelse(inferred_bin_ByAnt$timespan >= 0, "post", "pre")
    inferred_bin_ByAnt <- merge(inferred_bin_ByAnt, unique(inferred_bin[c("REP_treat", "TREATMENT")]))



    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ############################################
    # pathogen−induced changes relative to sham−induced changes

    ## copy from infer_full
    # BOTH FOR DATA ANALYSIS AND FOR PLOTS


    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


    ###### AGGREGATE VALUES FOR PRE.POST: TIME BINS | FOR STATS #####
    ### AGGREGATE TO THE LEVEL OF NEST (REP_treat)

    # calculate mean by REP
    inferred_count_bin_summ1 <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + timespan + REP_treat, FUN = mean, na.action = na.pass, inferred_bin_ByAnt)

    # alternative expand grid for paired columns (REP_treat and TREATMENT)
    ### IF NEEDED, MODIFY combos_bin2 TO  unique(inferred_bin[,c("timespan","PERIOD")])
    combos_bin1 <- unique(inferred_bin[, c("REP_treat", "TREATMENT")])
    combos_bin2 <- unique(inferred_bin[, "timespan", drop = FALSE])
    combos_bin3 <- merge(combos_bin1, combos_bin2)
    Counts_by_Behaviour_AllCombos1 <- merge(inferred_count_bin_summ1, combos_bin3, all = TRUE)

    ## replace the NAs with 0 counts
    Counts_by_Behaviour_AllCombos1[is.na(Counts_by_Behaviour_AllCombos1)] <- 0

    # add N events and SD to calculate the Std Err of the mean differences
    infer_bin_Nevents <- aggregate(Count_byAnt ~ PERIOD + timespan + TREATMENT, FUN = length, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
    colnames(infer_bin_Nevents)[match("Count_byAnt", colnames(infer_bin_Nevents))] <- "N_Count_REP"
    infer_bin_SD <- aggregate(cbind(Count_byAnt, duration) ~ PERIOD + timespan + TREATMENT, FUN = sd, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
    colnames(infer_bin_SD)[match("Count_byAnt", colnames(infer_bin_SD))] <- "SD_Count_byAnt"
    colnames(infer_bin_SD)[match("duration", colnames(infer_bin_SD))] <- "SD_duration"
    ## finally, get the mean & S.E. for each behav before/after  for barplots
    infer_bin_MEAN <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + timespan + TREATMENT, FUN = mean, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
    infer_bin_SE <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + timespan + TREATMENT, FUN = std.error, na.rm = T, na.action = na.pass, Counts_by_Behaviour_AllCombos1)
    # rename cols
    colnms <- c("Count_byAnt", "duration", "SUM_duration")
    colnames(infer_bin_MEAN)[match(colnms, colnames(infer_bin_MEAN))] <- paste0("Mean_", colnms)
    colnames(infer_bin_SE)[match(colnms, colnames(infer_bin_SE))] <- paste0("SE_", colnms)

    # MERGE everything
    infer_bin_full <- list(infer_bin_Nevents, infer_bin_SD, infer_bin_MEAN, infer_bin_SE)
    infer_bin_full <- Reduce(function(x, y) merge(x, y, all = TRUE), infer_bin_full)
    infer_bin_full[is.na(infer_bin_full)] <- 0
    # reorder levels to fix strange behaviour in plotting
    infer_bin_full$TREATMENT <- factor(infer_bin_full$TREATMENT, levels = c("Big Pathogen", "Big Sham", "Small Pathogen", "Small Sham"))



    if (TIME == "hour") {
      ## assign a time of the day!
      Time_dictionary <- data.frame(timespan = -36:35, time_of_day = rep(0:23, 3))
      infer_bin_full <- left_join(infer_bin_full, Time_dictionary, by = "timespan")
      inferred_bin_ByAnt <- left_join(inferred_bin_ByAnt, Time_dictionary, by = "timespan")
    } else if (TIME == "3h") {
      # Time_dictionary <- data.frame(timespan= -9:8, time_of_day= rep(seq(0, 23, by = 3),3))
      Time_dictionary <- data.frame(timespan = -36:35, time_of_day = rep(0:23, 3))
      Time_dictionary$timespan <- Time_dictionary$timespan / 3 # this works nicely
      infer_bin_full <- left_join(infer_bin_full, Time_dictionary, by = "timespan")
      inferred_bin_ByAnt <- left_join(inferred_bin_ByAnt, Time_dictionary, by = "timespan")
      infer_bin_full$timespan <- infer_bin_full$timespan * 3
      inferred_bin_ByAnt$timespan <- inferred_bin_ByAnt$timespan * 3
    } else if (TIME == "10min") {
      timespan_10 <- -216:215
      time_of_day_10 <- rep(seq(0, 24, by = 0.1666667), times = ceiling(length(timespan_10) / 144))
      Time_dictionary <- data.frame(timespan = timespan_10, time_of_day = time_of_day_10[1:length(timespan_10)])
      infer_bin_full <- left_join(infer_bin_full, Time_dictionary, by = "timespan")
      # convert notation from 10mins blocks to hours
      infer_bin_full$timespan <- infer_bin_full$timespan / 6
      inferred_bin_ByAnt <- left_join(inferred_bin_ByAnt, Time_dictionary, by = "timespan")
      inferred_bin_ByAnt$timespan <- inferred_bin_ByAnt$timespan / 6
    } else if (TIME == "h4") {
      Time_dictionary1 <- data.frame(timespan = c(-7:6), time_of_day = c(8, 12, 16, 20, 0, 4, 8, 12, 16, 20, 0, 4, 8, 12))
      infer_bin_full <- left_join(infer_bin_full, Time_dictionary1, by = "timespan")
      inferred_bin_ByAnt <- left_join(inferred_bin_ByAnt, Time_dictionary1, by = "timespan")
    }



    # save output
    if (TIME == "h4") {
      infer_bin_h4 <- infer_bin_full
      inferred_bin_h4_ByAnt <- inferred_bin_ByAnt
      tokeep <- c("tokeep", "infer_bin_h4", "inferred_bin_h4_ByAnt", tokeep)
    } else if (TIME == "3h") {
      infer_bin_3h <- infer_bin_full
      inferred_bin_3h_ByAnt <- inferred_bin_ByAnt # this is lacking the GAP
      tokeep <- c("tokeep", "infer_bin_1h", "inferred_bin_1h_ByAnt", tokeep)
    } else if (TIME == "hour") {
      ## add break for line plots
      # GAP <- expand.grid(PERIOD= "pre", timespan=c(-2,-1),time_of_day=c(10,11), TREATMENT=c("Small Pathogen","Small Sham","Big Pathogen","Big Sham"),
      #                  N_Count_REP = NA, SD_Count_byAnt = NA, SD_duration = NA, Mean_Count_byAnt=NA, Mean_duration=NA, Mean_SUM_duration = NA, SE_Count_byAnt=NA, SE_duration=NA, SE_SUM_duration=NA )
      # infer_bin_1h <- rbind(infer_bin_full,GAP)
      infer_bin_1h <- infer_bin_full
      inferred_bin_1h_ByAnt <- inferred_bin_ByAnt # this is lacking the GAP
      tokeep <- c("tokeep", "infer_bin_1h", "inferred_bin_1h_ByAnt", tokeep)
    } else if (TIME == "10min") {
      ## add break  for line plots
      ## redo proper breaks, it should be a repeat between -17 (first interval after -18 (18/6 = 3h)) and last before 0
      ## likely not needed anyway as it will be a barplot
      # GAP <- expand.grid(PERIOD= "pre", timespan=c(-3:0),time_of_day=Time_dictionary[Time_dictionary$timespan %in% c(-17:-1),"time_of_day"], TREATMENT=c("Small Pathogen","Small Sham","Big Pathogen","Big Sham"),
      #                   N_Count_REP = NA, SD_Count_byAnt = NA, SD_duration = NA, Mean_Count_byAnt=NA, Mean_duration=NA, Mean_SUM_duration = NA, SE_Count_byAnt=NA, SE_duration=NA, SE_SUM_duration=NA )
      # transforming GAP in new notation
      # GAP$timespan <- (GAP$timespan/6) + 0.1666667; GAP <- GAP[which(GAP$timespan<0),]
      # infer_bin_10min<- rbind(infer_bin_full,GAP)
      infer_bin_10min <- infer_bin_full
      inferred_bin_10min_ByAnt <- inferred_bin_ByAnt
      tokeep <- c("tokeep", "infer_bin_10min", "inferred_bin_10min_ByAnt", tokeep)
    }
  } # TIME

  # reorder levels
  infer_bin_1h$PERIOD <- factor(infer_bin_1h$PERIOD, levels = c("pre", "post"))
  infer_bin_3h$PERIOD <- factor(infer_bin_3h$PERIOD, levels = c("pre", "post"))
  infer_bin_h4$PERIOD <- factor(infer_bin_h4$PERIOD, levels = c("pre", "post"))
  infer_bin_10min$PERIOD <- factor(infer_bin_10min$PERIOD, levels = c("pre", "post"))

  inferred_bin_1h_ByAnt$PERIOD <- factor(inferred_bin_1h_ByAnt$PERIOD, levels = c("pre", "post"))
  inferred_bin_3h_ByAnt$PERIOD <- factor(inferred_bin_3h_ByAnt$PERIOD, levels = c("pre", "post"))
  inferred_bin_h4_ByAnt$PERIOD <- factor(inferred_bin_h4_ByAnt$PERIOD, levels = c("pre", "post"))
  inferred_bin_10min_ByAnt$PERIOD <- factor(inferred_bin_10min_ByAnt$PERIOD, levels = c("pre", "post"))

  # TRIMMED DATA- 4h by ANT
  inferred_bin_h4_ByAnt_trim <- inferred_bin_h4_ByAnt[which(inferred_bin_h4_ByAnt$time_of_day == 12), ]

  # TRIMMED DATA - 1h by REP
  infer_bin_1h_trim <- infer_bin_1h[which(infer_bin_1h$time_of_day >= 12 & infer_bin_1h$time_of_day <= 16), ]
  # remove anything after timespan = 4 to exclude next day!
  infer_bin_1h_trim <- infer_bin_1h_trim[which(infer_bin_1h_trim$timespan <= 4), ]

  # make clear that those are hours
  text.add <- ":00"
  infer_bin_1h_trim$time_of_day <- paste0(infer_bin_1h_trim$time_of_day, text.add)

  # TRIMMED DATA- 4h by REP
  infer_bin_h4_trim <- infer_bin_h4[which(infer_bin_h4$time_of_day == 12), ]
  # remove anything after timespan = 4 to exclude next day!
  infer_bin_h4_trim <- infer_bin_h4_trim[which(infer_bin_h4_trim$timespan <= 1), ]

  ###### DATA BY REP (mean on ants)
  inferred_bin_1h_ByRep <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat, FUN = mean, na.rm = T, na.action = na.pass, inferred_bin_1h_ByAnt)
  inferred_bin_3h_ByRep <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat, FUN = mean, na.rm = T, na.action = na.pass, inferred_bin_3h_ByAnt)
  inferred_bin_10min_ByRep <- aggregate(cbind(Count_byAnt, duration, SUM_duration) ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat, FUN = mean, na.rm = T, na.action = na.pass, inferred_bin_10min_ByAnt)
} # ID TIME_BINS


#####################################################################################
##################                 STATS             ################################
#####################################################################################

inferred_ByAnt$TREATMENT <- as.factor(inferred_ByAnt$TREATMENT)
inferred_ByAnt$PERIOD <- as.factor(inferred_ByAnt$PERIOD)
inferred_ByAnt$REP_treat <- as.factor(inferred_ByAnt$REP_treat)
inferred_ByAnt$Rec_Name <- as.factor(inferred_ByAnt$Rec_Name)

# Select variables for analysis
# numeric_variable_list  <- names(inferred_ByAnt) [ which (!names(inferred_ByAnt) %in% c( "PERIOD","TREATMENT","REP_treat","Rec_Name" ) )]
numeric_variable_list <- c("Count_byAnt", "duration", "SUM_duration")

for (VAR in numeric_variable_list) {

  # create Variable of interest, as bracket subsetting -inferred_ByAnt[,VAR]- is not liked by report()
  inferred_ByAnt$VARIABLE <- inferred_ByAnt[, VAR]

  print(paste0("########## DEPENDANT VAR:", VAR, " ##########"))
  # VAR <- "Count_byAnt"
  # hist(sqrt(VARIABLE))
  m1 <- lmer(sqrt(VARIABLE) ~ PERIOD * TREATMENT + (1 | REP_treat / Rec_Name), data = inferred_ByAnt) # the "/" is for the nesting #  + (1|time_of_day)
  # simplify model if the interaction is not significant
  m1 <- simplify_model(m1)

  output_lmer(m1)

  qqnorm(residuals(m1))
  qqline(residuals(m1))
  hist(residuals(m1))


  ID_model <- paste(deparse(substitute(inferred_ByAnt)), VAR, sep = "-")

  if (has_interaction(m1) == T) {
    # if model has interaction, create a pasted variable, rerun model on it and compute posthocs
    print("significant interaction of predictors: make new var in this format 'VAR1_VAR2' ")
    inferred_ByAnt$PERIOD_TREATMENT <- with(inferred_ByAnt, PERIOD_TREATMENT <- paste(PERIOD, TREATMENT, sep = "_"))
    m1 <- lmer(sqrt(VARIABLE) ~ PERIOD_TREATMENT + (1 | REP_treat / Rec_Name), data = inferred_ByAnt)
    # ID_model
    posthoc_list <- compute_posthocs(m1)
  } else if (has_interaction(m1) == F) {
    print("No significant interaction of predictors")
    # if model has been simplified, run compute_posthocs (looks at the modle to see if there are interactions)
    posthoc_list <- compute_posthocs(m1)
  }

  if (REPORT) {
    Report_stats <- c(Report_stats, list(paste(
      deparse(substitute(inferred_ByAnt)),
      VAR,
      "Skewness:", print(round(skewness(residuals(m1)), 2)),
      "Kurtosis:", print(round(kurtosis(residuals(m1)), 2)),
      m1 %>% report(),
      sep = " - "
    )))
  }
}

# when there is interaction use type 3 anova Anova(model, type"III")

#####################################################################################
############  PLOTS   ###############################################################
#####################################################################################

### GROOMING LOCATION
# mean by ant
Groom_location <- aggregate(REP_treat ~ TREATMENT + PERIOD + ant1.zones + Rec_Name + Ant_status, FUN = length, na.action = na.pass, inferred)
colnames(Groom_location)[match("REP_treat", colnames(Groom_location))] <- "Count_byAnt"
Groom_location$ant1.zones <- str_replace(Groom_location$ant1.zones, "1", "Nest Area")
Groom_location$ant1.zones <- str_replace(Groom_location$ant1.zones, "2", "Foraging Area")
# mean by nest
Groom_location_MEAN <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + ant1.zones + Ant_status, FUN = mean, na.rm = T, na.action = na.pass, Groom_location) # ; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
Groom_location_SE <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + ant1.zones + Ant_status, FUN = std.error, na.rm = T, na.action = na.pass, Groom_location)
colnames(Groom_location_SE)[match("Count_byAnt", colnames(Groom_location_SE))] <- "SD_Count"
## JOIN
Groom_location <- plyr::join(x = Groom_location_MEAN, y = Groom_location_SE, type = "full", match = "all")
Groom_location[is.na(Groom_location$SD_Count), "SD_Count"] <- 0
### Grooming Location
Groom_location_plot <- ggplot(Groom_location, aes(x = TREATMENT, y = Count_byAnt, fill = PERIOD)) +
  # geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar(aes(x = TREATMENT, ymin = Count_byAnt - SD_Count, ymax = Count_byAnt + SD_Count), position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  # facet_wrap(~ant1.zones) +
  facet_grid(Ant_status ~ ant1.zones) +
  labs(title = "Grooming Location", subtitle = expression(paste("Mean", phantom(.) %+-% phantom(.), "SE by replicate")), y = "Mean Freq by ant") +
  STYLE
# save plot
SavePrint_plot(
  plot_obj = Groom_location_plot,
  plot_name = "Groom_location_plot",
  # font_size_factor = 4,
  dataset_name = deparse(substitute(infer_full)),
  save_dir = DATADIR
)





### GROOMING PROVIDER
# mean by ant
Grooming_Act_prop <- aggregate(REP_treat ~ TREATMENT + PERIOD + Act_Name + Ant_status + AntTask_Act, FUN = length, na.action = na.pass, inferred)
colnames(Grooming_Act_prop)[match("REP_treat", colnames(Grooming_Act_prop))] <- "Count_byAnt"

# mean by nest
Grooming_Act_prop_MEAN <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + Ant_status + AntTask_Act, FUN = mean, na.rm = T, na.action = na.pass, Grooming_Act_prop) # ; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
Grooming_Act_prop_SE <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + Ant_status + AntTask_Act, FUN = std.error, na.rm = T, na.action = na.pass, Grooming_Act_prop)
colnames(Grooming_Act_prop_SE)[match("Count_byAnt", colnames(Grooming_Act_prop_SE))] <- "SD_Count"
## JOIN
Grooming_Act_prop <- plyr::join(x = Grooming_Act_prop_MEAN, y = Grooming_Act_prop_SE, type = "full", match = "all")
Grooming_Act_prop[is.na(Grooming_Act_prop$SD_Count), "SD_Count"] <- 0
### Grooming Location
Grooming_Act_prop_plot <- ggplot(Grooming_Act_prop, aes(x = TREATMENT, y = Count_byAnt, fill = PERIOD)) +
  # geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar(aes(x = TREATMENT, ymin = Count_byAnt - SD_Count, ymax = Count_byAnt + SD_Count), position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  # facet_wrap(~ant1.zones) +
  facet_grid(Ant_status ~ AntTask_Act) +
  labs(title = "Grooming performer", subtitle = expression(paste("Mean", phantom(.) %+-% phantom(.), "SE by replicate")), y = "Mean Freq by ant") #+
# STYLE
Grooming_Act_prop_plot

############################################################
############################################################
############################################################
#####  BARPLOTS FOR PRE-POST #####


warning("implement fixed bars!")

# ggplot(infer_full, aes_string(x = "PERIOD", y = MEAN, fill = "TREATMENT")) +
#   geom_errorbar(aes_string(ymin = paste0(MEAN,"-",SE), ymax = paste0(MEAN,"+",SE)),
#                 width=.6, position=position_dodge(0.9), size=1) +
#   geom_bar(position = position_dodge(0.9), stat="identity", width=0.8)+
#   STYLE +
#   colFill_TREATMENT +
#   labs(title = "Grooming in Sham and Pathogen treated colonies",
#        subtitle = paste("Mean",  "\u00b1", "SE by replicate (full period) "),
#        y = paste("Mean", suffix,"(sec) by ant",sep = " "))


# Open the cumulative PDF file
pdf(paste0(save_dir_plots, "all_plots_", Sys.Date(), ".pdf"), onefile = TRUE, width = 7, height = 4)
# Create a list of variable suffixes
suffixes <- c("Count_byAnt", "duration", "SUM_duration")
# Loop through the variables and create plots
for (suffix in suffixes) {
  MEAN <- paste0("Mean_", suffix)
  SE <- paste0("SE_", suffix)
  print(paste0("plotting and saving ", suffix))
  #### FULL PERIOD
  ## MEAN FREQUENCY,  MEAN DURATION & TOTAL DURATION
  plot_N <- paste0(MEAN, "_plot")

  # subset the posthoc list element
  posthoc_df <- posthoc_list[[grep(paste0("-", suffix, "-"), names(posthoc_list))]]
  # assign NA treatment to comply withh ggplot requirements for plotting geom_text in the case of no interaction present.
  # this solution is not generalisable
  if (!("TREATMENT" %in% names(posthoc_df))) {
    posthoc_df$TREATMENT <- NA
  }

  assign(
    plot_N,
    ggplot(infer_full, aes_string(x = "PERIOD", y = MEAN, fill = "TREATMENT")) +
      geom_errorbar(aes_string(ymin = paste0(MEAN, "-", SE), ymax = paste0(MEAN, "+", SE), color = "Ant_status"),
        position = position_dodge2(width = 0.8, preserve = "single")
      ) +
      geom_col(position = position_dodge2(width = 0.8, preserve = "single"), aes(color = Ant_status)) +
      STYLE +
      colFill_TREATMENT +
      scale_color_manual(values = c("black", "red")) +
      labs(
        title = "Grooming in Sham and Pathogen treated colonies",
        subtitle = paste("Mean", "\u00b1", "SE by replicate (full period) "),
        y = paste("Mean", suffix, "(sec) by ant", sep = " ")
      ) +
      if ("PERIOD_TREATMENT" %in% names(posthoc_df)) {
        geom_text(data = posthoc_df, aes(x = PERIOD, y = 70, group = TREATMENT, label = letters, fontface = "bold"), position = position_dodge2(width = 0.9, preserve = "single"))
      } else {
        geom_text(data = posthoc_df, aes(x = PERIOD, y = 70, group = 1, label = letters, fontface = "bold"), position = position_dodge2(width = 0.9, preserve = "single"))
      }
  )

  # save plot
  SavePrint_plot(
    plot_obj = get(plot_N),
    plot_name = plot_N,
    dataset_name = deparse(substitute(infer_full)),
    save_dir = DATADIR
  )

  #### 4H BINS | barplot
  if (exists(deparse(substitute(infer_bin_h4_trim)))) {
    plot_N_4H <- paste0(MEAN, "_plot")

    assign(
      plot_N_4H,
      ggplot(infer_bin_h4_trim, aes_string(x = "PERIOD", y = MEAN, fill = "TREATMENT")) +
        geom_errorbar(aes_string(ymin = paste0(MEAN, "-", SE), ymax = paste0(MEAN, "+", SE)),
          position = position_dodge2(width = 0.8, preserve = "single")
        ) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        STYLE +
        colFill_TREATMENT +
        labs(
          title = "Grooming in Sham and Pathogen treated colonies",
          subtitle = paste("Mean", "\u00b1", "SE by replicate in 4h BLOCK (since 12:00)"),
          y = paste("Mean", suffix, "(sec) by ant", sep = " ")
        )
    )
    # save plot
    SavePrint_plot(
      plot_obj = get(plot_N_4H),
      plot_name = plot_N_4H,
      dataset_name = deparse(substitute(infer_bin_h4_trim)),
      save_dir = DATADIR
    )
  }

  #### 1H BINS | barplot
  if (exists(deparse(substitute(infer_bin_h4_trim)))) {
    plot_N_1H <- paste0(MEAN, "_plot")

    assign(
      plot_N_1H,
      ggplot(infer_bin_1h_trim, aes_string(x = "PERIOD", y = MEAN, fill = "TREATMENT")) +
        geom_errorbar(aes_string(ymin = paste0(MEAN, "-", SE), ymax = paste0(MEAN, "+", SE)),
          position = position_dodge2(width = 0.8, preserve = "single")
        ) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        facet_wrap(~time_of_day) +
        STYLE +
        colFill_TREATMENT +
        labs(
          title = "Grooming in Sham and Pathogen treated colonies",
          subtitle = paste("Mean", "\u00b1", "SE by replicate per HOUR"),
          y = paste("Mean", suffix, "(sec) by ant", sep = " ")
        )
    )
    # save plot
    SavePrint_plot(
      plot_obj = get(plot_N_1H),
      plot_name = plot_N_1H,
      dataset_name = deparse(substitute(infer_bin_1h_trim)),
      save_dir = DATADIR
    )
  }


  #### DELTA POST PRE | barplot
  if (exists(deparse(substitute(infer_full_DELTA)))) {
    plot_N_DELTA <- paste0(MEAN, "_plot")

    assign(
      plot_N_DELTA,
      ggplot(infer_full_DELTA, aes_string(x = "TREATMENT", y = MEAN, fill = "TREATMENT")) +
        geom_errorbar(aes_string(ymin = paste0(MEAN, "-", SE), ymax = paste0(MEAN, "+", SE)),
          position = position_dodge2(width = 0.8, preserve = "single")
        ) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        STYLE +
        colFill_TREATMENT +
        labs(
          title = "Grooming in Sham and Pathogen treated colonies",
          subtitle = paste("Mean", "\u00b1", "pooled SE by replicate (full period)"),
          y = paste("\u0394", "mean", suffix, "\nPOST treatment relative to PRE")
        )
    )

    # save plot
    SavePrint_plot(
      plot_obj = get(plot_N_DELTA),
      plot_name = plot_N_DELTA,
      dataset_name = deparse(substitute(infer_full_DELTA)),
      save_dir = DATADIR
    )
  }

  ##### 3H BINS | scatterplot + line
  if (exists(deparse(substitute(infer_bin_3h)))) {
    plot_N_3H <- paste0(MEAN, "_plot")

    assign(
      plot_N_3H, ggplot(infer_bin_3h, aes_string(x = "timespan", y = MEAN, group = "TREATMENT", color = "TREATMENT")) +
        geom_vline(xintercept = 0, color = "red") +
        geom_jitter(width = 0.025) +
        # geom_point(size=1) +
        geom_smooth(data = subset(infer_bin_3h, PERIOD == "pre"), method = "lm") + # , formula = y ~ x + I(x^2)
        geom_smooth(data = subset(infer_bin_3h, PERIOD == "post"), method = "lm") + # , formula = y ~ x + I(x^2)
        STYLE_CONT +
        colFill_TREATMENT +
        colScale_TREATMENT +
        labs(
          title = "Grooming in Sham and Pathogen treated colonies",
          subtitle = paste("Mean", "\u00b1", "SE by replicate per 3 HOURS"),
          x = "time from treatment", y = paste("Mean", suffix, "(sec) by ant", sep = " ")
        ) #+
      # facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)
    )

    SavePrint_plot(
      plot_obj = get(plot_N_3H),
      plot_name = plot_N_3H,
      dataset_name = deparse(substitute(infer_bin_3h)),
      save_dir = DATADIR
    )
  }
}
# Close the cumulative PDF file
dev.off()







############################################################
#### TO THE LEVEL OF ANT
#####  SCATTER PLOTS FOR 4H BINS ##

## MEAN FREQUENCY

ggplot(inferred_bin_h4_ByAnt_trim, aes(x = PERIOD, y = Count_byAnt, colour = TREATMENT)) +
  geom_point(aes(colour = TREATMENT, alpha = 0.5, stroke = 0), position = position_jitterdodge()) +
  geom_violin(position = position_dodge(width = .75), fill = NA) + # make the outline of violin black
  stat_summary((aes(group = TREATMENT)), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, colour = "red", position = position_dodge(width = .75)) +
  STYLE +
  colFill_TREATMENT +
  colScale_TREATMENT +
  labs(title = "Frequency of Grooming in Sham and Pathogen treated colonies", subtitle = expression(paste("Individual datapoints with median (red) in 4h BLOCK (since 12:00)")), y = "Mean Freq by ant")


## MEAN DURATION
ggplot(inferred_bin_h4_ByAnt_trim, aes(x = PERIOD, y = duration, colour = TREATMENT)) +
  geom_point(aes(colour = TREATMENT, alpha = 0.5, stroke = 0), position = position_jitterdodge()) +
  geom_violin(position = position_dodge(width = .75), fill = NA) + # make the outline of violin black
  stat_summary((aes(group = TREATMENT)), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, colour = "red", position = position_dodge(width = .75)) +
  STYLE +
  colFill_TREATMENT +
  colScale_TREATMENT +
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies", subtitle = expression(paste("Individual datapoints with median (red) in 4h BLOCK (since 12:00)")), y = "Mean duration by ant")


## TOTAL DURATION
ggplot(inferred_bin_h4_ByAnt_trim, aes(x = PERIOD, y = SUM_duration, colour = TREATMENT)) +
  geom_point(aes(colour = TREATMENT, alpha = 0.5, stroke = 0), position = position_jitterdodge()) +
  geom_violin(position = position_dodge(width = .75), fill = NA) + # make the outline of violin black
  stat_summary((aes(group = TREATMENT)), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, colour = "red", position = position_dodge(width = .75)) +
  STYLE +
  colFill_TREATMENT +
  colScale_TREATMENT +
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies", subtitle = expression(paste("Individual datapoints with median (red) in 4h BLOCK (since 12:00)")), y = "Total duration by ant")

############################################################

### COUNT N OF EXPOSED RECEIVERS
Reps_N_exposed$N_received <- NA
for (REP.TREAT in unique(Reps_N_exposed$REP_treat)) {
  AntList <- as.numeric(unlist(strsplit(Reps_N_exposed[which(Reps_N_exposed$REP_treat == REP.TREAT), "N_ants"], ",")))
  GroomingRec <- inferred[which(inferred$REP_treat == REP.TREAT), "Rec_Name"]
  GroomingRec <- unique(gsub("ant_", "", GroomingRec))
  # How many exposed ants received grooming
  Reps_N_exposed[which(Reps_N_exposed$REP_treat == REP.TREAT), "N_received"] <- length(intersect(unique(GroomingRec), AntList))
}


### N of exposed ants in the colony
Reps_N_exposed$treat <- str_sub(Reps_N_exposed$REP_treat, -2, -1)
Mean_ants_exp <- aggregate(N_received ~ treat, FUN = mean, na.rm = T, na.action = na.pass, Reps_N_exposed)
SD_ants_exp <- aggregate(N_received ~ treat, FUN = sd, na.rm = T, na.action = na.pass, Reps_N_exposed)
colnames(SD_ants_exp)[match("N_received", colnames(SD_ants_exp))] <- "SD_received"
N.ants.exposed <- plyr::join(x = Mean_ants_exp, y = SD_ants_exp, type = "full", match = "all")
data.frame(treat = N.ants.exposed$treat, N_exposed = sprintf("%.2f \U00B1 %.2f", N.ants.exposed$N_received, N.ants.exposed$SD_received))





##########################################################################################
### EXTRA UNRELATED PLOTS FOR PRESENTATIONS ##############################################

if (Presentation_plots) {

  #######################################################################
  ###### PLOTTING A SINGLE GROOMING INTERACTION #########################

  # Load the relevant libraries - do this every time
  library(plyr)
  library(lubridate)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library("scales")
  library("viridis")
  library(RSBID) # samp;ing

  Grooming_int <- read.table("/home/cf19810/Documents/presentation_plot_data/Act28_Rec29_frames12357-12586-ROW15_uniqueID115_Grooming", header = T, stringsAsFactors = F, sep = ",")
  Grooming_int <- Grooming_int %>% dplyr::select(-contains(c(".x", ".y", ".angle", ".zone", ".time_sec", "tim-interval", "dt_FRAME", "time_interval")))

  for (variable in names(Grooming_int)) {
    if (variable != "frame") {
      Grooming_int[variable] <- as.numeric(rescale(Grooming_int[, variable]))
    }
  }

  Grooming_int <- Grooming_int[180:229, ]
  Grooming_int$frame <- 1:length(Grooming_int$frame)

  # long format
  Grooming_int_long <- melt(setDT(Grooming_int), id.vars = c("frame"), variable.name = "Vars")

  ## heatmap
  ggplot(Grooming_int_long, aes(frame, Vars)) +
    geom_tile(aes(fill = value), colour = "white", na.rm = T) +
    scale_fill_viridis(option = "B") +
    # guides(fill=guide_legend(title="Total Incidents")) +
    theme_bw() +
    theme_minimal() +
    labs(
      title = "",
      x = "frame", y = "variable"
    ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", plot.title.position = "plot")

  #################################################
  #### PLOTTING ALL HIT MISS DATA FOR A SUBSET ####

  Train_HitMiss <- read.table("/home/cf19810/Documents/presentation_plot_data/Training_R9SP_HitMiss_Auto_Grooming.txt", header = T, stringsAsFactors = F, sep = ",")
  Train_HitMiss <- Train_HitMiss %>% select(-contains(c("REPLICATE", "PERIOD", "BEH", "ROW", "Name", "unique_interaction_id", "ant", "frame", "duration_sec", "pair", "disagreement")))

  ### get labels of Relief selected var
  RELIEF_selected <- c("prop_time_undetected_REC", "sum_moved_distance_px_ACT", "prop_time_undetected_ACT", "transfmean_speed_pxpersec_REC", "transfmean_abs_Body_Rotation_REC", "transfSD_abs_jerk_PxPerSec3_REC", "stDev_turnAngle_REC", "transfmean_abs_accel_pxpersec2_REC", "transfmean_abs_jerk_PxPerSec3_REC", "stDev_Body_Rotation_REC", "transfmean_abs_ang_Velocity_Body_REC", "stDev_ang_Velocity_Body_ACT", "transfmean_abs_ang_Velocity_Body_ACT", "transfmean_abs_Body_Rotation_ACT", "transfSD_abs_Body_Rotation_REC", "stDev_Body_Rotation_ACT", "transfSD_abs_ang_Velocity_Body_ACT", "chull_area_REC", "sum_moved_distance_px_REC", "transfSD_abs_Body_Rotation_ACT", "transfmean_speed_pxpersec_ACT", "chull_area_ACT", "transfSD_abs_ang_Velocity_Movement_REC", "transfmean_abs_accel_pxpersec2_ACT", "stDev_turnAngle_ACT", "root_mean_square_deviation_px_ACT", "root_mean_square_deviation_px_REC")

  # scaling for plotting
  for (variable in names(Train_HitMiss)) {
    if (variable != "Hit") {
      Train_HitMiss[variable] <- as.numeric(rescale(Train_HitMiss[, variable]))
    }
  }

  Train_HitMiss <- Train_HitMiss[complete.cases(Train_HitMiss), ]
  Train_HitMiss$Hit <- as.factor(Train_HitMiss$Hit)
  Train_HitMiss_sampled <- SBC(Train_HitMiss, "Hit") # Under-Sampling Based on Clustering (SBC)

  # add colour highlight
  a <- ifelse(names(Train_HitMiss) %in% RELIEF_selected, "black", "gray")
  Alpha <- ifelse(names(Train_HitMiss) %in% RELIEF_selected, "1", "0")

  # Train_HitMiss_sampled$unique_interaction_id <- 1:length(Train_HitMiss_sampled$unique_interaction_id)

  # assign numeration by hit
  Train_HitMiss_sampled <- Train_HitMiss_sampled %>%
    group_by(Hit) %>%
    dplyr::mutate(id = row_number())

  # long format
  Train_HitMiss_long <- melt(setDT(Train_HitMiss_sampled), id.vars = c("id", "Hit"), variable.name = "Vars")

  Alpha <- ifelse(Train_HitMiss_long$Vars %in% RELIEF_selected, "1", "0")
  levels(Train_HitMiss_long$Hit) <- list("Non Grooming" = "0", "Grooming" = "1")

  ## heatmap
  ggplot(Train_HitMiss_long, aes(id, Vars)) +
    geom_tile(aes(fill = value), colour = "white", na.rm = T) +
    scale_fill_viridis() +
    facet_wrap(~Hit) +
    # guides(fill=guide_legend(title="Total Incidents")) +
    theme_bw() +
    theme_minimal() +
    labs(
      title = "Heatmap of movement variables during a Grooming interaction",
      x = "interaction N", y = ""
    ) +
    theme(
      strip.text.x = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title.position = "plot"
    )

  ## heatmap2
  ggplot(Train_HitMiss_long, aes(id, Vars)) +
    geom_tile(aes(fill = value, alpha = Alpha), colour = "white", na.rm = T) +
    scale_fill_viridis() +
    facet_wrap(~Hit) +
    # guides(fill=guide_legend(title="Total Incidents")) +
    theme_bw() +
    theme_minimal() +
    labs(
      title = "Heatmap of movement variables during a Grooming interaction",
      x = "interaction N", y = ""
    ) +
    theme(
      strip.text.x = element_text(size = 18),
      axis.text.y = element_text(colour = a),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title.position = "plot"
    )

  +
    scale_x_discrete(labels = c("0" = "Non Grooming", "1" = "Grooming"))
}



####################################################################################
################ SCRAPS ############################################################

# ## COUNTS
# Xpos <- barplot( Count ~ period , Counts_AUTO_MEAN, beside=T, xlab="", ylab=" ", ylim=c(0,30)
#                  ,main="Auto classified")
#
# ##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
# segments(x0 = Xpos[1,],
#          x1 = Xpos[1,],
#          y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"],
#          y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="pre"],
#          lwd=2)
#
# segments(x0 = Xpos[2,],
#          x1 = Xpos[2,],
#          y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"],
#          y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"],
#          lwd=2)
# #
# # text(x = ((Xpos[1,]+Xpos[2,])/2),
# #      y = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"]+15,
# #      stars.pval(posthoc_FREQ_summary$p.value))
# mtext("comparison of detected grooming", line=0, side=3, outer=T, cex=1.5)
#




# inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
# inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
#
# #split REPS
#
# inferred_R3SP <- inferred[which(inferred$REPLICATE=="R3SP"),]
# inferred_R9SP <- inferred[which(inferred$REPLICATE=="R9SP"),]
#
# #bins of an hour
# inf_R3SP_bin <- table(cut(inferred_R3SP$T_start_UNIX, breaks="hour"))
# inf_R9SP_bin <- table(cut(inferred_R9SP$T_start_UNIX, breaks="hour"))
#
#
# barplot(inf_R9SP_bin)

#
# end <- inferred_R3SP[which(inferred_R3SP$T_stop_UNIX==max(inferred_R3SP$T_stop_UNIX)),"T_stop_sec"]
#
# myFrame <- as.data.frame(table(myTable))
#
# ggplot(inf_R3SP_bin, aes(, Y)) + geom_point() + geom_vline(xintercept = as.Date("2020-07-01"))
#
#
#
# #
# barplot(inf_R3SP_bin) + abline(v=end)
#
# +  abline(v =  ymd_hms(max(inferred_R3SP$T_start_UNIX)-4*3600, tz="GMT"))
#
#
# + abline(v = as.POSIXct(strptime("2021-03-17 07:56:42 GMT", format="%Y-%m-%d %H:%M:%OS")))
#
# + abline(v =  max(inferred_R3SP$T_start_UNIX)-24*3600)
#
#
# exp_R3SP_time <-  as.POSIXct( "2021-03-15 12:11:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
# #exp_R3SP_time <-max(inferred_R3SP$T_start_UNIX)-24*3600 #more formally correct
# exp_R9SP_time <-  as.POSIXct( "2021-04-26 11:26:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
# #exp_R9SP_time <-max(inferred_R9SP$T_start_UNIX)-24*3600 #more formally correct

# inferred_R3SP$PERIOD_new <- NA
# inferred_R3SP$PERIOD_new <- ifelse(inferred_R3SP$T_start_UNIX < exp_R3SP_time, "pre", "post")
# #n occurrences
#

#
# ##############
# # Load manual annotations
#
# annotations <- read.csv(paste(DATADIR,"/annotations_TRAINING_DATASET.csv",sep = ""), sep = ",")
# #transform zulu time in GMT
# annotations$T_start_UNIX <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
# annotations$T_stop_UNIX  <- as.POSIXct(annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
# #assign time in sec to avoid issues on time management and matching
#
# #SELECT ONLY the exposed nurses
# annotations_R3 <-  annotations[which(annotations$treatment_rep=="R3SP" & annotations$Receiver %in% c(5,17)),]
# annotations_R9 <-  annotations[which(annotations$treatment_rep=="R9SP" & annotations$Receiver %in% c(23,29,32)),]
# annotations <- rbind(annotations_R3,annotations_R9)
#
# ## count the number of observations of each behaviour - WARNING; some behavs not observed e.g. before the treatment (period), so will need to account for that (next step)
# Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + period + treatment_rep, FUN=length, na.action=na.pass, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
# ## calculate mean durations for each behaviour
# Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + period + treatment_rep, FUN=mean, na.rm=T, na.action=na.pass, annotations)
# ## merge counts & durations carefully
# Counts_by_Behaviour_CLEAN    <-  plyr::join(x=Counts_by_Behaviour_CLEAN, y=Durations_by_Behaviour_CLEAN, type = "full", match = "all")
#
# ## create a data frame with all combinations of the conditioning variables
# all_combos <- expand.grid ( Behaviour=unique(annotations$Behaviour), period=unique(annotations$period), treatment_rep=unique(annotations$treatment_rep))
#
# ## add the missing cases
# Counts_by_Behaviour_AllCombos <- plyr::join (x = Counts_by_Behaviour_CLEAN , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )
#
# ## Focus only on a few important behaviours
# Counts_by_Behaviour_AllCombos$Behaviour <- as.character(Counts_by_Behaviour_AllCombos$Behaviour)  ## naughty R
# Counts_by_Behaviour_AllCombos <- Counts_by_Behaviour_AllCombos[which(Counts_by_Behaviour_AllCombos$Behaviour %in% c("G")),]
#
# ## replace the NAs with 0 counts
# Counts_by_Behaviour_AllCombos$Count[which(is.na(Counts_by_Behaviour_AllCombos$Count))] <- 0
# ## finally, get the mean & S.E. for each behav before/after  for barplots
# Counts_by_Behaviour_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=na.pass, Counts_by_Behaviour_AllCombos)
# Counts_by_Behaviour_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=std.error, na.rm=T, na.action=na.pass, Counts_by_Behaviour_AllCombos)
#
#
#
# ## show the mean counts for each behav | stage
# pdf(file=paste(DATADIR,"Grooming_Auto_Man__pre-post.pdf", sep = ""), width=5, height=8)
# par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.5,0.1), mgp=c(1.3,0.3,0), tcl=-0.2,oma=c(0,0,2,0))
#
# ## COUNTS
# Counts_by_Behaviour_MEAN$period <- factor(Counts_by_Behaviour_MEAN$period , levels = c("pre","post"))
#
# Xpos <- barplot( Count ~ period , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour count" , ylim=c(0,30)
#                  ,main="Manual annotation")
# ##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
# segments(x0 = Xpos[1,],
#          x1 = Xpos[1,],
#          y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"],
#          y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="pre"],
#          lwd=2)
#
# segments(x0 = Xpos[2,],
#          x1 = Xpos[2,],
#          y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"],
#          y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"],
#          lwd=2)
# #
# # text(x = ((Xpos[1,]+Xpos[2,])/2),
# #      y = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"]+15,
# #      stars.pval(posthoc_FREQ_summary$p.value))
#
#
#
# #########





# ###################################################################
# ####################### GAMM FITTING #############################
#
# inferred_bin_1h_ByRep$REP_treat <- as.factor(inferred_bin_1h_ByRep$REP_treat)
# inferred_bin_10min_ByRep$REP_treat <- as.factor(inferred_bin_10min_ByRep$REP_treat)
#
#
# #### FIT GAMMM ########
# #GAM can capture complex relationships by fitting a non-linear smooth function through the data, while controlling how wiggly the smooth can get
#
# ## should I use discrete time windows?                        currently: 1h bins
# ## should I do it on means by ant or with individual ants?    currently: means by ant
# ## given that the pre-post effect is given, I did subset for post only
#
# ## some guidance:
# ## https://stats.stackexchange.com/questions/403772/different-ways-of-modelling-interactions-between-continuous-and-categorical-pred
# ## http://r.qcbs.ca/Workshops/workshop08/workshop08-en/workshop08-en.html#41
#
# ############################# 1h blocks
# ## test duration / post exposure
# inferred_bin_1h_ByRep_post <- inferred_bin_1h_ByRep[which(inferred_bin_1h_ByRep$PERIOD=="post"),]
# hist(log10(inferred_bin_1h_ByRep_post$duration+0.0001), breaks =100)
#
# VAR <- "duration"
#
# #transform the var before to avoid issues with model terms extraction
# # m1 models the observations as coming from either a smooth timespan effect depending on which TREATMENT the observation comes from (the TREATMENT parametric term just sets the mean density for each species and is needed as discussed above), plus a random intercept.
# # Taken together, the curves for individual colonies REP arise from shifted versions of the TREATMENT specific curves, with the amount of shift given by the random intercept. This model assumes that all colonies REP have the same shape of smooth as given by the smooth for the particular TREATMENT that colonies REP comes from.
# # k is the number of basis functions, aka the maximum possible degrees of freedom allowed for a smooth term in the model (set as high as possible)
# #  s(x1, by = x2) is the interaction term for smoothed var * categorical var
# m1 <- mgcv::gam(Box_Cox(inferred_bin_1h_ByRep_post[,VAR]) ~ TREATMENT + s(timespan,m = 1, by = TREATMENT, k = 25) + s(REP_treat, bs = 're'),
#                 data =  inferred_bin_1h_ByRep_post, method = 'REML')
# # version without interaction | can't be trusted, the plot fitted values shows altered relationships
# m2 <- mgcv::gam(sqrt(duration) ~ TREATMENT +  s(timespan) + s(REP_treat, bs = 're'),
#                 data =  inferred_bin_1h_ByRep_post, method = 'REML')
# #version with an added random effect | can't be trusted, the plot fitted values shows altered relationships
# #global or average smooth effect of x on y (the s(x) term) plus a smooth difference term (the second s(x, by = f, m = 1) term). As the penalty here is on the first derivative (m = 1) for this difference smoother, it is penalising departure from a flat line, which when added to the global or average smooth term (s(x)) reflects a deviation from the global or average effect.
# m3 <- gam(sqrt(duration) ~ TREATMENT + s(timespan) + s(timespan, by = TREATMENT, m = 1, k = 25) + s(REP_treat, bs = 're'),
#           data = inferred_bin_1h_ByRep_post, method = 'REML')
#
#
#
# # Plot the smooth terms: understand the effect of the timespan variable on the response.
# plot(m1, pages = 1, shade = TRUE, shade.col = "lightblue")
# # Plot residuals vs. fitted values: the residuals should be randomly scattered around zero
# plot(fitted(m1),  residuals(m1), xlab = "Fitted values", ylab = "Residuals")
# abline(h = 0, col = "red")
# #QQ-plot:normality of the residuals.
# qqnorm(resid(m1))
# qqline(resid(m1))
# # Skewness-kurtosis plot
# descdist(resid(m1))
# # estimated density
# # Plot the estimated density
# plot(density(inferred_bin_1h_ByRep_post$duration), main = "Estimated Density", xlab = "Duration", ylab = "Density")
#
# #plot fitted values and their confidence intervals
# plotme<-ggemmeans(m1,terms=c("timespan","TREATMENT"),rg.limit=26000)
# plot(plotme,colors= myColors_Treatment) + labs(title="predicted values, period POST, 1h aggregation", x="timespan (hours)",y=paste("transformed",VAR))
# summary(m1)
#
#
#
# #########################################################
# ################# 10min blocks ##########################
#
# ## test duration / post exposure
# inferred_bin_10min_ByRep_post <- inferred_bin_10min_ByRep[which(inferred_bin_10min_ByRep$PERIOD=="post"),]
# hist(Box_Cox(inferred_bin_10min_ByRep_post$duration+0.0001), breaks =100)
#
# #mod_VAR_list <- list()
#
# for (VAR in colnames(inferred_bin_10min_ByRep_post)){
#   if (VAR %in% numeric_variable_list) {
#
#     print( ggplot(inferred_bin_10min_ByRep_post, aes(x = (inferred_bin_10min_ByRep_post[,VAR]))) +
#              geom_histogram(position = "identity", bins = 30) + facet_wrap(~TREATMENT)  +
#              xlab(VAR) )
#
#     #Counts_by_Behaviour_AllCombos1$period = relevel(Counts_by_Behaviour_AllCombos1$period, ref="pre")
#     print(paste(VAR,"######################"),sep=" ")
#
#
#
#
#     # m1 models the observations as coming from either a smooth timespan effect depending on which TREATMENT the observation comes from (the TREATMENT parametric term just sets the mean density for each species and is needed as discussed above), plus a random intercept.
#     # Taken together, the curves for individual colonies REP arise from shifted versions of the TREATMENT specific curves, with the amount of shift given by the random intercept. This model assumes that all colonies REP have the same shape of smooth as given by the smooth for the particular TREATMENT that colonies REP comes from.
#     # k is the number of basis functions, aka the maximum possible degrees of freedom allowed for a smooth term in the model (set as high as possible)
#     #  s(x1, by = x2) is the interaction term for smoothed var * categorical var
#     # m = 1 penalising departure from straight line
#     m1 <- mgcv::gam(Box_Cox(inferred_bin_10min_ByRep_post[,VAR]) ~ TREATMENT + s(timespan,m = 1, by = TREATMENT, k = 30) + s(REP_treat, bs = 're'),
#                     data =  inferred_bin_10min_ByRep_post, method = 'REML')
#     # # m = 2 penalising departure from a straight line with constant slope
#     # m2 <- mgcv::gam(Box_Cox(VAR) ~ TREATMENT + s(timespan, m = 2, by = TREATMENT, k = 30) + s(REP_treat, bs = 're'),
#     #                 data =  inferred_bin_10min_ByRep_post, method = 'REML')
#     # #version with an added random effect | can it be trusted?
#     # #global or average smooth effect of x on y (the s(x) term) plus a smooth difference term (the second s(x, by = f, m = 1) term). As the penalty here is on the first derivative (m = 1) for this difference smoother, it is penalising departure from a flat line, which when added to the global or average smooth term (s(x)) reflects a deviation from the global or average effect.
#     # m3 <- gam(Box_Cox(VAR) ~ TREATMENT + s(timespan) + s(timespan, by = TREATMENT, k = 30) + s(REP_treat, bs = 're'),
#     #           data = inferred_bin_10min_ByRep_post, method = 'REML')
#     #
#     # AIC(m1,m2,m3)
#
#     # Plot the smooth terms: understand the effect of the timespan variable on the response.
#     plot(m1, pages = 1, shade = TRUE, shade.col = "lightblue")
#     # Plot residuals vs. fitted values: the residuals should be randomly scattered around zero
#     plot(fitted(m1),  residuals(m1), xlab = "Fitted values", ylab = "Residuals")
#     abline(h = 0, col = "red")
#     #QQ-plot:normality of the residuals.
#     qqnorm(resid(m1))
#     qqline(resid(m1))
#     # Skewness-kurtosis plot
#     descdist(resid(m1))
#     # estimated density
#     # Plot the estimated density
#     plot(density(inferred_bin_10min_ByRep_post[,VAR]), main = "Estimated Density", xlab = VAR, ylab = "Density")
#
#     #plot fitted values and their confidence intervals
#     plotme<-ggemmeans(m1,terms=c("timespan","TREATMENT"),rg.limit=26000)
#     plot(plotme,colors= myColors_Treatment) + labs(title="predicted values, period POST, 10min aggregation", x="timespan (hours)",y=paste("transformed",VAR))
#     summary(m1)
#
#     #  Chi-square test to compare lm.1 and lm.2 (i.e. it tests whether reduction in the residual sum of squares are statistically significant or not). Note that this makes sense only if lm.1 and lm.2 are nested models.
#     anova(m1,m2, m3, test = "Chisq")
#     # the more complex model is not significantly better at encapsulating data variation
#
#
#
#
#
#
#
#     # Give the rows meaningful names:
#     mod_VAR_list <- c(mod_VAR_list,list(lt))
#     names(mod_VAR_list)[length(mod_VAR_list)] <- paste(summary(m1)$call$data,VAR,sep = "-")
#
#   }}
#
#
#
# # emmeans_res <- emmeans(m1, ~ timespan + TREATMENT,type="response")
# # # Plot the results
# # ggplot(as.data.frame(emmeans_res), aes(x = TREATMENT, y = response, color = timespan)) +
# #   geom_point() +
# #   geom_errorbar(aes(ymin = response - SE, ymax = response + SE), width = 0.2) +
# #   labs(title = "Emmeans plot on the original response scale") +
# #   theme_bw()
#
#
#
# # # Fit the GAMM
# # gamm_model <- gamm4((log10(rel_conc_imputed + GENE_cost)) ~ Treatment, random = ~ (1 | Colony), data = GENE_data)
# # # Extract residuals
# # gamm_residuals <- resid(gamm_model$gam, type = "pearson")
# # # Plot residuals vs. fitted values
# # plot(fitted(gamm_model$gam), gamm_residuals, xlab = "Fitted values", ylab = "Pearson residuals")
# # abline(h = 0, lty = 2, col = "red")
# # # QQ plot of residuals
# # qqnorm(gamm_residuals)
# # qqline(gamm_residuals)
# # # Summary of the model
# # summary(gamm_model$gam)
# # # Summary of the random effects
# # summary(gamm_model$lme)
# # # Check residuals for the GAMM model
# # plot(gamm_model$lme, resid(., type = "pearson"), xlab = "Fitted values", ylab = "Pearson residuals")
# # abline(h = 0, lty = 2, col = "red")
# # qqnorm(resid(gamm_model$lme, type = "pearson"))
# # qqline(resid(gamm_model$lme, type = "pearson"))



###############################################
###############################################
###############################################



# ###### VALUES FOR PRE.POST: TIME BINS #####
# ### TO THE LEVEL OF ANT
#
# inferred_bin_h4_ByAnt_trim$TREATMENT <- as.factor(inferred_bin_h4_ByAnt_trim$TREATMENT)
# inferred_bin_h4_ByAnt_trim$PERIOD <- as.factor(inferred_bin_h4_ByAnt_trim$PERIOD)
# inferred_bin_h4_ByAnt_trim$REP_treat<- as.factor(inferred_bin_h4_ByAnt_trim$REP_treat)
# inferred_bin_h4_ByAnt_trim$Rec_Name<- as.factor(inferred_bin_h4_ByAnt_trim$Rec_Name)
# inferred_bin_h4_ByAnt_trim$time_of_day<- as.factor(inferred_bin_h4_ByAnt_trim$time_of_day)
# inferred_bin_h4_ByAnt_trim$timespan<- as.factor(inferred_bin_h4_ByAnt_trim$timespan)
#
# # Select variables for analysis
# numeric_variable_list1  <- names(inferred_bin_h4_ByAnt_trim) [ which (!names(inferred_bin_h4_ByAnt_trim) %in% c( "PERIOD","TREATMENT","REP_treat","Rec_Name" ) )]
#
# for (VAR in colnames(inferred_bin_h4_ByAnt_trim)){
#   if (VAR %in% numeric_variable_list1) {
#
#     print(ggplot(inferred_bin_h4_ByAnt_trim, aes(x = sqrt(inferred_bin_h4_ByAnt_trim[,VAR]))) +
#       geom_histogram(position = "identity", bins = 30) + facet_wrap(~TREATMENT) +
#       xlab(VAR) )
#
#     #Counts_by_Behaviour_AllCombos1$period = relevel(Counts_by_Behaviour_AllCombos1$period, ref="pre")
#     print(paste("###################### LMER OF",VAR,"######################"),sep=" ")
#     m1 <- lmer(sqrt(inferred_bin_h4_ByAnt_trim[,VAR]) ~ PERIOD * TREATMENT + (1|REP_treat/Rec_Name) , data = inferred_bin_h4_ByAnt_trim) # the "/" is for the nesting #  + (1|time_of_day)
#
#     print(paste("###################### TEST NORMALITY OF",VAR," RESIDUALS ######################"),sep=" ")
#     print(summary(m1))
#     test_norm(residuals(m1)) #test residuals' normality. null hypothesis for the Shapiro-Wilk test is that a variable is normally distributed
#
#     # POST-HOCs for TREATMENT pre-post
#     print(paste("###################### POST-HOCs for TREATMENT pre-post of ",VAR,"######################"),sep=" ")
#     posthoc_Treatment <- emmeans(m1, specs = trt.vs.ctrlk ~ PERIOD | TREATMENT) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2”
#     posthoc_Treatment_summary<-summary(posthoc_Treatment$contrasts)
#     print(posthoc_Treatment_summary)
#
#     # POST-HOCs for PERIOD by treatment
#     print(paste("###################### POST-HOCs for PERIOD by treatment of ",VAR,"######################"),sep=" ")
#     posthoc_Period <- emmeans(m1, specs = pairwise ~ TREATMENT | PERIOD) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2”
#     posthoc_Period_summary<-summary(posthoc_Period$contrasts)
#     print(posthoc_Period_summary)
#
#     par(mfrow=c(1,2))
#     plot(m1)
#     qqnorm(residuals(m1))
#     qqline(residuals(m1))
#     hist(residuals(m1))
#     #anova(m1)
#   }}





# ############# PLOTS STYLE ################
# STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),theme_bw(),
#               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Montserrat") ),
#               scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
# )
#
# STYLE_continous <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
#                         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Montserrat")),
#                         theme_bw()
# )
#
#
# STYLE_NOVIR <- list(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Montserrat")),
#                     theme_bw(),
#                     scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
# )


# # replicate folder
# for (REP.n in 1:length(files_list)) {
#   # REP.n <- 1    #temp
#   REP.folder      <- files_list[REP.n]
#   REP.files       <- list.files(REP.folder, pattern = "CapDef3.myr")
#   REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
#
#   #replicate file
#   for (REP.FILES in REP.filefolder) {
#     #get substring in variable until R9BS_
#     # REP.FILES <- REP.filefolder[1]
#     REP_treat <- sub("\\_.*", "", basename(REP.FILES))
#     #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
#
#     # REP.FILES <-  REP.filefolder[1]   #temp
#     print(REP.FILES) ##}}
#     #open experiment
#     exp <- fmExperimentOpen(REP.FILES)
#     # exp.Ants <- exp$ants
#     exp_end <- fmQueryGetDataInformations(exp)$end
#
#     ########## GET EXPOSED ANTS # AW 17June2022
#     e.Ants <- exp$ants
#     Exposed_list <- vector()
#     if (REP_treat %in% Reps_N_exposed$REP_treat) {
#     for (ant in e.Ants){
#       #exclude dead ants
#       if (FALSE %in% ant$getValues("IsAlive")[,"values"] ) { #if FALSE doesn't appear
#         # assign DEAD flag in inferred
#         inferred[which(inferred$REP_treat==REP_treat & inferred$Rec_Name==paste0("ant_",ant$ID)),"dead"] <- "dead"
#       }#exclude dead
#
#       if (TRUE %in% ant$getValues("Exposed")[,"values"]) {
#         exposed <-ant$ID
#         Exposed_list <- c(Exposed_list, exposed)
#         # assign EXPOSED flag in inferred
#         inferred[which(inferred$REP_treat==REP_treat & inferred$Rec_Name==paste0("ant_",ant$ID)),"exposed"] <- "exposed"
#       } # list exposed
#     }
#     explist<- data.frame(Exposed_list,stringsAsFactors = F)
#     #Collapse
#     explist2 <- data.frame(val=paste0(explist$Exposed_list,collapse = ', '),stringsAsFactors = F)
#
#       Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP_treat) ,"N_ants"] <- explist2
#     }
#
#     for (ROW in 1:nrow(inferred)) {
#       if (inferred[ROW,"REP_treat"]==REP_treat) {
#         #add end time minus
#         inferred[ROW,"end_time"] <-  as.POSIXct(exp_end ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
#         # inferred[ROW,"return_time"] <-  as.POSIXct( exp_end  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
#         #assign a new common time for analyses
#         #inferred[ROW,"time_since_start"] <-
#       }
#     }
#     rm(list=(c("exp"))) #remove experiment
#   }}
#
# #subtract end window
# inferred$end_time <- inferred$end_time - window_shift
# ## force unix_time
# inferred$end_time <- as.POSIXct(inferred$end_time,  origin="1970-01-01", tz="GMT" )
# inferred$return_time <- inferred$end_time - 24*3600
#
# ## WINDOW SHIFT GETS US THE APPROX TIME, ASSIGN THE FINAL RETURN TIME
# ##substitute time with hand-collected return time for increased precision
# #paste the day of inferred$return_time with the time of ReturnTime_mins
# inferred$Date_only <- as.character(as.Date(inferred$return_time, format = "%Y-%m-%d"))
# inferred <- left_join(inferred, ReturnTime_mins, by = "REP_treat")
#
#
#
# inferred$return_time  <- paste(inferred$Date_only,inferred$Return.hour,"GMT")
# inferred$return_time  <- as.POSIXct(inferred$return_time,  origin="1970-01-01", tz="GMT" )
# #the hand-collected return time is the time at which the return operation is complete. As sometimes it could have taken a few minutes, subtract 5 mins to avoid loosing grooming events
# inferred$return_time  <- inferred$return_time - (60*10) #mins
#
# ## add a common time to all rows
# # Negative before exposure, positive after
# inferred$time_stop_since_treat <- as.numeric(difftime(inferred$T_stop_UNIX,inferred$return_time, units = "secs"))
#
# inferred$PERIOD <- NA
# # Assign Pre and Post labels
# # The time between the nurses sampling from the nest and the return time is labeled as EXPOSURE_GAP (3h window)
# for (ROW in 1:nrow(inferred)) {
#   #since these are times in seconds, it is ok to  have strict ">" instead of ">="
# if(inferred[ROW,"time_stop_since_treat"] > 0){ inferred[ROW,"PERIOD"] <- "post"
# }else if ( inferred[ROW,"time_stop_since_treat"] < (-3*3600) & inferred[ROW,"time_stop_since_treat"] > (-27*3600)) {  inferred[ROW,"PERIOD"] <- "pre"  } else{  inferred[ROW,"PERIOD"] <- "EXPOSURE_GAP"}
# }
