rm(list=ls())

########################################################################
###  INPUT TO DEFINE BY USER############################################
####Please fill in the path to the data folder, the code folder, and the c++ executables folder, as in the exemple below
# Define a function to choose the data path
# choose_data_path <- function() {
#   # Ask user which data they would like to process
#   choice <- as.character(readline(prompt = "Which data would you like to process?\nAll interactions (1)\nGrooming interactions (2)\n"))
#   # Set data_path depending on user's choice
#   if (choice == "1") {
#     data_path <- "/media/lsartori/LS_1/circadian_rhythm_2022_experiment/Lasius-Bristol_pathogen_experiment/main_experiment"
#     print("ALL INTERACTIONS")
#   } else if (choice == "2") {
#     data_path <- "/media/lsartori/LS_1/circadian_rhythm_2022_experiment/Lasius-Bristol_pathogen_experiment/main_experiment_grooming"
#     print("GROOMING INTERACTIONS")
#   } else {
#     stop("Invalid choice entered. Please choose either '1' or '2'.")
#   }
#   return(data_path)
# }
# 
# # Call the function to choose the data path
# data_path <- choose_data_path()

data_path <- "/media/lsartori/LS_2/circadian_rhythm_2022_experiment/Lasius-Bristol_pathogen_experiment/main_experiment" ##NORMAL INTERACTIONS
# data_path <- "/media/lsartori/LS_1/circadian_rhythm_2022_experiment/Lasius-Bristol_pathogen_experiment/main_experiment_grooming" ##GROOMING INTERACTIONS


code_path  <- "/media/lsartori/LS_2/circadian_rhythm_2022_experiment/scripts/Social_Network_Plasticity_Exp_2022_LS/1_data_post_processing/source"
executables_path <- "~/executables"
FRAME_RATE <- 8 #AW & LS
########################################################################
###   END INPUT ########################################################
########################################################################

source(paste(code_path,"/libraries.R",sep=""))
source(paste(code_path,"/functions_and_parameters.R",sep=""))

to_keep <- c(ls(),"to_keep")

#####Run analysis programs ####
# source(paste(code_path,"/1_trackconverter.R",sep=""))
# clean()
# source(paste(code_path,"/2_define_deaths.R",sep=""))
# clean()
# source(paste(code_path,"/3_apply_rotation_to_datfiles.R",sep=""))
# clean()
# source(paste(code_path,"/4_retagged_ant_modifications.R",sep=""))
# clean()
# source(paste(code_path,"/5_zoneconverter_nest.R",sep=""))
# clean()
# source(paste(code_path,"/6_time_investment.R",sep=""))
# clean()
# source(paste(code_path,"/7_trajectory.R",sep=""))
# clean()
# source(paste(code_path,"/8_process_trajectory_files.R",sep=""))
# clean()
# source(paste(code_path,"/9_interaction_detection.R",sep=""))
# clean()
# source(paste(code_path,"/10_process_interaction_files.R",sep=""))
# clean()

# TO USE AW/LS:
#source(paste(code_path,"/11_randomise_interactions.R",sep="")) # LS: Done
#clean()
#source(paste(code_path,"/12_simulate_transmission.R",sep="")) # LS: Done
# clean()
#source(paste(code_path,"/13_network_analysis.R",sep="")) # LS: Done (both for duration & number of events)
# clean()
 source(paste(code_path,"/14_summarise_interactions.R",sep=""))
# clean()

# heatmaps are optional (might be interesting?) -> long run time
# source(paste(code_path,"/15_heatmaps_individual.R",sep=""))
# clean()
# source(paste(code_path,"/16_heatmaps_groups.R",sep=""))
# clean()
# source(paste(code_path,"/17_brood_location.R",sep=""))
# clean()
# source(paste(code_path,"/18_process_heatmaps.R",sep=""))
# clean()