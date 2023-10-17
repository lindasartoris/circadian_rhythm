library(ggplot2)
library(dplyr)


# Define paths
path1 <- "/media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/intermediary_analysis_steps/full_interaction_lists/PostTreatment/observed"
path2 <- "/media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/intermediary_analysis_steps/full_interaction_lists/PreTreatment/observed"

# Get a list of all files in each path
files1 <- list.files(path1, full.names = TRUE)
files2 <- list.files(path2, full.names = TRUE)

# Read each file as tab-separated values
data1 <- lapply(files1, read.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data2 <- lapply(files2, read.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# Combine data frames and select columns
combined_data <- bind_rows(data1, data2) %>%
  select(Startframe, period, REP_treat)

# Print the first few rows of the combined data frame
head(combined_data)
data1 <- NULL
data2 <- NULL


# # Combine data frames and select columns
# combined_data <- bind_rows(data1, data2) %>%
#   select(Startframe, period, REP_treat)
# 
# # Group by REP_treat and period, and add a column for max Startframe
# combined_data <- combined_data %>%
#   group_by(REP_treat, period) %>%
#   mutate(max_Startframe = max(Startframe))
# 
# ##
combined_data$max_Startframe <- 360000
combined_data$min_Startframe <- 1


combined_data_1 <-  combined_data[which(combined_data$REP_treat=="R8BP"),]



# Define the bin width
bin_width <- 1500

# Split the data by REP_treat
split_data <- split(combined_data_1, combined_data_1$REP_treat)

# # Iterate over each REP_treat and create a plot
# for (i in names(split_data)) {
#   # Calculate the frequency of Startframe in each bin
#   plot_data <- split_data[[i]] %>%
#     group_by(bin = cut(Startframe, breaks = seq(0, 36000, by = bin_width))) %>%
#     summarise(freq = n()) %>%
#     mutate(REP_treat = i)
#   
#   # Create a plot with the frequency on the y-axis and the bin on the x-axis
#   plot <- ggplot(plot_data, aes(x = bin, y = freq)) +
#     geom_col() +
#     ggtitle(paste("REP_treat:", i)) +
#     xlab("Startframe (bins of 1500)") +
#     ylab("Frequency")
#   
#   # Print the plot
#   print(plot)
# }

