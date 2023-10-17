#Test size differences

USER <- "supercompAdriano"

if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  #SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}
###### PACKAGES
library(lme4)
library(report)

###### SOURCES
source(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/FUNCTIONS_Analysis_and_Styling.R",sep="/"))

###### LOAD MEASURES
measures <- read.table(file.path(DATADIR,"ant_length_per_Colony.txt"),header=T,stringsAsFactors = F, sep=",")

# extract the last two characters and create a new variable
measures$Size <- ifelse(substr(measures$REP_treat, nchar(measures$REP_treat)-1, nchar(measures$REP_treat)-1) == "B", "Big", "Small")

table(measures$Size)

# plotting
ggplot(measures, aes(x = Size, y = worker_length_mm)) +
  geom_violin(trim = TRUE,width =1, alpha = 0.6) +
  geom_boxplot(colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+  #geom_jitter(height = 0.1, width = 0.15) +
  labs(x = "Treatment Size", y = "Worker length (mm)", color = "Species",caption = "mean of 2 measures per ant") +
  STYLE

model <- lmer(worker_length_mm ~ Size + (1 | REP_treat), data = measures)

output_lmer(model)
report(model)




# Create a minimal reproducible sample

min_sample <- measures %>%
  group_by(TS, REP_treat) %>%
  slice_sample(n = 1) %>%
  ungroup()
min_sample$worker_length_mm <- round(min_sample$worker_length_mm,0)
dput(as.data.frame(min_sample[,c("worker_length_mm", "TS", "REP_treat","Size")]))
