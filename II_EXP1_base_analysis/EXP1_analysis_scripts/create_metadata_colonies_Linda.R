#script to add the time stamps to the colony metadata file
# exp_start	start of the recording (Day 1, 8 AM/PM)
# exp_end	end of the experiment (Day 2, 6:30 AM/PM)
# pre_start	start of the pre-treatment network (Day 1, 9:30 AM/PM)
# pre_end	end of the pre-treatment network (Day 1/2, 6:30 AM/PM) 
# middle_start	start of the corresponding circadian period in the middle of pre- and post-treatment (Day 1/2, 9:30 AM/PM)
# middle_end	stop of the corresponding circadian period in the middle of pre- and post-treatment (Day 2, 6:30 AM/PM)
# post_start	start of the post-treatment network (Day 2, 9:30 AM/PM)
# post_end	end of the post-treatment network (Day 2/3, 6:30 AM/PM) 

#some replicates were not started able to start before or at 8 AM/PM 
#these have to be taken care of individually:
#DF_haydon_c16 (2022-11-08T08:04:42.205220638Z)
#NF_haydon_c35 2022-11-22T20:07:34.308587587Z
#NS_haydon_c48 2022-12-06T20:34:03.213602667Z
#DS__karla_c13 2022-11-01T09:53:04.400663158Z --> pre_start and pre_end can be shifted by 23 min
#DF_karla_c45 2022-11-22T08:04:43.331678602Z
#NF_karla_c60 2022-11-29T20:01:58.962857794Z
#NS_karla_67 2022-12-06T20:32:07.169920500Z 
#DS_polyakov_46 2022-11-21T08:00:50.012365123Z
#DF_polyakov_c69 2022-11-14T08:00:37.129063959Z
#DF_prideaux_c42 2022-12-05T08:01:56.472733603Z
#DF_trojan_c20 2022-12-05T08:03:41.367948407Z

#one replicate had a later post_start time because of re-tagging: 09:40:22 - removed and frozen at 18:38 so set post times to 9:38 - 18-38

rm(list=ls())

usr <- "lsartori"
WORKDIR <- paste("/media",usr,"LS_2/circadian_rhythm_2022_experiment",sep="/")
metadata_colonies_file <- paste(WORKDIR,"scripts/EXP1_base_analysis", "metadata_colonies_for_R.csv", sep = "/")
metadata_colonies <- read.table(metadata_colonies_file, header = T, stringsAsFactors = F, sep = ",") ### body length information
metadata_colonies$exp_start <- as.POSIXct(metadata_colonies$exp_start, format = "%Y-%m-%dT%H:%M:%OSZ", origin = "1970-01-01", tz = "GMT")      
head(metadata_colonies)
str(metadata_colonies)


for (i in 1:nrow(metadata_colonies)) {
  if (metadata_colonies$treatment[i] %in%c("DF", "DS")) {
    print(paste("day ",i))
    metadata_colonies$post_start[i] <- as.POSIXct(as.POSIXct(paste(format(metadata_colonies$exp_start[i] + (25.5*3600), "%Y-%m-%d"), "09:30:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))
    # metadata$exp_start_time <- as.POSIXct(paste(format(time_treatment_start - 86400, "%Y-%m-%d"), "09:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
    
  } else {
    print(paste("night ",i))
    metadata_colonies$post_start[i] <- as.POSIXct(paste(format(metadata_colonies$exp_start[i]+ (25.5*3600), "%Y-%m-%d"), "21:30:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
    # metadata$exp_start_time <- as.POSIXct(paste(format(time_treatment_start - 86400, "%Y-%m-%d"), "09:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
    
  }
}

str(metadata_colonies)

#back from seconds to GMT
metadata_colonies$post_start <- as.POSIXct(metadata_colonies$post_start, origin = "1970-01-01", tz = "GMT")
str(metadata_colonies)
#for all others
metadata_colonies$post_end <- metadata_colonies$post_start + 9*3600       # until end of tracking 6:30 AM/PM 
metadata_colonies$exp_end <- metadata_colonies$post_end
metadata_colonies$middle_end <- metadata_colonies$post_start - 3*3600 
metadata_colonies$middle_start <- metadata_colonies$middle_end - 9*3600 
metadata_colonies$pre_end <- metadata_colonies$middle_start - 3*3600 
metadata_colonies$pre_start <- metadata_colonies$pre_end - 9*3600 

#for c13
metadata_colonies[metadata_colonies$colony=="c13", "pre_start"] <- paste(format(metadata_colonies$exp_start[metadata_colonies$colony=="c13"], "%Y-%m-%d"), "10:00:00")
metadata_colonies[metadata_colonies$colony=="c13", "pre_end"] <- metadata_colonies$pre_start[metadata_colonies$colony=="c13"] + 9*3600 

#for c52
metadata_colonies[metadata_colonies$colony=="c52", "post_start"] <- paste(format(metadata_colonies$post_start[metadata_colonies$colony=="c52"], "%Y-%m-%d"), "09:38:00")
metadata_colonies[metadata_colonies$colony=="c52", "post_end"] <- metadata_colonies$post_start[metadata_colonies$colony=="c52"] + 9*3600 
metadata_colonies[metadata_colonies$colony=="c52", "exp_end"] <- metadata_colonies$post_start[metadata_colonies$colony=="c52"] + 9*3600 

output_name <- file.path(paste0(WORKDIR,"/scripts/EXP1_base_analysis/metadata_colonies.csv"))
write.table(metadata_colonies,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
metadata_colonies_file <- paste(WORKDIR,"scripts/EXP1_base_analysis", "metadata_colonies.csv", sep = "/")

#did it work?
metadata_colonies <- read.table(metadata_colonies_file, header = T, stringsAsFactors = F, sep = ",")
str(metadata_colonies)

metadata_colonies$pre_start[1]
