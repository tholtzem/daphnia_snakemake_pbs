# use this script to read depth, plot summary stats and calculate depth divers filters

library(tidyverse)

getwd()
#setwd("/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs")


args <- commandArgs(trailingOnly=T)
basedir <- "depth/HiC/minid95/" # Make sure to edit this to match your $BASEDIR

### ---------------------------------------------------------------------- ###
### I. Read depth
### ---------------------------------------------------------------------- ###

bam_list <- sort(read_lines(paste0(args[1]))) # go through file line by line and sort
bam_list
df <- read_csv(args[1], col_names=F, show_col_types = F) # also import as df to get the number of rows
df

for (i in 1:nrow(df)){
    bamfile = bam_list[i]
    bamfile
    # Compute depth stats
    depth <- read_tsv(paste0(basedir, bamfile), col_names = F, show_col_types = F)$X1
    mean_depth <- mean(depth)
    sd_depth <- sd(depth)
    mean_depth_nonzero <- mean(depth[depth > 0])
    mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
    median <- median(depth)
    presence <- as.logical(depth)
    proportion_of_reference_covered <- mean(presence)
      
  # Bind stats into dataframe and store sample-specific per base depth and presence data
    if (i==1){
      output <- data.frame(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
      total_depth <- depth
      total_presence <- presence
    } else {
      output <- rbind(output, cbind(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered))
      total_depth <- total_depth + depth
      total_presence <- total_presence + presence
    }
}

output %>%
  mutate(across(where(is.numeric), round, 3))

# Write summary stats to file
write.table(output, args[2], sep ="\t", quote = F, row.names = F)

