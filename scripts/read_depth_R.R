library(tidyverse)

getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/")

basedir <- "depth/" # Make sure to edit this to match your $BASEDIR
bam_list <- read_lines(paste0("list/depth_dedup.list"))

#bam_list

getwd()
for (i in 1:167){
    bamfile = bam_list[i]
    # Compute depth stats
    depth <- read_tsv(paste0(basedir, bamfile), col_names = F)$X1
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
                
write.table(output,"depth_statistics.txt", sep ="\t", quote = F)

#set genome coordinates for plots, the depth files are a list of depth values for each site. Genomic coordinates for scaffolds can be e.g. calculated from the *.fai file (first column = scavvold name, second column = scaffold length). For example: gal1= 1:2950711; gal2 = 2950712-5869635...) 
#dgal1
chr <- "dgal1"
coord_start <- 1
coord_end <- 2950711
total_depth <- total_depth[coord_start:coord_end]
total_presence <- total_presence[coord_start:coord_end]

#Plot the depth distribution

pdf(file=paste0(basedir,chr,"_depth_distribution.pdf"))
#tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
tibble(total_depth = total_depth, position = coord_start:coord_end)  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point(size = 0.1)
dev.off()

# Total depth per site across all individuals 
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)

pdf(file=paste0(basedir,chr,"_depth_summary.pdf"))
total_depth_summary %>%
  ggplot(aes(x = log(total_depth), y = n)) +
  geom_point()
dev.off()
pdf(file=paste0(basedir,chr,"_presence_summary.pdf"))
total_presence_summary %>%
  ggplot(aes(x = total_presence, y = n)) +
  geom_col()
dev.off()

