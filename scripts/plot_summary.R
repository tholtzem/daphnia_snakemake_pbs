library(tidyverse)

getwd()
#setwd("/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/")

args <- commandArgs(trailingOnly=T)
basedir <- "depth/HiC/minid95/" # Make sure to edit this to match your $BASEDIR
### ---------------------------------------------------------------------- ###
### II. Look at summary stats
### ---------------------------------------------------------------------- ###

stats <- read.table(file = args[1], sep = '\t', header = TRUE)
#stats <- output

# subset dataframe stats
df1 <- subset(stats, mean_depth > 1)
df1

dfu1 <- subset(stats, mean_depth < 1)
dfu1$bamfile

df_u10 <- subset(stats, mean_depth < 10)
df_u10$bamfile

df10 <- subset(df1, mean_depth > 10)
df10$bamfile

df20 <- subset(df1, mean_depth > 20)
df20$bamfile

df10_u20 <- subset(df10, mean_depth < 20)
df10_u20$bamfile

# Bar lot of mean read depth per sample
#pdf("depth/plots/references_depth_hist_dfAll.pdf")
pdf(args[2])
barplot(stats$mean_depth, names=stats$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
#title("Mean read depth per sample")
title("Mean read depth per sample (all samples)")
dev.off()


# Bar lot of mean read depth per sample
#pdf("depth/plots/references_depth_hist_df10.pdf")
pdf(args[3])
barplot(df10$mean_depth, names=df10$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
#title("Mean read depth per sample")
title("Mean read depth per sample (df10)")
dev.off()

# Check distributions of all samples using boxplots
#R's boxplot function uses the standard rule to indicate an observation as a potential outlier
#if it falls more than 1.5 times the IQR (Inter-Quartile Range, calculated as Q3-Q1) below Q1 or above Q3.
#The potential outliers are plotted with circles and the Whiskers (lines that extend from Q1 and Q3 typically to the minimum and maximum) are shortened to only go as far as observations that are within 1.5*IQR of the upper and lower quartiles. 
#The box part of the boxplot is a box that goes from Q1 to Q3 and the median is displayed as a line somewhere inside the box
#pdf("depth/plots/references_depth_df1_boxplot.pdf")
pdf(args[4])
boxplot(df1$mean_depth)
title("Distribution of mean read depth (df1)")
dev.off()

#pdf("depth/plots/references_depth_df10_boxplot.pdf")
pdf(args[5])
boxplot(df10$mean_depth)
title("Distribution of mean read depth (df10)")
dev.off()

# summary
df1_sum <- summary(df1$mean_depth)
df1_sum

meanDepth_allInd <- mean(df1$mean_depth)
medianDepth_allInd <- median(df1$mean_depth)

IQRT<-IQR(df1$mean_depth)
IQRT
#IQRT<-(26.92935-19.27842)
#IQRT

# number of individuals
N = nrow(df1)
N

#1.5 times the interquartile range from the mean (median) were excluded as they are expected to be enriched for paralogs
## mean
MaxDepth1 = (meanDepth_allInd + IQRT*1.5)*N
MinDepth1  = (meanDepth_allInd - IQRT*1.5)*N
## mediam
MaxDepth2 = (medianDepth_allInd + IQRT*1.5)*N
MinDepth2  = (medianDepth_allInd - IQRT*1.5)*N
# d + 3*sqrt(d)
HengLi1_max <- (meanDepth_allInd + 3*sqrt(meanDepth_allInd))*N
HengLi2_max <- (medianDepth_allInd + 3*sqrt(meanDepth_allInd))*N

HengLi1_min <- (meanDepth_allInd - 3*sqrt(meanDepth_allInd))*N
HengLi2_min <- (medianDepth_allInd - 3*sqrt(meanDepth_allInd))*N

# alternative approach

## sum mean_depth for all individuals
total <- mean(df1$mean_depth)
# sum the square root (SQRT) of the mean depth
SQRT <- sqrt(df1$mean_depth)
# sum IQRT for all individuals
iqr <- IQR(df1$mean_depth)

## Filter1
maxFilter1 <- (total + 3*SQRT)
minFilter1 <- (total - 3*SQRT)
## Filter2
maxFilter2 <- (total + iqr*1.5)
minFilter2 <- (total - iqr*1.5)
## Filter3 - 2 times the mean
maxFilter3 <- (2*total)


filters <- data.frame(maxFilter1, minFilter1, maxFilter2, minFilter2, maxFilter3, N)
colnames(filters) <- c("HengLi_max", "HengLi_min", "IQR_max", "IQR_min", "2x_mean", "Number_of_samples")
output <- rbind(filters)
#write.table(filters, file = paste0("depth/stats/", "references_depthFilter.list"), sep = "\t", quote = F, row.names = F)
write.table(filters, file = args[6], sep = "\t", quote = F, row.names = F)
############################################################
#non-zero

dfNonZero <- summary(df1$mean_depth_nonzero)
dfNonZero

# Bar lot of mean read depth per sample
## samples having zero depth excluded
#pdf("depth/plots/references_depth_NonZero_hist_perSample.pdf")
pdf(args[7])
barplot(df1$mean_depth_nonzero, names=df1$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth (df10)", line = 3)
title("Mean read depth (non-zero) per sample (df10)")
dev.off()

##########################################################
# write samples having a mean read depth higher than 10 to a text file, without using quotes
#realigned <- paste0(df10$bamfile,".minq20.realigned.bam")
#Extract Characters Before Pattern using sub()
## first pattern, then replacement of the pattern, then string
dedupBAM_df1 <-sub(".depth.gz", "", df1$bamfile)
write.table(dedupBAM_df1, file = args[8], sep = "\n",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

dedupBAM_df10 <-sub(".depth.gz", "", df10$bamfile)
write.table(dedupBAM_df10, file = args[9], sep = "\n",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

dedupBAM_df_u10 <-sub(".depth.gz", "", df_u10$bamfile)
write.table(dedupBAM_df_u10, file = args[10], sep = "\n",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
