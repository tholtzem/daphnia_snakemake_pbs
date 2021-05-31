rm(list=c(ls()))
setwd("/home/tania/github/daphnia_snakemake_pbs/")
getwd()

library(Sushi)
library(devtools)

ALG1= read.table( "bedtools/ALG1.test.bed", sep='\t')
genome = read.table( "ref/dgal_ra.genome")
head(ALG1)
head(genome)

plotManhattan(bedfile = ALG1, pvalues = ALG1[,5], col = SushiColors(6), genome = genome, cex=0.75 )

ALG1[,5]
