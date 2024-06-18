'''
This script creates df for random hg38 fragments
Randomly select 1000 sites
Find percentage in Exon and
percentage in Transcription Unit
'''


library(dplyr)
library(GenomicRanges)
library(stringr)
source("utilities.R")

# Read in random fragments
load("data/hg38.randomFragments.1000000.100.RData")
randomSite <- names(hg38.randomFragments.1000000.100)

# Create grange object
chromosome <- str_extract(randomSite, "chr\\w+")
position <- str_extract(randomSite, "\\d+-\\d+")

position_split <- str_split(position, "-", simplify = TRUE)

df <- data.frame(
  chromosome = chromosome,
  start = as.numeric(position_split[, 1]),
  end = as.numeric(position_split[, 2])
)

saveRDS(df, file = "data/hg38RandomGR.rds")

inExon <- readRDS("data/hg38Ref/genomeAnnotations/hg38.exons.rds")
inTU <- readRDS("data/hg38Ref/genomeAnnotations/hg38.TUs.rds")

hg38.dash <- findRandomVal(df, inExon, inTU)

saveRDS(hg38.dash, file = "reference/hg38.dash.rds")

