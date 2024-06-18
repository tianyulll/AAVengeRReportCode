'''
This script creates random mm9 fragments
Randomly select 1000 sites
Find percentage in Exon and
percentage in Transcription Unit
'''

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
source("utilities.R")

random_sites <- get_random_sites("data/mm9Ref/referenceGenomes/blat/mm9.2bit", 10000, 100)

saveRDS(random_sites, file = "reference/mouseRandom10000.rds")

df <- readRDS("reference/mouseRandom10000.rds")
inExon <- readRDS("data/mm9Ref/genomeAnnotations/mm9.exons.rds")
inTU <- readRDS("data/mm9Ref/genomeAnnotations/mm9.TUs.rds")

mm9.dash <- findRandomVal(as.data.frame(df), inExon, inTU)

saveRDS(mm9.dash, file = "reference/mm9.dash.rds")
