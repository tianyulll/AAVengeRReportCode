'''
This script creates df for random hg38 fragments
Randomly select 1000 sites
Find percentage in Exon and
percentage in Transcription Unit
'''


library(dplyr)
library(GenomicRanges)
library(stringr)

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

findRandomVal <- function(df, inExon, inTU){
  # Create random sites gr object
  set.seed(1)
  df_subset <- df[sample(nrow(df), 1000), ]
  randomSite.gr <- makeGRangesFromDataFrame(df_subset)
  # is the integration in exon
  isExon <- findOverlaps(randomSite.gr, inExon, type="within", select = "arbitrary") %>%
    sapply(function(x){ifelse(is.na(x), 0, 1)})
  # is the integration in TU
  isTU <- findOverlaps(randomSite.gr, inTU, type="within", select = "arbitrary") %>%
    sapply(function(x){ifelse(is.na(x), 0, 1)})
  return(c(sum(isExon)/length(isExon)*100, sum(isTU)/length(isTU)*100))
}

hg38.dash <- findRandomVal(df, inExon, inTU)

saveRDS(hg38.dash, file = "reference/hg38.dash.rds")

