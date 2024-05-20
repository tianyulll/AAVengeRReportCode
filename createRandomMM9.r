'''
This script creates random mm9 fragments
Randomly select 1000 sites
Find percentage in Exon and
percentage in Transcription Unit
'''

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

get_random_sites <- function(twobit_file, num_sites, site_length) {
  set.seed(0)
  # Load the genome
  genome <- import(twobit_file)
  
  # Get chromosome lengths
  chr_lengths <- seqlengths(genome)
  
  # Remove the mitochondrial chromosome (chrM) if it exists
  chr_lengths <- chr_lengths[names(chr_lengths) != "chrM"]
  
  # Initialize a list to store selected sites
  selected_sites <- vector("list", num_sites)
  
  for (i in 1:num_sites) {
    # Randomly select a chromosome
    chr_name <- sample(names(chr_lengths), 1)
    chr_length <- chr_lengths[[chr_name]]
    
    # Ensure the start position allows for the full site length
    start_pos <- sample(0:(chr_length - site_length), 1)
    end_pos <- start_pos + site_length - 1
    
    selected_sites[[i]] <- GRanges(seqnames=chr_name, ranges=IRanges(start=start_pos, end=end_pos))
  }
  
  # Combine selected sites into a single GRanges object
  selected_sites <- do.call(c, selected_sites)
  return(selected_sites)
}


random_sites <- get_random_sites("data/mm9Ref/referenceGenomes/blat/mm9.2bit", 10000, 100)

saveRDS(random_sites, file = "reference/mouseRandom10000.rds")

df <- readRDS("reference/mouseRandom10000.rds")
inExon <- readRDS("data/mm9Ref/genomeAnnotations/mm9.exons.rds")
inTU <- readRDS("data/mm9Ref/genomeAnnotations/mm9.TUs.rds")

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

mm9.dash <- findRandomVal(as.data.frame(df), inExon, inTU)

saveRDS(mm9.dash, file = "reference/mm9.dash.rds")
