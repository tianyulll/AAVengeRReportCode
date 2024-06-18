library(tidyr)
library(ggplot2)
library(dplyr)
library(argparse)
library(vegan)
library(RColorBrewer)
library(RMariaDB)
library(readxl)
library(rtracklayer)

# find the random in exon
# and in gene percentange 
# from the random integration sites
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

# creates random integration sites from a species
# whole genome reference 
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

reconstructDBtable <- function(o, tmpDirPath){
  r <- tibble()
  
  if(nrow(o) > 0){
    r <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
      f <- tmpFile()
      writeBin(unserialize(y$data[[1]]), file.path(tmpDirPath, paste0(f, '.xz')))
      system(paste0('unxz ', file.path(tmpDirPath, paste0(f, '.xz'))))
      d <- readr::read_tsv(file.path(tmpDirPath, f))
      invisible(file.remove(file.path(tmpDirPath, f)))
      d
    }))
  }
  
  r
}

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

pullDBsubjectSites <- function(dbConn, trial, subject, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from sites where trial = '", trial, "' and subject = '", subject, "'"))
  reconstructDBtable(o, tmpDirPath)
}

pullDBTrialSites <- function(dbConn, trial, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from sites where trial = '", trial, "'"))
  reconstructDBtable(o, tmpDirPath)
}

## get AAVengeR results
read_data <- function(file) {
  file_extension <- tools::file_ext(file)
  
  if (file_extension == "rds") {
    message("reading rds as df")
    data <- readRDS(file)
  } else if (file_extension %in% c("xls", "xlsx")) {
    message("reading excel as df")
    data <- read_excel(file)
  } else {
    message("name given - trying pulling from database")
    dbConn <- dbConnect(RMariaDB::MariaDB(), group = 'AAVengeR')
    data <- suppressMessages(pullDBTrialSites(dbConn, file, './'))
    dbDisconnect(dbConn) 
  }
  
  if ( nrow(data) == 0 ) {
    stop("reading NOT successful")
  } else(
    message("Success! Number of rows read in:", nrow(data))
  )
  
  return(data)
}

## run a summary on aavenger output df
# getSummary <- function(df) {
#   
# }

## run a summary of rearrangement
