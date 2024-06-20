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

## get available trials
getAvailTrial <- function(){
  dbConn <- dbConnect(RMariaDB::MariaDB(), group = 'AAVengeR')
  o <- dbGetQuery(dbConn, paste0("SELECT DISTINCT trial FROM sites"))
  dbDisconnect(dbConn)
  print(o)
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

# calculate length from rearrangement (repLeaderSeqMap column)
getRearrangeLength <- function(x) {
  if (!grepl(";", x)) {
    return(0)
  } else {
    matches <- str_match_all(x, "\\.\\.(\\d+)\\[")
    values <- unlist(matches)
    last_value <- tail(values, n=1)
    return(as.numeric(last_value))
  }
}

## version 1 
## no rearrangment sites are skipped
getRearrangeDf <- function(df) {
  
  o <- df %>% 
    select(sample, repLeaderSeqMap) %>%
    mutate(breaks = str_count(repLeaderSeqMap, ";")) %>%
    mutate(breaks = ifelse(is.na(breaks), 0, breaks)) %>% 
    mutate(boolBreak = ifelse(breaks == 0, 0, 1)) %>%
    mutate(rearrangeLength = sapply(repLeaderSeqMap, getRearrangeLength)) %>%
    group_by(sample) %>%
    summarise(totalBreaks = sum(breaks),
              totalLength = sum(rearrangeLength), totalBreakBool = sum(boolBreak)) %>% 
    mutate("break%" = round(totalBreaks / totalLength*100, 2),
           "breakBool%" = round(totalBreakBool / totalLength*100, 2))
  
  return(o)
}

## run a summary on aavenger output df
# getSummary <- function(df) {
#   
# }

## run a summary of rearrangement


# Create ITR remnant plots
plotRemnant <- function(df, outDir){
  
  buildAAVremnantPlots_ITRdumbellTip1 <- 125
  buildAAVremnantPlots_ITRdumbellTip2 <- 147
  remnant_colors <- c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')
  
  x <- lapply(split(df, df$sample), function(x){
    
    range <- seq(0, buildAAVremnantPlots_ITRlength, buildAAVremnantPlots_NTbinSize)
    
    d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$repLeaderSeqMap, '\\.\\.\\d+'))) + buildAAVremnantPlots_ITRseqStart,
                bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (buildAAVremnantPlots_NTbinSize/2),
                r = stringr::str_count(x$repLeaderSeqMap, ';'),
                r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n()) %>% ungroup() %>% 
      mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', expression(">= 5"))))) %>%
      droplevels() # drop unused factor levels
    
    range2 <- (range * buildAAVremnantPlots_NTbinSize) - buildAAVremnantPlots_NTbinSize
    
    p <- ggplot(d, aes(bin, n, fill = r2)) + 
      theme_bw() +
      geom_col() +
      scale_fill_manual(name = 'Recombinations', 
                        values = rev(remnant_colors[1:length(levels(d$r2))]), 
                        drop = FALSE) +
      scale_x_continuous(breaks = range,
                         labels = range2, 
                         limits = c(buildAAVremnantPlots_NTbinSize, 
                                    cut(buildAAVremnantPlots_ITRlength, breaks = c(-Inf, range, Inf), labels = FALSE) - (buildAAVremnantPlots_NTbinSize/2))) +
      scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) + 
      #geom_vline(xintercept = cut(buildAAVremnantPlots_ITRdumbellTip1, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
      #geom_vline(xintercept = cut(buildAAVremnantPlots_ITRdumbellTip2, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
      ggtitle(paste0(x$sample[1], ' | ', x$info[1], '\n', formatC(n_distinct(x$posid), format="d", big.mark=","), ' sites')) + 
      labs(x = 'ITR position', y = 'Integrations') +
      guides(fill=guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +
      theme(text = element_text(size=16), plot.title = element_text(size = 14),
            legend.title = element_text(size=12),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="bottom", plot.margin=grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
      geom_point(data = tibble(x = cut(buildAAVremnantPlots_ITRseqStart, breaks = c(-Inf, range, Inf), labels = FALSE) - (buildAAVremnantPlots_NTbinSize/2)), aes(x, 0), 
                 size = 7, shape="\u27A1", inherit.aes = FALSE) +
      coord_cartesian(clip = "off")
    
    
    if (args$saveimg) {
      ggsave(file.path(outDir,"reportPlots/remnantPlots", paste0(x$trial[1], '-', x$subject[1], '-', x$sample[1], '.png')), 
             p, dpi = 300, width = 10, height = 7, units = 'in', create.dir = T)
    }
    
    return(p)
    
  })
  
  return(x)
}
