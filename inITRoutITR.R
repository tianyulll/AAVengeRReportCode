library(dplyr)
library(data.table)
library(ggplot2)
library(rcompanion)
library(ggmosaic)

# ITR width: 96
# Position of vector in vector plasmid.
vectorStart <- 1346
vectorEnd <- 3404
ITRwidth <- 96

# Read in result from anchorReadRearrangments module.
result <- readRDS('anchorReadRearrangements/readRearrangements.rds')

# Select reads with one or more rearrangements.
r <- data.table(result[grepl(';', result$rearrangement),] %>% select(qname, rearrangement, subject, sample))
rm(result)

# Parse each rearrangement line-by-line.
# Takes ~ 5 minutes w/ Berkely data set.
n <- 1 # Counter
w <- 5 # Wiggle

r <-  rbindlist(lapply(split(r, 1:nrow(r)), function(x){
  if(n %% 1000 == 0) message(n, '/', nrow(r))
  n <<- n+1
  a <- unlist(strsplit(x$rearrangement, ';'))[2]
  b <- stringr::str_extract(a, '\\d+[\\+\\-]\\d+')
  x$strand <- stringr::str_extract(b, '[\\+\\-]')
  d <- unlist(strsplit(b, x$strand))
  p <- as.integer(ifelse(x$strand == '+', d[1], d[2]))
  x$pos <- ifelse(p >= (vectorStart - w) & p <= (vectorEnd + w), p - vectorStart, NA)
  x
}))

r <- r[! is.na(r$pos),]
saveRDS(r, 'AnchorReadMaps.rds')


