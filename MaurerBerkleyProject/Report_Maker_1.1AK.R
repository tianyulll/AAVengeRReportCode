options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(RMySQL)  # loads DBI
library(lubridate)
library(gtools) 
library(wordcloud) 
library(ggplot2) 
library(reldist) 
library(vegan) 
library(argparse)
library(dplyr)
library(stringr)
library(readr)
library(tidyverse)
library(tidyr)
library(IntegrationFeatureHeatmap)
library(gtools)
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation
library("readxl") # for reading in excel
library(forcats)
library(rcompanion)

# Set up input arguments
# 
#-----------------------------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--Experiment", type = "character", default = "AAV", help = "Experimental identifier")
parser$add_argument("--WorkingDirectory", type = "character", default = "/Users/bushmangradstudent/Desktop/BushmanLab/Report_Maker", 
                    help = "Directory where your reportmaker.R file is located")
parser$add_argument("--outputDir", type = "character", default = "/Users/bushmangradstudent/Desktop/BushmanLab/Report_Maker/Output", help = "Output directory")
parser$add_argument("--dataDir", type = "character", default = "/Users/bushmangradstudent/Desktop/BushmanLab/Report_Maker/Data", help = "Direction where data is")
parser$add_argument("--rmdFile", type="character", default="/Users/bushmangradstudent/Desktop/BushmanLab/Report_Maker/AAV_Baseline.rmd", help="Rmd file used to build reports")
parser$add_argument("--SubjectNameConversions", type = "character", default = "Subject_Key.csv", help = "File containing subject name conversions")

args <- parser$parse_args()

#what is the Bogus Control Called
Bogus <- "AAVBogusControl_230119"
nBogus <- "Artifacts"

#Delineate which GTSPs are in which experimental group
G1 <- c(5786,5791)
G2 <- c(5787,5792)
G3 <- c(5788,5793)
G4 <- c(5789,5794)
G5 <- c(5790,5795)


# name the groups
nG1 <- "HDRiCrude"
nG2 <- "DMSOCrude"
nG3 <- "Negative"
nG4 <- "HDRiPure"
nG5 <- "DMSOPure"


# Read in files
#-----------------------------------------------------------------------------------------------
#subject_Key <- read.csv("Subject_Key.csv")
sites_data <- read.csv("Data/ Berkley_AAV_20231003.csv")
anchor_read <- read.csv("Data/AAV_anchor_read_rearrangements.csv")

# Filter data
# Filtering based on sonic lengths and reads (both should be >2)
 #filtered_data <- sites_data %>% dplyr::filter(sonicLengths > 2 & reads > 2) #removing reads with fewer than 2 soniclengths or reads
#-----------------------------------------------------------------------------------------------

#Messed-upness from repLeaderSeqMap
#-----------------------------------------------------------------------------------------------
test_data <- sites_data

test_data$sample<-gsub("GTSP","",as.character(test_data$sample))
test_data$sample<-gsub(Bogus,"1111", as.character(test_data$sample))
test_data$sample<-gsub(Bogus,"1111", as.character(test_data$sample))

#Delineate which GTSPs are in which experimental group
test_data$Drug <- ""
test_data$Drug[test_data$sample %in% G1] <- "HDRiCrude"
test_data$Drug[test_data$sample %in%  G2] <- "DMSOCrude"
test_data$Drug[test_data$sample %in%  G3] <- "Negative"
test_data$Drug[test_data$sample %in%  G4] <- "HDRiPure"
test_data$Drug[test_data$sample %in%  G5] <- "DMSOPure"

for (i in 1:nrow(test_data)) {
  if(test_data[i, "Drug"] == ""){
    test_data[i, "Drug"] = "Positive Control"
  }
}

Test <- select(test_data, Drug, subject, sample)

if(!"rearrangement" %in% colnames(test_data)){
  test_data <- test_data %>%
    mutate(rearrangement = repLeaderSeqMap)
}

r1 <- test_data %>%
  select(trial, subject, sample, rearrangement, Drug) %>%
  mutate(last = sapply(strsplit(rearrangement, ";"), tail, 1))

store_length <- c() #storing our values
for (i in 1:nrow(r1)){ #initialize for loop
  row <- r1[i, ] #grabbing the row we're on
  column <- row$last #grabbing the actual value we want
  if (is.na(column) == TRUE){ #if we have an NA
    store_length <- c(store_length, 0) #we are going to set the value to 0
  } else if (column == "") {
    store_length <- c(store_length, 0)
  } else if (lengths(column) == 0) {
    store_length <- c(store_length, 0)
  } else { #otherwise
    length <- regmatches(column, regexpr("\\.\\.([0-9]+)\\[", row$last)) #grab the number between the '..' and '['
    length <- stringr::str_replace(length, '..', '') #replace ..
    length <- stringr::str_replace(length, '\\[', '') #replace [
    length <- as.numeric(length) #change the value to a number rather than a string
    store_length <- c(store_length, length) #store the value
  }
}

r1$length <- store_length #add a new column with our length

r2 <- r1 %>% #we're going to replace all the NA values in our database in the last and repLeaderSeqMap
  mutate(rearrangement = as.character(rearrangement)) %>%
  mutate_at(c("rearrangement", "last"), ~replace_na(.,"0")) 

r3 <- r2 %>%
  mutate(breaks = str_count(rearrangement,  ";")) %>%
  mutate(count = 1) %>%
  mutate(boolean_breaks = grepl(";", rearrangement)) %>%
  group_by(subject, Drug) %>%
  summarise(across(c(length, breaks, count, boolean_breaks), sum)) %>%
  mutate(percent_rearrangements_per_length = (breaks/length) * 1000) %>%
  mutate(boolean_rearrangement_percent = (boolean_breaks/length)*1000) %>%
  mutate(percent_number = (breaks/count)*100) %>%
  mutate(percent_boolean_number = (boolean_breaks/count)*100) %>%
  select(subject, Drug, length, breaks, percent_rearrangements_per_length, boolean_rearrangement_percent, percent_number, percent_boolean_number) %>%
  mutate_if(is.numeric, round, digits = 2)

boolean_percent <- r3

boolean_percent <- boolean_percent %>%
  # select(subject, sample, length, breaks, percent_rearrangements_per_length, boolean_rearrangement_percent, percent_boolean_number) %>%
  dplyr::rename("Weighted Rearrangements/Length (%)" = percent_rearrangements_per_length) %>%
  dplyr::rename("Boolean Rearrangements/Length (%)" = boolean_rearrangement_percent)

boolean_graph <- ggplot(data = boolean_percent) +
  geom_bar(mapping =
             aes(x = fct_reorder(subject, Drug),
                 y = percent_boolean_number,
                 fill = Drug),
           stat = "identity") +
  geom_text(aes(x = subject,
                y = percent_boolean_number,
                label = percent_boolean_number),
            vjust = -0.4) +
  xlab("Sample") + ylab("% of Sequences Containing At Least One Rearrangement") +
  theme_linedraw() +
  theme(axis.text.x = element_text(face = "bold", angle = 65, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral") +
  ylim(0,100)

ggsave(filename = 'boolean_graph.png', plot = boolean_graph, width = 8, height = 6, path = "Output", device = 'png')
#-----------------------------------------------------------------------------------------------


# Messed-upness from anchorReadRearrangements -- using the UMI rearranged % given in the data -- not making any alterations
#-----------------------------------------------------------------------------------------------
ARR_Data <- anchor_read

ARR_Data <- ARR_Data %>% filter(trial == "pRep")

ARR_Data$sample<-gsub("GTSP","",as.character(ARR_Data$sample))
ARR_Data$sample<-gsub(Bogus,"1111", as.character(ARR_Data$sample))
ARR_Data$sample<-gsub(Bogus,"1111", as.character(ARR_Data$sample))

#Delineate which GTSPs are in which experimental group
ARR_Data$Drug <- ""
ARR_Data$Drug[ARR_Data$sample %in% G1] <- "HDRiCrude"
ARR_Data$Drug[ARR_Data$sample %in%  G2] <- "DMSOCrude"
ARR_Data$Drug[ARR_Data$sample %in%  G3] <- "Negative"
ARR_Data$Drug[ARR_Data$sample %in%  G4] <- "HDRiPure"
ARR_Data$Drug[ARR_Data$sample %in%  G5] <- "DMSOPure"

for (i in 1:nrow(ARR_Data)) {
  if(ARR_Data[i, "Drug"] == ""){
    ARR_Data[i, "Drug"] = "Positive Control"
  }
}

ARR_Data$rearrangement <- ""
for (i in 1:nrow(ARR_Data)) {
  ARR_Data[i, "rearrangement"] <- as.numeric(sub("%", "", ARR_Data[i, "percentUMIsRearranged"]))
}

ARR_Data$rearrangement <- as.numeric(ARR_Data$rearrangement)
anchor_graph <- ggplot(data = ARR_Data) +
  geom_bar(mapping =
             aes(x = fct_reorder(subject, Drug),
                 y = rearrangement,
                 fill = Drug),
           stat = "identity") +
  geom_text(aes(x = subject,
                y = rearrangement,
                label = percentUMIsRearranged),
            vjust = -0.4) +
  xlab("GTSP Number") + ylab("UMIs Rearranged (%)") +
  theme_linedraw() +
  theme(axis.text.x = element_text(face = "bold", angle = 65, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral") +
  ylim (0,100)
  

ggsave(filename = 'anchor_read_graph.png', plot = anchor_graph,width = 8, height = 6, path = "Output", device = 'png')

#-----------------------------------------------------------------------------------------------


# Combined Graphs with Error Bars
#-----------------------------------------------------------------------------------------------
aData <- ARR_Data %>% select(Drug, rearrangement)
bData <- boolean_percent %>% 
  ungroup() %>%
  select(Drug, percent_boolean_number)

for (i in 1:nrow(aData)) {
  if (grepl("DMSO", aData[i, "Drug"])) {
    aData[i, "Drug"] <- "DMSO"
  } else if (grepl("HDRi", aData[i, "Drug"])) {
    aData[i, "Drug"] <- "HDRi"
  }
}

for (i in 1:nrow(bData)) {
  if (grepl("DMSO", bData[i, "Drug"])) {
    bData[i, "Drug"] <- "DMSO"
  } else if (grepl("HDRi", bData[i, "Drug"])) {
    bData[i, "Drug"] <- "HDRi"
  }
}

aData <- aData %>%
  group_by(Drug) %>%
  summarize(rearrangement = list(rearrangement))

bData <- bData %>%
  group_by(Drug) %>%
  summarize(rearrangement = list(percent_boolean_number))

aData$sd <- 0
for (i in 1:nrow(aData)) {
  aData[i, "sd"] <- sd(unlist(aData[i, "rearrangement"]))
}

bData$sd <- 0
for (i in 1:nrow(bData)) {
  bData[i, "sd"] <- sd(unlist(bData[i, "rearrangement"]))
}

#-----------------------------------------------------------------------------------------------

# Table with Chao1 estimates
#-----------------------------------------------------------------------------------------------

chao_data <- sites_data %>% select(trial, subject, sample, UMIs, sonicLengths, reads)

chao_data$sample<-gsub("GTSP","",as.character(chao_data$sample))
chao_data$sample<-gsub(Bogus,"1111", as.character(chao_data$sample))
chao_data$sample<-gsub(Bogus,"1111", as.character(chao_data$sample))

#Delineate which GTSPs are in which experimental group
chao_data$Drug <- ""
chao_data$Drug[chao_data$sample %in% G1] <- "HDRiCrude"
chao_data$Drug[chao_data$sample %in%  G2] <- "DMSOCrude"
chao_data$Drug[chao_data$sample %in%  G3] <- "Negative"
chao_data$Drug[chao_data$sample %in%  G4] <- "HDRiPure"
chao_data$Drug[chao_data$sample %in%  G5] <- "DMSOPure"

for (i in 1:nrow(chao_data)) {
  if(chao_data[i, "Drug"] == ""){
    chao_data[i, "Drug"] = "Positive Control"
  }
}

group_test <- chao_data %>%
  mutate(count = 1) %>%
  group_by(subject, Drug, sample) %>%
  summarize(ChaoLengths = list(sonicLengths), UMIs = sum(UMIs), "Total Reads" = sum(reads), "Unique Sites" = sum(count), "Inferred Cells" = sum(sonicLengths))

stored_data <- c()
for (i in 1:nrow(group_test)) {
  alpha <- pull(group_test[i, "ChaoLengths"])
  alpha <- unlist(alpha)
  chao <- vegan::estimateR(alpha)
  stored_data <- c(stored_data, chao["S.chao1"])
}

group_test$Chao1 <- stored_data

chao_table <- select(group_test, subject, Drug, UMIs, "Total Reads", "Unique Sites", Chao1, "Inferred Cells")

#Chromosomal Site Locations
#-----------------------------------------------------------------------------------------------

# Set up the integration sites. Make sure all unwanted sites and yeast positive control sites have been removed from input file
rawtable <- read_excel("Data/Berkley_AAV_20231003_2.xlsx")
x_table <- list(rawtable$sample, rawtable$posid)
table <- as.data.frame(x_table, col.names = c('sample', 'posid'))
table <- table %>%
  separate(posid, c("chromosome", "start"))
table$sample<-gsub("GTSP","",as.character(table$sample))
table$sample<-gsub(Bogus,"1111", as.character(table$sample))
table$sample<-gsub(Bogus,"1111", as.character(table$sample))
table <- table %>% filter(!sample == "AAVPosControl_CL_230815")
len <- c(nrow(table))
sites <- structure(list(chromosome = c(table$chromosome),
                        start = as.integer(c(table$start)), 
                        gtsp = as.integer(c(table$sample))), .Names = c("chromosome", "start", "gtsp"), row.names = c(NA, len), class = "data.frame")

#Delineate which GTSPs are in which experimental group
sites$Drug <- ""
sites$Drug[sites$gtsp %in% G1] <- "HDRiCrude"
sites$Drug[sites$gtsp %in%  G2] <- "DMSOCrude"
sites$Drug[sites$gtsp %in%  G3] <- "Negative"
sites$Drug[sites$gtsp %in%  G4] <- "HDRiPure"
sites$Drug[sites$gtsp %in%  G5] <- "DMSOPure"

# create a color key for the plot
group.colors <- c(HDRiCrude = "yellow", DMSOCrude = "red", Negative = "blue", HDRiPure = "green", "DMSOPure" = "orange")
# hg38 chromosome sizes
chrom_sizes <- structure(list(V1 = c( "chr1", "chr2", "chr3", "chr4", 
                                      "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                                      "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                      "chr20", "chr21", "chr22", "chrX", "chrY"), V2 = c( 248956422L, 
                                                                                          242193529L, 198295559L, 190214555L, 181538259L, 170805979L, 159345973L, 
                                                                                          145138636L, 138394717L, 133797422L, 135086622L, 133275309L, 114364328L, 
                                                                                          107043718L, 101991189L, 90338345L, 83257441L, 80373285L, 58617616L, 
                                                                                          64444167L, 46709983L, 50818468L, 156040895L, 57227415L)), .Names = c("V1", 
                                                                                                                                                               "V2"), class = "data.frame", row.names = c(NA, -24L))



# hg38 centromere locations
centromeres <- structure(list(X.bin = c(23L, 20L, 2L, 1L, 14L, 16L, 1L, 14L, 
                                        1L, 1L, 10L, 1L, 15L, 13L, 1L, 1L, 11L, 13L, 1L, 1L, 1L, 12L, 
                                        10L, 10L), chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", 
                                                             "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
                                                             "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                                                             "chr19", "chr20", "chr21", "chr22","chrX", "chrY"), chromStart = c(121700000L, 91800000L, 87800000L, 48200000L, 46100000L, 58500000L, 
                                                                                                                                58100000L, 43200000L, 42200000L, 38000000L, 51000000L, 33200000L, 
                                                                                                                                16500000L, 16100000L, 17500000L, 35300000L, 22700000L, 15400000L, 
                                                                                                                                24200000L, 25700000L, 10900000L, 13700000L, 58100000L, 10300000L), chromEnd = c(125100000L, 96000000L, 94000000L, 51800000L, 51400000L, 62600000L, 
                                                                                                                                                                                                                62100000L, 47200000L, 45500000L, 41600000L, 55800000L, 37800000L, 
                                                                                                                                                                                                                18900000L, 18200000L, 20500000L, 38400000L, 27400000L, 21500000L, 
                                                                                                                                                                                                                28100000L, 30400000L, 13000000L, 17400000L, 63800000L, 10600000L), ix = c(1270L, 
                                                                                                                                                                                                                                                                                          770L, 784L, 447L, 452L, 628L, 564L, 376L, 411L, 583L, 105L, 341L, 
                                                                                                                                                                                                                                                                                          447L, 304L, 3L, 3L, 3L, 354L, 192L, 125L, 410L, 275L, 22L, 3L
                                                                                                                                                                                                                ), n = c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                                                                                                                                                                                         "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N"
                                                                                                                                                                                                                ), size = c(3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
                                                                                                                                                                                                                            3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
                                                                                                                                                                                                                            3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
                                                                                                                                                                                                                            3000000L, 3000000L, 3000000L, 3000000L, 3000000L), type = c("centromere", 
                                                                                                                                                                                                                                                                                        "centromere", "centromere", "centromere", "centromere", "centromere", 
                                                                                                                                                                                                                                                                                        "centromere", "centromere", "centromere", "centromere", "centromere", 
                                                                                                                                                                                                                                                                                        "centromere", "centromere", "centromere", "centromere", "centromere", 
                                                                                                                                                                                                                                                                                        "centromere", "centromere", "centromere", "centromere", "centromere", 
                                                                                                                                                                                                                                                                                        "centromere", "centromere", "centromere"), bridge = c("no", "no", 
                                                                                                                                                                                                                                                                                                                                              "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", 
                                                                                                                                                                                                                                                                                                                                              "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no"
                                                                                                                                                                                                                                                                                        )), .Names = c("X.bin", "chrom", "chromStart", "chromEnd", "ix", 
                                                                                                                                                                                                                                                                                                       "n", "size", "type", "bridge"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                     -24L))


# set the column names for the datasets
# IMPORTANT: fields common across datasets should have the same name in each
colnames(chrom_sizes) <- c("chromosome", "size")
colnames(centromeres) <- c('bin', "chromosome", 'start', 'end',
                           'ix', 'n', 'size', 'type', 'bridge')

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX", "chrY")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23, 24)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)
sites[["chromosome"]] <- factor(x = sites[["chromosome"]], 
                                levels = chrom_order)
centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)



#-----------------------------------------------------------------------------------------------
#ITRinITR and ITRinVector

# ITR width: 96
# Position of vector in vector plasmid.
vectorStart <- 1346
vectorEnd <- 3404
ITRwidth <- 96
r <- readRDS("Data/AnchorReadMaps.rds")
r$sample2<-gsub("GTSP","",as.character(r$sample))

#Delineate which GTSPs are in which experimental group
r$Drug <- ""
r$Drug[r$sample2 %in% G1] <- nG1
r$Drug[r$sample2 %in%  G2] <- nG2
r$Drug[r$sample2 %in%  G3] <- nG3
r$Drug[r$sample2 %in%  G4] <- nG4
r$Drug[r$sample2 %in%  G5] <- nG5


distribution <- ggplot(r, aes(pos)) +
  theme_bw() +
  geom_density() +
  geom_vline(xintercept = c(0, ITRwidth, vectorEnd-vectorStart-ITRwidth, vectorEnd-vectorStart), color = 'red') +
  facet_wrap(~Drug, ncol = 2, scales = 'free_y') +
  ylab("Density of Reads") +
  xlab("Vector Position (bp)") 

ggsave(filename = 'distribution.png', plot = distribution, width = 8, height = 6, path = "Output", device = 'png')


x <- c(nG1,nG2,nG4,nG5)
inITR <- c(1:ITRwidth,(vectorEnd-vectorStart-ITRwidth):(vectorEnd-vectorStart))


ratios <- c()
for (i in x) {
  r1 <- r %>% 
    filter(r$Drug == i)
  iITR = r1 [r1$pos %in% inITR]
  oITR = r1 [!r1$pos %in% inITR]
  withinITR <- nrow(iITR)
  outsideITR <- nrow(oITR)
  values <- c(withinITR,outsideITR)
  ratios = rbind(ratios,values)
  r1 <- r
}

ratios <- cbind(x,ratios)
rownames(ratios) <- x

colnames(ratios) <- c("Drug","ITRtoITR", "ITRtoVector")

ratios <- as.data.frame(ratios)
ratios$ITRtoITR <- as.numeric(ratios$ITRtoITR)
ratios$ITRtoVector <- as.numeric(ratios$ITRtoVector)

Crude <-dplyr::filter(ratios,grepl('Crude',Drug))
Pure <-dplyr::filter(ratios,grepl('Pure',Drug))

Crude <- Crude[-1]
Pure <- Pure [-1]

Crude <- as.matrix(Crude)
Pure <- as.matrix(Pure)
Total <- Crude + Pure
rownames(Total) <- c("HDRi","DMSO")

png("Output/MosaicCrude.png",width=3.25,height=3.25,units="in",res=1200)
par(cex.main= 0.8)
mosaicplot(Crude,
           main = "Mosaic Plot of Crude Samples",
           color = TRUE)
dev.off()

png("Output/MosaicPure.png",width=3.25,height=3.25,units="in",res=1200)
par(cex.main= 0.8)
mosaicplot(Pure,
           main = "Mosaic Plot of Pure Samples",
           color = TRUE)
dev.off()


png("Output/MosaicTotal.png",width=3.25,height=3.25,units="in",res=1200)
par(cex.main= 0.8)
mosaicplot(Total,
           main = "Mosaic Plot by Drug Condition",
           color = TRUE)
dev.off()

p.values <- c(chisq.test(Crude)$p.value,chisq.test(Pure)$p.value, chisq.test(Total)$p.value)
cramers <- c(cramerV(Crude),cramerV(Pure),cramerV(Total))
statstable <- cbind(p.values,cramers)
rownames(statstable) <- c("Crude Prep Comparison", "Pure Prep Comparison", "Total Drug Comaprison")
colnames(statstable) <- c("P-value by chisquared","CramersV")




#Heatmaps
#-----------------------------------------------------------------------------------------------


# Make PDF
#-----------------------------------------------------------------------------------------------
chaoTable <- chao_table
rmarkdown::render("AAV_Baseline.Rmd",
                  output_dir = "Output",
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = paste0("Anna Maurer HDR inhibition study"),
                                'chromosome' = centromeres,
                                'chrom_sizes' = chrom_sizes,
                                'sites' = sites,
                                'chaoTable' = chao_table))

