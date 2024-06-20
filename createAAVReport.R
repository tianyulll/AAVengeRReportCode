suppressMessages({
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(argparse)
  library(vegan)
  library(RColorBrewer)
  library(stringr)
  source("utilities.R")
})

# parameters
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "input data file in rds/excel, or name in database")
parser$add_argument("-o" ,"--outputDir", default = "output", help = "Output directory")
parser$add_argument("-t",  "--reportTitle", default = "AAV_report", help = "file name for the report")
parser$add_argument("--piNote", help = "path to text file for summary notes")
parser$add_argument("-m", "--meta", help = "path to meta data table, must have column sample and info")
parser$add_argument("--species", default = "human", required = F, help = "choose from human or mice")
parser$add_argument("-f", "--filter", required = F,  help = "only keep sites with reads greater than this threshold")
parser$add_argument("--saveimg", action = 'store_true', help = "if provided, plots will NOT be saved separately")

# these parameters will be removed in future versions
parser$add_argument("--itrStart", type="integer", default = 57, help = "itr seq start position for remnant plot")
parser$add_argument("--itrLength", type="integer", default = 197, help = "itr length for remnant plots")
parser$add_argument("--ntBinSize", type="integer", default = 3, help = "bin size for remnant plot")

args <- parser$parse_args()



# Handle zero input
if (is.null(args$input)) {
  message("Input not given; pulling available info from db")
  getAvailTrial()
  stop("Please run pipeline with one of the trials")
}
# read input as df
df <- read_data(args$input)

# to be removed in future versions
buildAAVremnantPlots_ITRlength <-  args$itrLength
buildAAVremnantPlots_ITRseqStart <- args$itrStart
buildAAVremnantPlots_NTbinSize <- args$ntBinSize

# Parse meta data
# use subject info if not provided
if (is.null(args$meta)) {
  meta <- df %>% select(sample, subject) %>% unique() %>% dplyr::rename(info = subject)
} else {
  meta <- readRDS(args$meta)
}

# run a filter
if (!is.null(args$filter)) {
  df <- df %>%
    filter(reads > args$filter)
}

# Run a summary
summary <- df %>% 
  select(trial, subject, sample, sonicLengths, reads) %>%
  mutate(count = 1) %>%
  group_by(sample, subject) %>%
  summarize(ChaoLengths = list(sonicLengths),"total reads" = sum(reads), 
            "Unique Sites" = sum(count), "inferred cell" = sum(sonicLengths)) %>%
  mutate("Chao1" = vegan::estimateR(unlist(ChaoLengths))["S.chao1"]) %>%
  select(- ChaoLengths) %>%
  dplyr::rename(patientID = subject) %>%
  left_join(y = meta, by = "sample") # merge with meta df

# Join meta-info with aavenger
meta.summary <- summary %>% 
  select(sample, info) %>%
  mutate(description = paste(sample, info, sep = "-"))

df <- df %>%
  left_join(meta, by = "sample")

# Make rearrangment summary
rearrangement <- getRearrangeDf(df)


# Create abundance plot with top 10 most abundant site
abundance <- df %>%
  # tidy-up long geneName
  mutate(nearestGene = ifelse(nchar(nearestGene) > 10, substr(nearestGene, 1, 10), nearestGene)) %>% 
  select(sample, sonicLengths, nearestGene, posid, info) %>%
  mutate(abundantCloneName = paste0(posid, "\n", nearestGene, ":", sonicLengths)) %>% 
  select(-nearestGene, -posid) %>%
  group_by(sample) %>%
  mutate(totalClone = sum(sonicLengths)) %>%
  slice_max(order_by = sonicLengths, n = 10, with_ties = F) %>%
  mutate(sonicPercent = sonicLengths / totalClone) %>%
  mutate(fileName = sample) %>%
  mutate(sample = paste(sample, info, sep = "\n")) %>%  # combine with meta-info
  select(-info)

abundantPlot <- lapply(split(abundance, abundance$sample), function(tmp){
  "
  Iterate through the abundnace dataframe by samples.
  Calculate the rest clonotypes as low abundance.
  Save the raw plots and write on-disk.
  "
  # add low abundance
  totalClone <- tmp$totalClone[[1]]
  fileName = tmp$fileName[1]
  
  tmpRow <- tibble(sample = tmp$sample[1], sonicLengths = tmp$totalClone[1] - sum(tmp$sonicLengths),
                   abundantCloneName = "Low abund", 
                   sonicPercent = (tmp$totalClone[1] - sum(tmp$sonicLengths)) / tmp$totalClone[1])
  tmp <- bind_rows(tmpRow, tmp)
  
  # order by abundance
  tmp2 <- subset(tmp, abundantCloneName != "Low abund")
  tmp2 <- tmp2[order(tmp2$sonicPercent, decreasing = TRUE),]
  tmp$abundantCloneName <- factor(tmp$abundantCloneName, 
                                  levels = c("Low abund", unique(tmp2$abundantCloneName)))
  
  p <- ggplot(tmp, aes(x = sample, y = sonicPercent, fill = abundantCloneName), ) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c('#eeeeee', brewer.pal(10, "Paired"))) +
    scale_y_continuous(labels = scales::percent) +
    annotate('text', x=1:length(totalClone), y=1.04, label=totalClone, size=2.7, hjust=0.5) +
    #geom_text(aes(label = sonicLengths), position = position_stack(vjust = 0.5),size = 3) +
    theme(panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_line(linewidth = 0.2),
          legend.text = element_text(size = 8),
          legend.title = element_blank())
  
  if (args$saveimg) {
    ggsave(filename = file.path(args$outputDir,"reportPlots/abundancePlots", paste0(fileName, '.png')), 
           plot = p, dpi = 300, create.dir = T)
  }
  
  return(p)
})

message("Abundance plots created!")


remnant.plot <- suppressMessages(plotRemnant(df, args$outputDir))
message("remnant plots created!")

# Create Gene distribution df
gene.dist <- df %>% 
  select(sample, inGene, inExon) %>%
  group_by(sample) %>%
  mutate("In Exon%" = round(mean(inExon)*100, digits = 1), 
         "In Transcription Unit%" = round(mean(inGene)*100, digits = 1)) %>%
  select(-inGene, -inExon) %>%
  unique() %>%
  left_join(y = meta.summary, by = "sample") %>%
  arrange(sample)


# Find random Value
if (args$species == "human") {
  message("Using human as default species...")
  dash <- readRDS(file = "reference/hg38.dash.rds")
} else if (args$species == "mice") {
  dash <- readRDS(file = "reference/mm9.dash.rds")
} else (
  stop("provided species is not included in reference")
)

# Create report
message("knitting reports...")

rmarkdown::render("AAV_report.Rmd",
                  output_dir = args$outputDir, output_file = paste0(args$reportTitle, ".pdf"),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = args$reportTitle))

message("done")
