library(tidyr)
library(ggplot2)
library(dplyr)
library(argparse)
library(vegan)
library(RColorBrewer)

# parameters
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "input data file in RDS", required = T)
parser$add_argument("-o" ,"--outputDir", default = "output", help = "Output directory")
parser$add_argument("-r",  "--reportTitle", default = "AAV_report", help = "file name for the report")
parser$add_argument("--itrStart", type="integer", default = 57, help = "itr seq start position for remnant plot")
parser$add_argument("--itrLength", type="integer", default = 197, help = "itr length for remnant plots")
parser$add_argument("--ntBinSize", type="integer", default = 3, help = "bin size for remnant plot")
parser$add_argument("--piNote", help = "path to text file for summary notes")
parser$add_argument("--meta", help = "path to meta data table, must have column sample and info")
parser$add_argument("-s",  "--species", default = "human", help = "choose from human or mice")
args <- parser$parse_args()

# Read inputs
df <- readRDS(args$input)
buildAAVremnantPlots_ITRlength <-  args$itrLength
buildAAVremnantPlots_ITRseqStart <- args$itrStart
buildAAVremnantPlots_NTbinSize <- args$ntBinSize

# Parse meta data
if (is.null(args$meta)) {
  meta <- df %>% select(sample) %>% unique()
  meta$info <- meta$sample
} else {
  meta <- readRDS(args$meta)
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
  rename(patientID = subject) %>%
  left_join(y = meta, by = "sample") # merge with meta df

# Join meta-info with aavenger
meta.summary <- summary %>% 
  select(sample, info) %>%
  mutate(description = paste(sample, info, sep = "-"))

df <- df %>%
  left_join(meta, by = "sample")
  
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
  return(p)
  ggsave(file.path(args$outputDir,"reportPlots/abundancePlots", paste0(tmp$sample[1], '.png')), 
         p, dpi = 300, create.dir = T)
})

# Create ITR remnant plots
plotRemnant <- function(df, outDir){
  
  buildAAVremnantPlots_ITRdumbellTip1 <- 125
  buildAAVremnantPlots_ITRdumbellTip2 <- 147
  remnant_colors <- c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')
  
  x <- lapply(split(df, df$sample), function(x){
    message('sample: ', x$sample[1])
    
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
    
    ggsave(file.path(outDir,"reportPlots/abundancePlots", paste0(x$trial[1], '-', x$subject[1], '-', x$sample[1], '.png')), 
           p, dpi = 300, width = 10, height = 7, units = 'in', create.dir = T)
    p    
  })
  return(x)
}

remnant.plot <- plotRemnant(df, args$outputDir)

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
if (args$s == "human") {
  dash <- readRDS(file = "reference/hg38.dash.rds")
} else if (args$s == "mice") {
  dash <- readRDS(file = "reference/mm9.dash.rds")
} else (
  stop("this species is not included in reference")
)

# Create report
rmarkdown::render("AAV_report.Rmd",
                  output_dir = args$outputDir, output_file = paste0(args$r, ".pdf"),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = args$r))
