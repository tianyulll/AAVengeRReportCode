library(tidyr)
library(ggplot2)
library(dplyr)
library(argparse)
library(vegan)

#library(tinytex) for latex engine

# parameters
parser <- ArgumentParser()
#parser$add_argument("--trial", help = "Trial identifier")
#parser$add_argument("--patient", help = "Patient identifier")
parser$add_argument("-i", "--input", help = "input data file in RDS")
parser$add_argument("-o" ,"--outputDir", default = "output", help = "Output directory")
parser$add_argument("-r",  "--reportTitle", default="AAV_report", help = "file name for the report")
args <- parser$parse_args()

# Read inputs
df <- readRDS("DeJongMice.rds")

# Run a summary
# Maybe include abundance?
summary <- df %>% 
  select(trial, subject, sample, sonicLengths, reads) %>%
  mutate(count = 1) %>%
  group_by(sample, subject) %>%
  summarize(ChaoLengths = list(sonicLengths),"total reads" = sum(reads), 
            "Unique Sites" = sum(count), "inferred cell" = sum(sonicLengths)) %>%
  mutate("Chao1" = vegan::estimateR(unlist(ChaoLengths))["S.chao1"]) %>%
  select(- ChaoLengths)

# Create ITR remnant plots
plotRemnant <- function(df, outDir){
  buildAAVremnantPlots_NTbinSize <- 3
  buildAAVremnantPlots_ITRlength <-  197
  buildAAVremnantPlots_ITRseqStart <- 57
  buildAAVremnantPlots_ITRdumbellTip1 <- 125
  buildAAVremnantPlots_ITRdumbellTip2 <- 147
  buildAAVremnantPlots_plotOutputWidthInches <- 10
  
  x <- lapply(split(df, df$sample), function(x){
    message('sample: ', x$sample[1])
    
    range <- seq(0, buildAAVremnantPlots_ITRlength, buildAAVremnantPlots_NTbinSize)
    
    d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$repLeaderSeqMap, '\\.\\.\\d+'))) + buildAAVremnantPlots_ITRseqStart,
                bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (buildAAVremnantPlots_NTbinSize/2),
                r = stringr::str_count(x$repLeaderSeqMap, ';'),
                r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n()) %>% ungroup() %>% 
      mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', expression(">= 5")))))
    
    range2 <- (range * buildAAVremnantPlots_NTbinSize) - buildAAVremnantPlots_NTbinSize
    
    p <- ggplot(d, aes(bin, n, fill = r2)) + 
      theme_bw() +
      geom_col() +
      scale_fill_manual(name = 'Recombinations', 
                        values = rev(c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')), 
                        drop = FALSE) +
      scale_x_continuous(breaks = range,
                         labels = range2, 
                         limits = c(buildAAVremnantPlots_NTbinSize, 
                                    cut(buildAAVremnantPlots_ITRlength, breaks = c(-Inf, range, Inf), labels = FALSE) - (buildAAVremnantPlots_NTbinSize/2))) +
      scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) + 
      geom_vline(xintercept = cut(buildAAVremnantPlots_ITRdumbellTip1, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
      geom_vline(xintercept = cut(buildAAVremnantPlots_ITRdumbellTip2, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
      ggtitle(paste0(x$subject[1], ' | ', x$sample[1], ' | ', formatC(n_distinct(x$posid), format="d", big.mark=","), ' sites')) + 
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
    ggsave(file.path(outDir, paste0(x$trial[1], '-', x$subject[1], '-', x$sample[1], '.png')), p, dpi = 300, width = buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
    p    
  })
  return(x)
}

remnant.plot <- plotRemnant(df, paste0(args$outputDir, "/","reportPlots/remnantPlots"))

#
hm <- df %>% 
  select(sample, inGene, inExon) %>%
  group_by(sample) %>%
  mutate(ExonPer = mean(inExon)*100, GenePer = mean(inGene)*100) %>%
  select(-inGene, -inExon) %>%
  unique()

# Create report
rmarkdown::render("AAV_report.Rmd",
                  output_dir = args$outputDir,
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = args$r))
