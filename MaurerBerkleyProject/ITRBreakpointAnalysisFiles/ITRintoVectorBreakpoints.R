library(dplyr)
library(ggplot2)

test_data <- sites

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
  select(sample2, rearrangement, Drug) %>%
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
  mutate(breaks = str_count(rearrangement,  ";") -1) %>%
  mutate(count = 1) %>%
  mutate(boolean_breaks = breaks > 1) %>%
  group_by(sample2, Drug) %>%
  summarise(across(c(length, breaks, count, boolean_breaks), sum))  %>%
  mutate(percent_rearrangements_per_length = (breaks/length) * 1000) %>%
  mutate(boolean_rearrangement_percent = (boolean_breaks/length)*1000) %>%
  mutate(percent_number = (breaks/count)*100) %>%
  mutate(percent_boolean_number = (boolean_breaks/count)*100) %>%
  select(sample2, Drug, length, breaks, percent_rearrangements_per_length, boolean_rearrangement_percent, percent_number, percent_boolean_number) %>%
  mutate_if(is.numeric, round, digits = 2)

boolean_percent <- r3

#boolean_percent <- boolean_percent %>%
  # select(subject, sample, length, breaks, percent_rearrangements_per_length, boolean_rearrangement_percent, percent_boolean_number) %>%
  #dplyr::rename("Weighted Rearrangements/Length (%)" = percent_rearrangements_per_length) %>%
  #dplyr::rename("Boolean Rearrangements/Length (%)" = boolean_rearrangement_percent) % > %
  

boolean_graph <- ggplot(data = boolean_percent) +
  geom_bar(mapping =
             aes(x = fct_reorder(sample2, Drug),
                 y = percent_boolean_number,
                 fill = Drug),
           stat = "identity") +
  geom_text(aes(x = sample2,
                y = percent_boolean_number,
                label = percent_boolean_number),
            vjust = -0.4) +
  xlab("Sample") + ylab("% of Sequences Containing At Least One Rearrangement") +
  theme_linedraw() +
  theme(axis.text.x = element_text(face = "bold", angle = 65, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral") +
  ylim(0,100)

ggsave(filename = 'boolean_graph.png', plot = boolean_graph, width = 8, height = 6, path = "Output", device = 'png')


###Making the New Breakpoint Graphs###
x <- lapply(split(sites, sites$sample), function(x){
  message('sample: ', x$sample[1])
  
  if(x$sample[1] == 'GTSP5773') browser()
  
  range <- seq(0, 2035, 10)
  
  d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$rearrangement, '\\.\\.\\d+'))) + 34,
              bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (3/2),
              r = stringr::str_count(x$rearrangement, ';') -1,
              r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n(), per = ((n /sum(d$n))) * 100) %>% ungroup() %>% mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', '≥ 5')))) 
  
 
   ##Formatting for Labels##
  e <- d 
  e <- e %>% group_by(bin) %>% summarise(per2 = sum(per), n2 = sum(n))
  d$per2 <- ""
  d$n2 <- " "
  d$per2 <- as.double(d$per2)
  d$n2 <- as.double(d$n2)
  d <- dplyr::rows_update(d, e, by = c('bin'))
  d$per2[duplicated(d$per2)] <- ""
  d$n2[duplicated(d$n2)] <- ""
  d$per2 <- as.numeric(d$per2)
  d$n2 <- as.numeric(d$n2)
  d <- d %>% mutate_if(is.numeric, round, digits = 2)
  range2 <- range * 10
  
  p <- ggplot(d, aes(bin, n, fill = r2)) + 
    theme_bw() +
    geom_col() +
    scale_fill_manual(name = 'Recombinations After Break Into Vector', 
                      values = rev(c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')), 
                      drop = FALSE) +
    scale_x_continuous(breaks = range,
                       labels = range2) +
    scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    geom_text(aes(x = bin , y = n2, label = per2), size = 3, color = "black", vjust = -0.5) +
    ggtitle(paste0(x$Drug[1], ' | ', x$sample[1], ' | ', formatC(nrow(x), format="d", big.mark=","), ' reads')) + 
    labs(x = 'Vector position', y = 'Reads') +
    guides(fill=guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +
    theme(text = element_text(size=12), plot.title = element_text(size = 14),
          legend.title = element_text(size=8),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="bottom", plot.margin=grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
    coord_cartesian(clip = "off")
 
  z <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
    geom_point(color = "white") +
    labs(x = "test", y = "test")
  
  p2 <- z + 
    annotate("rect", xmin = cut(0, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(130, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("red"),colour = 'black', size = 0.5) +
    
    annotate("rect", xmin = cut(130, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(252, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("white"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = cut(222, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(344, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("lightblue"),colour = 'black', size = 0.5) +
    
    annotate("rect", xmin = cut(345, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(648, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("lightgray"),colour = 'black', size = 0.5) +
    
    annotate("rect", xmin = cut(649, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(852, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("darkgray"),colour = 'black', size = 0.5) +
   
    annotate("rect", xmin = cut(853, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(898, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("white"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = cut(899, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(1594, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("pink"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = cut(1595, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(1739, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("white"), colour = 'black', size = 0.5) +
    annotate("rect", xmin = cut(1740, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(1861, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("tan"),colour = 'black', size = 0.5) +
   
    annotate("rect", xmin = cut(1862, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(1929, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("white"),colour = 'black', size = 0.5) + 
    annotate("rect", xmin = cut(1930, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
             xmax = cut(2059, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
             fill = c("red"),colour = 'black', size = 0.5)  +
    
    theme_void()
  
  
  
  p3 <- p2/p + plot_layout(heights = c(1, 10))
   
   
   

  ### ggsave(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.pdf')), p, width = opt$buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
  ggsave(file.path("VectorRearrangements", paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.png')), p3, dpi = 300, width = 10, units = 'in')
  #saveRDS(p, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.rds')))
  p    
})
#saveRDS(x, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, 'plots.rds'))

#q(save = 'no', status = 0, runLast = FALSE) 
