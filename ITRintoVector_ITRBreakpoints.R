library(dplyr)
library(ggplot2)
library(reshape2)

sites <- read_rds("AnchorReadMaps.rds")
test_data <- sites

test_data$sample<-gsub("GTSP","",as.character(test_data$sample))
#Delineate which GTSPs are in which experimental group
G1 <- c(5786,5791)
G2 <- c(5787,5792)
G3 <- c(5788,5793)
G4 <- c(5789,5794)
G5 <- c(5790,5795)

#Delineate which GTSPs are in which experimental group
test_data$Drug <- ""
test_data$Drug[test_data$sample %in% G1] <- "HDRiCrude"
test_data$Drug[test_data$sample %in%  G2] <- "DMSOCrude"
test_data$Drug[test_data$sample %in%  G3] <- "Negative"
test_data$Drug[test_data$sample %in%  G4] <- "HDRiPure"
test_data$Drug[test_data$sample %in%  G5] <- "DMSOPure"


###Making the New Breakpoint Graphs###
x <- lapply(split(sites, sites$sample), function(x){
  message('sample: ', x$sample[1])
  
  if(x$sample[1] == 'GTSP5773') browser()
  
  range <- seq(0, 130, 3)
  
  b <- nrow(x)
  
  d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$rearrangement, '\\.\\.\\d+'))) + 34,
              bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (3/2),
              r = stringr::str_count(x$rearrangement, ';') -1,
              r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n(),per = ((n/b)) * 100) %>% ungroup() %>% mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', '≥ 5')))) 
  
  
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
  d$per2 <- lapply(d$per2, paste0, '%')
  
  
  range2 <- (range * 3)-3 
  
  
  p <- ggplot(d, aes(bin, n, fill = r2)) + 
    theme_bw() +
    geom_col() +
    scale_fill_manual(name = 'Recombinations after Vector Transition', 
                      values = rev(c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')), 
                      drop = FALSE) +
    scale_x_continuous(breaks = range,
                       labels = range2, 
                       limits = c(3, 
                                  cut(130, breaks = c(-Inf, range, Inf), labels = FALSE) - (3/2))) +
    scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) + 
    geom_text(aes(x = bin , y = n2+b*0.002, label = per2), size = 3, color = "black",hjust = 0, angle = 90) +
    geom_vline(xintercept = cut(71, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
    geom_vline(xintercept = cut(93, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
    ggtitle(paste0(x$Drug[1], ' | ', x$sample[1], ' | ', formatC(nrow(x), format="d", big.mark=","), ' reads')) + 
    labs(x = 'ITR position', y = '# of Vector Positions') +
    guides(fill=guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +
    theme(text = element_text(size=16), plot.title = element_text(size = 14),
          legend.title = element_text(size=12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="bottom", plot.margin=grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
    geom_point(data = tibble(x = cut(34, breaks = c(-Inf, range, Inf), labels = FALSE) - (3/2)), aes(x, 0), 
               size = 7, shape="\u27A1", inherit.aes = FALSE) +
    coord_cartesian(clip = "off") 
  
  z <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
    geom_point(color = "white") +
    labs(x = "test", y = "test") 
  
  
  p2 <- z +
    annotate("rect", xmin = 0, ymin = -Inf,
             xmax = (52/15), ymax = Inf,
             fill = c("red"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = (52/15), ymin = -Inf,
             xmax = (64/15), ymax = Inf,
             fill = c("lightblue"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = (64/15), ymin = -Inf,
             xmax = (75/15), ymax = Inf,
             fill = c("blue"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = (75/15), ymin = -Inf,
             xmax = (85/15), ymax = Inf,
             fill = c("lightgreen"),colour = 'black', size = 0.5) +
    annotate("rect", xmin = (85/15), ymin = -Inf,
             xmax = (96/15), ymax = Inf,
             fill = c("darkgreen"),colour = 'black', size = 0.5) +
    #annotate("rect", xmin = (97/15), ymin = -Inf,
    #xmax = (124/15), ymax = Inf,
    #fill = c("pink"),colour = 'black', size = 0.5) +
    # annotate("rect", xmin = cut(0-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
    #          xmax = cut(60-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
    #          fill = c("red"),colour = 'black', size = 0.5) +
    # annotate("rect", xmin = cut(60-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
    #          xmax = cut(71-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
    #          fill = c("lightblue"),colour = 'black', size = 0.5) +
    # annotate("rect", xmin = cut(71-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
    #          xmax = cut(83-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
  #          fill = c("blue"),colour = 'black', size = 0.5) +
  # annotate("rect", xmin = cut(83-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
  #          xmax = cut(93-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
  #          fill = c("lightgreen"),colour = 'black', size = 0.5) +
  # annotate("rect", xmin = cut(93-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
  #          xmax = cut(104-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
  #          fill = c("darkgreen"),colour = 'black', size = 0.5) +
  # annotate("rect", xmin = cut(104-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymin = -Inf,
  #          xmax = cut(130-3, breaks = c(-Inf, range, Inf), labels = FALSE), ymax = Inf,
  #          fill = c("pink"),colour = 'black', size = 0.5) +
  theme_void()
  
  
  
  p3 <- p/p2 + plot_layout(heights = c(10, 1))
   
   

 
  ggsave(file.path("ITRRearrangements", paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.png')), p3, dpi = 300, width = 10, units = 'in')
 
  p    
})

#q(save = 'no', status = 0, runLast = FALSE) 
