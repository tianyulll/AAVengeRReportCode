library(dplyr)
library(ggplot2)
library(reshape2)
library(writexl)

sites <- read_rds("AnchorReadMaps.rds")

#Making Blank DataFrame#

df <- data.frame(matrix(ncol = 9, nrow = 5))
colnames(df) <- c("Regions", "GTSP5787","GTSP5792","GTSP5790","GTSP5795","GTSP5786","GTSP5791","GTSP5789","GTSP5794")
df$Regions <- c("Red","Blue1","Blue2","Green1","Green2")
samples <- c("GTSP5787","GTSP5792","GTSP5790","GTSP5795","GTSP5786","GTSP5791","GTSP5789","GTSP5794")
#samples <- c("GTSP5786")


#Running Through Each of the Samples#
for (i in samples) {
x <- sites %>% filter(sample == i)
range <- seq(0, 130, 3)
b <- nrow(x)
d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$rearrangement, '\\.\\.\\d+'))) + 34,
            bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (3/2),
            r = stringr::str_count(x$rearrangement, ';') -1,
            r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n(),per = ((n/b)) * 100) %>% ungroup() %>% mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', '≥ 5')))) 


##Data Frames Made to Get Percentages##
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

##Seperating Based on ITR Breakpoint Region##
e1 <- e %>% filter(e$bin < 21)
e2 <- e %>% filter(e$bin > 21 & e$bin < 25)
e3 <- e %>% filter(e$bin > 25 & e$bin < 28)
e4 <- e %>% filter(e$bin > 28 & e$bin < 32) 
e5 <- e %>% filter(e$bin > 32)

##Combining Percents and Adding to DataFrame##

percents <- c(sum(e1$per2),sum(e2$per2),sum(e3$per2),sum(e4$per2),sum(e5$per2))
df[i] <- percents

}


groups <- c("", "DMSOCrude 1","DMSOCrude 2","HDRiCrude 1","HDRiCrude 2","DMSOPure 1","DMSOPure 2","HDRiPure 1","HDRiPure 2")
df <- rbind(groups,df)
write_xlsx(df, "ITRRearrangements/percents_output.xlsx")
