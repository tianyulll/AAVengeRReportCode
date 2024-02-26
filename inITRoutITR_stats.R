library(dplyr)
library(data.table)
library(ggplot2)
library(rcompanion)
library(ggmosaic)

#Delineate which GTSPs are in which experimental group
G1 <- c(5786,5791)
G2 <- c(5787,5792)
G3 <- c(5788,5793)
G4 <- c(5789,5794)
G5 <- c(5790,5795)
#G4 <- c("5331,5332,5333,5334")

# name the groups
nG1 <- "HDRiCrude"
nG2 <- "DMSOCrude"
nG3 <- "Negative"
nG4 <- "HDRiPure"
nG5 <- "DMSOPure"
#nG4 <- "Control"

# ITR width: 96
# Position of vector in vector plasmid.
vectorStart <- 1346
vectorEnd <- 3404
ITRwidth <- 96
r <- readRDS('AnchorReadMaps.rds')
r$sample2<-gsub("GTSP","",as.character(r$sample))
#Delineate which GTSPs are in which experimental group
r$Drug <- ""
r$Drug[r$sample2 %in% G1] <- nG1
r$Drug[r$sample2 %in%  G2] <- nG2
r$Drug[r$sample2 %in%  G3] <- nG3
r$Drug[r$sample2 %in%  G4] <- nG4
r$Drug[r$sample2 %in%  G5] <- nG5


t <- ggplot(r, aes(pos)) +
  theme_bw() +
  geom_density() +
  geom_vline(xintercept = c(0, ITRwidth, vectorEnd-vectorStart-ITRwidth, vectorEnd-vectorStart), color = 'red') +
  facet_wrap(~Drug, ncol = 2, scales = 'free_y')

ggsave(filename = 'distribution.png', plot = t,width = 8, height = 6, path = "Output", device = 'png')


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


