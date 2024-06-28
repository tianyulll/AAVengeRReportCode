library(aavAnalysis)
library(argparse)

# parameters
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "input data file in rds/excel, or name in database")
parser$add_argument("-o" ,"--outputDir", default = "output", help = "Output directory")
parser$add_argument("-t",  "--reportTitle", default = "AAV_report", help = "file name for the report")
parser$add_argument("--piNote", help = "path to text file for summary notes")
parser$add_argument("-m", "--meta", help = "path to meta data table, must have column sample and info")
parser$add_argument("--species", default = "human", required = F, help = "choose from human or mice")
parser$add_argument("-f", "--filter", default = 0, required = F,  help = "only keep sites with reads greater than this threshold")
parser$add_argument("--saveimg", action = 'store_true', help = "if provided, plots will be saved separately as png")
parser$add_argument("--savetsv", action = 'store_true', help = "if provided, tables will be saved separately as tsv")


# these parameters will be removed in future versions
parser$add_argument("--itrStart", type="integer", default = 57, help = "itr seq start position for remnant plot")
parser$add_argument("--itrLength", type="integer", default = 197, help = "itr length for remnant plots")
parser$add_argument("--ntBinSize", type="integer", default = 3, help = "bin size for remnant plot")

args <- parser$parse_args()

# to be removed in future versions
buildAAVremnantPlots_ITRlength <-  args$itrLength
buildAAVremnantPlots_ITRseqStart <- args$itrStart
buildAAVremnantPlots_NTbinSize <- args$ntBinSize

# Handle zero input
if (is.null(args$input)) {
  message("Input not given; pulling available info from db")
  aavAnalysis::getAvailAAVengeR()
  stop("Please run pipeline with one of the trials")
}

# read input as df
df <- aavAnalysis::getAAVengerData(args$input)

# Process input
df <- aavAnalysis::getProcessedDf(df, meta = args$meta, minreads = args$filter)

# Run a summary
df.summary <- aavAnalysis::getSummary(df)

message("summary created!")

# Abundance analysis
df.abund <- aavAnalysis::getAbundanceDf(df)
abundancePlots <- aavAnalysis::plotListAbundance(df.abund,
                                                 saveimg = args$saveimg, outputDir = args$o)

message("Abundance plots created!")

# Rearrangment summary
df.rearrangement <- aavAnalysis::getRearrangeDf(df)

rearrangePlot <- aavAnalysis::plotRearrangment(df.rearrangement)

# ITR remnant plots
remnantPlots <- aavAnalysis::plotListRemnant(df = df, saveimg = args$saveimg, outDir = args$o,
                                             buildAAVremnantPlots_ITRseqStart = buildAAVremnantPlots_ITRseqStart,
                                             buildAAVremnantPlots_ITRlength = buildAAVremnantPlots_ITRlength,
                                             buildAAVremnantPlots_NTbinSize = buildAAVremnantPlots_NTbinSize)

message("remnant plots created!")


# Create Gene distribution df

## Find random Value
if (args$species == "human") {
  message("Using human as default species...")
  dash <- readRDS(file = "reference/hg38.dash.rds")
} else if (args$species == "mice") {
  dash <- readRDS(file = "reference/mm9.dash.rds")
} else (
  stop("provided species is not included in reference")
)

gene.df <- aavAnalysis::getGeneDf(df)
gene.plot <- aavAnalysis::plotRandomGene(df = gene.df,
                                         random_exon = dash[[1]], random_tu = dash[[2]])

## saving tables separately?
if (args$savetsv) {

  message("saving tables...")

  tableDir = file.path(args$outputDir, "tables")
  if (!file.exists(tableDir)) {
    dir.create(tableDir)
    message(paste("Directory created", tableDir))
  }

  write.table(summary, sep = '\t', row.names = F, quote = F,
              file = file.path(tableDir, "summaryTable.tsv"))
  write.table(rearrangement, sep = '\t', row.names = F, quote = F,
              file = file.path(tableDir, "rearrangementTable.tsv"))
}


# Create report
message("knitting reports...")

suppressMessages(rmarkdown::render("AAV_report.Rmd",
                  output_dir = args$outputDir, output_file = paste0(args$reportTitle, ".pdf"),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = args$reportTitle)))

message("done")
