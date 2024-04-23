# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. 
`reference` folder contains required reference data for calculating random integration site distribution. \
`sampleReport` contains a sample output by running the toolkit on tak981 dataset. \
`createRandomPercentage.r` is not required. It was used to build reference dataset.

Usage: `Rscript createAAVReport.R -i DeJongMice.rds -o test -r test1`

Parameters: 

`-i --input:	expect RDS input (required)` \
`-o --outputDir:	directory to write outputs ` \
`-r --reportTitle:	report title and output file title` \
`--itrStart:	itr seq start position for remnant plot` \
`--itrLength:	itr length for remnant plots` \
`--piNote:  path to notes from pi`

`AAV_report.rmd` The default report rmd file. \
`libraries.tsv` Main libraries required for this repo. \
`full_depend.tsv` Full dependencies info.

A few assumptions of the toolkit: the input file must be AAVengeR outputs in RDS. The meta-must have column sample and info in RDS. All other parameters are optional. By default the outputs will write into the current working directory.
