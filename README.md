# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. 

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
