# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. 
`reference` folder contains required reference data for calculating random integration site distribution. \
`sampleReport` contains a sample output by running the toolkit. \
`createRandom**.r` is not required. It was used to build randomization reference for species.

Usage: `Rscript createAAVReport.R -i aavenger_output.rds -o outputDir -r reportTitle -s Species`

Parameters: 

`-i --input:	expect RDS input (required)` \
`-o --outputDir:	directory to write outputs ` \
`-r --reportTitle:	report title and output file title` \
`--itrStart:	itr seq start position for remnant plot` \
`--itrLength:	itr length for remnant plots` \
`--piNote:  path to notes from pi` \
`--meta:  meta data of samples` \
`-s:  species of the data `

`AAV_report.rmd` The default report rmd file. \
`libraries.tsv` Main libraries required for this repo. \

A few assumptions of the toolkit: the input file must be AAVengeR outputs in RDS. The meta-must have column sample and info in RDS. All other parameters are optional. By default the outputs will write into the current working directory.
