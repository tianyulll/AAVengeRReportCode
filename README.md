# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. 

Usage: `Rscript createAAVReport.R -i DeJongMice.rds -o test -r test1`

Parameters: 

`-i --input:	expect RDS input (required)` \
`-o --outputDir:	directory to write outputs ` \
`-r --reportTitle:	report title and output file title` \
`--itrStart:	itr seq start position for remnant plot` \
`--itrLength:	itr length for remnant plots`

`AAV_report.rmd` The default report rmd file. \
`dependencies.tsv` A list of dependent libraries.



The SonjaBaselines are files that were created by Sonja

In the MaurerBerkleyProject, Report_Maker1.1AK and AAV_Baseline.Rmd are the first versions updated by Aradhana
