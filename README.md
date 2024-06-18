# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. \
Default report contains summary tables, remnant plots of integration sites, abundance plots of clones, in exon and in gene percentage plots.

### Overview
`reference` folder contains required reference data for calculating random integration site distribution. \
`sampleReport` contains a sample output by running the toolkit. \
`createRandom**.r` is not required. It was used to build randomization reference for species.
`AAV_report.rmd` The default report rmd file. \
`libraries.tsv` Main libraries required for this repo. 

## Usage
1. `git clone` this repository. 
2. Install libraries required in "libraries.tsv"
3. `cd` into this repo
4. `Rscript createAAVReport.R -i aavenger_output.rds -o outputDir -r reportTitle -s Species`

example: `Rscript createAAVReport.R -i AAVHelaTopo -o sampleReport -t topoDefault --meta topoisomerase_meta.rds`

### Parameters: 

`-i --input:	accepts rds/excel, or project name in database` \
`-o --outputDir:	directory to write outputs ` \
`-t --reportTitle:	report title and output file title` \
`--piNote:  path to notes from pi` \
`-m --meta:  meta data of samples` \
`-s --species:  species of the data ` \
`-f --filters: a threshold for filtering site based on reads` \
`--itrStart:	itr seq start position for remnant plot` \
`--itrLength:	itr length for remnant plots` 

Default behaviors of the software: 

It will write temporary files in working directory. All plots and reports will write in the output directory. Only `input` is required. However it is recommended to provide more parameters. The default species is human. Currently we only support human and mouse.m 

Meta data can be provided as a RDS file containing sample and info columns. sample columns must correspond to GTSP ids. If not provided, the subject column will be used as meta data information.

The input can be a local rds or excel file from AAVengeR outputs. 
It also accepts a string of the project name in AAVengeR database.

## Contact
[Tianyu](mailto:tianyu.lu@pennmedicine.upenn.edu)
[Aradhana]()
