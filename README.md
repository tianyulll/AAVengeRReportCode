# AAVengeRReportCode

This tool kit generates a standarized, reproducible report from AAVengeR outputs. \
Default report contains summary tables, remnant plots of integration sites, abundance plots of clones, in exon and in gene percentage plots.

### Overview
`reference` folder contains required reference data for calculating random integration site distribution. \
`sampleReport` contains a sample output by running the toolkit. \
`script`  It was used to build randomization reference for different species. \
`AAV_report.rmd` The default report format. \
`libraries.tsv` Main libraries required for this repo. 

## Usage
1. `git clone` this repository. 
2. Install libraries required in "libraries.tsv"
3. `cd` into this repo
4. `sh init.sh` will initialize credential information to access database.
5. Runn script: `Rscript createAAVReport.R -i aavenger_output.rds -o outputDir -r reportTitle -s Species`. 
You can find a full list of available parameters below.
6. To get a list of available aavenger report, you can run the script with no inputs:   ``Rscript createAAVReport.R`

example: `Rscript createAAVReport.R -i AAVHelaTopo -o sampleReport -t topoDefault --meta topoisomerase_meta.rds`

### Parameters: 

`-i --input:	accepts rds/excel, or project name in database` \
`-o --outputDir:	directory to write outputs ` \
`-t --reportTitle:	report title and output file name` \
`--piNote:  path to notes from pi` \
`-m --meta:  meta data contain column sample and info in rds` \
`--species:  species of the samples ` \
`-f --filters: a threshold for filtering site based on reads` \
`--saveimg: if flagged, plots will be saved as png `


Behavior of the pipeline: 

It will write temporary files in working directory. All plots and reports will write in the output directory. `input` is required to run the pipeline. However it is recommended to provide more parameters. The default species is human. Currently we only support human and mouse. 

Meta data can be provided as a RDS file containing sample and info columns. sample columns must correspond to GTSP ids. If not provided, the subject column will be used as meta data information.

The input can be a local rds or excel file from AAVengeR outputs. 
It also accepts a string of the project name in AAVengeR database. If you don't know the name of the project, a dry run of the script will search for existing trial names and print.

## Contact

[Tianyu](mailto:tianyu.lu@pennmedicine.upenn.edu) \
[Aradhana]()
