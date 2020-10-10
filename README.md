[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/danielzmbp/wsgups/master)

# WSGUPS
A pipeline for Whole-genome Screening of Genes Under Positive Selection.
## Requirements
* Conda
* Snakemake
## Instructions
* Download cds_from_genomic and protein files, as well as the genomic .gff and place them in the folder samples
* Make sure all files have congruent file names
* Run snakemake with the command:
`snakemake --snakefile Snakefile --use-conda -j <cpu cores to use>`
* Use absrel_results notebook to analyze the results
* Run interproscan on .faa files
    * Use `interproscansnake` file
* Run GO enrichment
    * Using go_BP notebook
