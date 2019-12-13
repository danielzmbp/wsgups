[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/danielzmbp/wsgups/master)

# WSGUPS
A pipeline for Whole-genome Screening of Genes Under Positive Selection.
## Instructions
* Download cds_from_genomic and protein files, as well as the .gff
* Use the program proteinortho with the protein and .gff to find the ortholog groups
    * Synteny option on and diamond blast algorithm
    * `proteinortho -project=fungi -cpus=2 -p=diamond -synteny .`
* Pool all protein and cds sequences into one file using cat
    * AA.faa
    * NT.fna
* Rename files so that only the index is the fast name
* Use pillars notebook to make all the family files (make sure the column slicing is not masking any species)
    * families/family_*.fna or faa
* Run snakemake with “fubarsnake" file to align them, make the tree and calculate fubar.
* Run snakemake with “statfubar” to get the final_results folder with the final stats
* Run snakemake with "absrelsnake" file to perform absrel
* Run snakemake with "statabsrel" to get the final statistics
* Download families and use “absrel_results” notebooks to analyse the results
* Run interproscan on .faa files
    * Use `interproscansnake file`
* Run GO enrichment
    * Using go_BP notebook
