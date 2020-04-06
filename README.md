[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/danielzmbp/wsgups/master)

# WSGUPS
A pipeline for Whole-genome Screening of Genes Under Positive Selection.
## Instructions
* Download cds_from_genomic and protein files, as well as the genomic .gff  
* Use the program proteinortho with the protein and .gff to find the ortholog groups
    * Synteny option on and diamond blast algorithm
    * `proteinortho -project=<project-name> -cpus=<cpu-threads> -p=diamond -synteny *.faa`
* Rename files so that only the index is the fasta header
* Pool all protein and cds sequences into one file using cat
    * `cat *.faa > AA.faa`
    * `cat *.fna > NT.fna`
* Use "pillars.py" script to make all the family files (make sure the column slicing is not masking any species)
    * families/family_*.fna or faa
* Perform **FUBAR** analysis
    * Run snakemake with “fubarsnake.smk" file to align them, make the tree and calculate fubar.
    * Run script "fubar.py" in scripts folder in order to move positive reults to a new folder called families_fubar
* Perform **aBSREL** analysis
    * Run snakemake with "absrelsnake.smk" file to perform absrel
    * Run script "absrel.py" in scripts folder in order to move positive results to a new folder called families_absrel
* Use absrel_results notebooks to analyze the results
* Run interproscan on .faa files
    * Use `interproscansnake` file
* Run GO enrichment
    * Using go_BP notebook
