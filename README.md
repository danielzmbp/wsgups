# WSGUPS
A pipeline for Whole-genome Screening of Genes Under Positive Selection.
## Instructions
* Download cds_from_genomic and protein files, as well as the .gff
* Use the program proteinortho with the protein and .gff to find the ortholog groups
    * Synteny option on
    * 
* Pool all protein and cds sequences into one file using cat 
    * AA.faa
    * NT.fna
* Rename files so that only the index is the fast name
* Use pillars notebook to make all the family files (make sure the column slicing is not masking any species)
    * families/family_*.fna or faa 
* Run snakemake with “fubarsnake" file to align them, make the tree and calculate fubar.
    * Change the rule hyphen to match the current folder
* Run snakemake with “statfubar” to get the final_results folder with the final stats
* Use fubar notebook to filter the families predicted under positive selection
    * Run through python directly on the server
* Use snakemake "absrelsnake" file to perform absrel
    * Change the folder to match current
* Use absrel notebook to filter families under absrel
    * Run through ipython directly on the server
* Download families and use “absrel_results” notebooks to analyse the results
* Run interproscan on .faa files



