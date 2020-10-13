# WSGUPS
A pipeline for Whole-genome Screening of Genes Under Positive Selection.

## Requirements
* Conda
* Snakemake > 5.0

## Instructions
* Download coding DNA (extension should be .fna) and protein files (extension should be .faa) and place them in the folder samples
* Make sure .faa and .fna files have the same name for the same species 
* All fasta headers should be unique and have no spaces
* Run snakemake from the wsgups folder with the command:
`snakemake --snakefile Snakefile --use-conda -j <cpu cores to use> -k`
* Use absrel_results notebook to analyze the results

### Optional
* Optionally, to run proteinortho using the -synteny parameter to get more accurate orthogroups add .gff3 files to the samples folder with the same for each taxa.
* To do the GO enrichment you need the annotated proteins in .tsv in the interproscan output format
    * To run de novo analysis use the interproscansnake file
* Run GO enrichment
    * Using go_BP notebook
    
## License
This module is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).

