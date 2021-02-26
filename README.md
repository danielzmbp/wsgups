# WSGUPS
A snakemake pipeline for Whole-genome Screening of Genes Under Positive Selection (WSGUPS).

## Requirements
* Conda
* Snakemake > 5.0

## Instructions
* Download coding DNA (extension should be .fna) and protein files (extension should be .faa) and place 
them in the folder samples
* Add InterProScan annotated *.tsv files to samples folder with the same name for each taxa.
* Include .gff3 files in the samples folder with the same name for each taxa.
* Make sure that .faa, .fna, .tsv and .gff3 files have the same name for the same species 
* All fasta headers should be unique and have no spaces
* Run snakemake from the wsgups folder with the command:
`snakemake --snakefile Snakefile --use-conda -j <cpu cores to use> -k`
* Use absrel_results notebook to analyze the results

## Output files
* `all_ann.tsv`: the main result file. It's tab separated table which include all annotations with annotated positive selection resuls
* `gene_species_assoc.csv`: gene and species in csv format.
* `go_mapping.csv`: GO associations file.

## License
This pipeline is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).

## Preprint
Please find the preprint [here](https://www.biorxiv.org/content/10.1101/2021.01.12.426341v1).
