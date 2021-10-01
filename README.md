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
* `all_ann.tsv`: the main result file. It's a tab separated table which include all annotations with annotated positive selection resuls
* `gene_species_assoc.csv`: gene and species in csv format.
* `go_mapping.csv`: GO associations file.

## License
This pipeline is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).

## Publication
If you find this useful please consider citing our publication found [here](https://www.mdpi.com/2076-0817/10/7/807).

Gómez-Pérez, D.; Kemen, E. Predicting Lifestyle from Positive Selection Data and Genome Properties in Oomycetes. *Pathogens* **2021**, 10, 807. https://doi.org/10.3390/pathogens10070807 

`````
@Article{pathogens10070807,
AUTHOR = {Gómez-Pérez, Daniel and Kemen, Eric},
TITLE = {Predicting Lifestyle from Positive Selection Data and Genome Properties in Oomycetes},
JOURNAL = {Pathogens},
VOLUME = {10},
YEAR = {2021},
NUMBER = {7},
ARTICLE-NUMBER = {807},
URL = {https://www.mdpi.com/2076-0817/10/7/807},
PubMedID = {34202069},
ISSN = {2076-0817},
ABSTRACT = {As evidenced in parasitism, host and niche shifts are a source of genomic and phenotypic diversification. Exemplary is a reduction in the core metabolism as parasites adapt to a particular host, while the accessory genome often maintains a high degree of diversification. However, selective pressures acting on the genome of organisms that have undergone recent lifestyle or host changes have not been fully investigated. Here, we developed a comparative genomics approach to study underlying adaptive trends in oomycetes, a eukaryotic phylum with a wide and diverse range of economically important plant and animal parasitic lifestyles. Our analysis reveals converging evolution on biological processes for oomycetes that have similar lifestyles. Moreover, we find that certain functions, in particular carbohydrate metabolism, transport, and signaling, are important for host and environmental adaptation in oomycetes. Given the high correlation between lifestyle and genome properties in our oomycete dataset, together with the known convergent evolution of fungal and oomycete genomes, we developed a model that predicts plant pathogenic lifestyles with high accuracy based on functional annotations. These insights into how selective pressures correlate with lifestyle may be crucial to better understand host/lifestyle shifts and their impact on the genome.},
DOI = {10.3390/pathogens10070807}
}
