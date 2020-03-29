TAX, = glob_wildcards("../notrans/{tax}.faa")

rule final:
    input:
        expand("annotations/{tax}.tsv", tax= TAX)
rule interpro:
    input:
        "../notrans/{tax}.faa"
    output:
        "annotations/{tax}.tsv"
    shell:
        "/home-link/bbpgo01/interproscan-5.36-75.0/interproscan.sh -o {output} -goterms -i {input} -f TSV -iprlookup -pa"
        
