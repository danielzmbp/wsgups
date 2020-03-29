import os

wd = os.getcwd()

# hyphy_dir=  "/home-link/bbpgo01/.conda/envs/python3/lib/hyphy"

FAM, = glob_wildcards("families_fubar/family_{fam}.faa")

rule final:
    input:
        expand("families_fubar/family_{fam}.aln.codon.ABSREL.json", fam=FAM)

rule hyphy:
    input:
        tree = "families_fubar/family_{fam}.tree",
        align = "families_fubar/family_{fam}.aln.codon"
    output:
        json = "families_fubar/family_{fam}.aln.codon.ABSREL.json",
        log = "families_fubar/family_{fam}.aln.codon.ABSREL.log"
    conda:
        "envs/hyphy.yaml"
    shell:
        "tmpThing=$(find " + wd +
        " -name aBSREL.bf|tail -n 1);(echo " + wd +
        "{input.align}; echo " + wd +
        "{input.tree})|hyphy $tmpThing > {output.log} || touch {output.log} {output.json} "
