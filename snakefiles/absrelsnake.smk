FAM, = glob_wildcards("families_fubar/family_{fam}.faa")

rule final:
    input:
        expand("families_fubar/codon_alns/{fam}.aln.codon.ABSREL.json", fam=FAM)

rule hyphy:
    input:
        tree = "families_fubar/trees/{fam}.tree",
        align = "families_fubar/codon_alns/{fam}.aln.codon"
    output:
        json = "families_fubar/codon_alns{fam}.aln.codon.ABSREL.json",
        log = "families_fubar/logs/{fam}.aln.codon.ABSREL.log"
    conda:
        "envs/hyphy.yaml"
    shell:
        """
        hyphy absrel --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log} || touch {output.log} {output.json}
        """
