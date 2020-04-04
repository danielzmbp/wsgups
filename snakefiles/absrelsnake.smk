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
        """
        tmpThing=$(find /.. -name aBSREL.bf|tail -n 1);(echo {input.align}; echo {input.tree};|
        hyphy $tmpThing > {output.log} || touch {output.log} {output.json}
        """
