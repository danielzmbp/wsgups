import pandas as pd
from Bio import SeqIO
import glob
import os
import re
from shutil import copyfile
from snakemake.utils import min_version

min_version("5")

(G,) = glob_wildcards("samples/{g}.faa")


localrules:
    final,
    proteinortho,
    make_families,
    clean_fna,
    mafft,
    codonaln,
    aggregate_fams,
    fasttree,
    fubar,
    move_fubar,
    final_stats,
    absrel_stats,
    absrel,


rule final:
    input:
        "final_results/all_ann.csv",


rule proteinortho:
    input:
        expand("samples/{g}.faa", g=G),
    output:
        "proteinortho/protein_families.poff.tsv",
        "AA.faa",
        "NT.fna",
    conda:
        "envs/proteinortho.yaml"
    shell:
        "proteinortho -clean -synteny -singles -nograph -project=protein_families samples/*.faa && "
        "mv protein_families* proteinortho/ && cat samples/*.faa > AA.faa && "
        "cat samples/*.fna > NT.fna"


checkpoint make_families:
    input:
        "proteinortho/protein_families.poff.tsv",
    output:
        directory("families/faas"),
    params:
        cutoff=5,  # change for minimum family size to analyze
    script:
        "scripts/pillars.py"


rule clean_fna:
    input:
        "families/fnas/{fam}.fna",
    output:
        temp("families/cleaned_fnas/{fam}.fna.cleaned"),
    run:
        new_sequences = []
        for record in SeqIO.parse(input[0], "fasta"):
            if len(record.seq) % 3 != 0:
                record.seq += "N"
                new_sequences.append(record)
            else:
                new_sequences.append(record)
        SeqIO.write(new_sequences, output[0], "fasta")


rule mafft:
    input:
        "families/faas/{fam}.faa",
    output:
        temp("families/alns/{fam}.aln"),
    conda:
        "envs/mafft.yaml"
    shell:
        "mafft --auto --thread 1 {input} > {output}"


rule fasttree:
    input:
        "families/alns/{fam}.aln",
    output:
        "families/trees/{fam}.tree",
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -nosupport {input} > {output} || true"


rule codonaln:
    input:
        pro_align="families/alns/{fam}.aln",
        nucl_seq="families/cleaned_fnas/{fam}.fna.cleaned",
    output:
        alignment="families/codon_alns/{fam}.aln.codon",
    conda:
        "envs/pal2nal.yaml"
    shell:
        "pal2nal.pl {input.pro_align} {input.nucl_seq} -output fasta > {output.alignment}"


rule fubar:
    input:
        tree="families/trees/{fam}.tree",
        align="families/codon_alns/{fam}.aln.codon",
    output:
        json=temp("families/codon_alns/{fam}.aln.codon.FUBAR.json"),
        log="families/logs/{fam}.aln.codon.FUBAR.log",
        cache=temp("families/codon_alns/{fam}.aln.codon.FUBAR.cache"),
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy fubar --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log} ||"
        "touch {output.log} {output.json} {output.cache}"


def aggregate_fams(wildcards):
    checkpoint_output = checkpoints.make_families.get(**wildcards).output[0]
    return expand(
        "families/logs/{fam}.aln.codon.FUBAR.log",
        fam=glob_wildcards(os.path.join(checkpoint_output, "{fam}.faa")).fam,
    )


rule aggregate_fams:
    input:
        aggregate_fams,
    output:
        "final_results/fams_fubar.txt",
    run:
        with open(output[0], "w") as out:
            for currentFile in input:
                with open(currentFile) as f:
                    for line in f:
                        if "## FUBAR" in line:
                            if "no" not in line:
                                result = re.search("inferred(.*)sites", line)
                                out.write(
                                    currentFile.split("/")[-1].split(".")[0]
                                    + " "
                                    + result.group(1)
                                    + "\n"
                                )  # family num_selected_sites


checkpoint move_fubar:
    input:
        "final_results/fams_fubar.txt",
    output:
        directory("families_fubar"),
    run:
        families = pd.read_csv(input[0], sep="\s+", index_col=False, header=None).iloc[:, 0]
        families_in_dir = glob.glob("families/**/*")
        dirs = ["fnas", "faas", "logs", "trees", "codon_alns", "alns", "cleaned_fnas"]

        for d in dirs:
            if not os.path.exists("families_fubar/" + d):
                os.makedirs("families_fubar/" + d)

        for i in range(0, len(families_in_dir)):
            if int(families_in_dir[i].split(".")[0].split("/")[-1]) in list(families):
                copyfile(
                    families_in_dir[i],
                    "families_fubar/" + families_in_dir[i].split("/", 1)[1],
                )


rule absrel:
    input:
        tree="families_fubar/trees/{i}.tree",
        align="families_fubar/codon_alns/{i}.aln.codon",
    output:
        json="families_fubar/codon_alns/{i}.aln.codon.ABSREL.json",
        log="families_fubar/logs/{i}.aln.codon.ABSREL.log",
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy absrel --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log}"


def aggregate_fubar(wildcards):
    checkpoint_output = checkpoints.move_fubar.get(**wildcards).output[0]
    return expand(
        "families_fubar/logs/{i}.aln.codon.ABSREL.log",
        i=glob_wildcards(os.path.join(checkpoint_output, "trees/{i}.tree")).i,
    )


rule absrel_stats:
    input:
        aggregate_fubar,
    output:
        "final_results/fams_absrel.txt",
    run:
        with open(output[0], "w") as out:
            for currentFile in input:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0** branches" not in line:
                                out.write(
                                    currentFile.split("/")[-1].split(".")[0] + "\n"
                                )


checkpoint move_absrel:
    input:
        "final_results/fams_absrel.txt",
    output:
        directory("families_absrel"),
    run:
        families = pd.read_csv(input[0], sep="\s+", index_col=False, header=None)[0]
        families_in_dir = glob.glob("families_fubar/**/*")
        dirs = ["fnas", "faas", "logs", "trees", "codon_alns", "alns", "cleaned_fnas"]

        for d in dirs:
            if not os.path.exists("families_absrel/" + d):
                os.makedirs("families_absrel/" + d)

        for i in range(0, len(families_in_dir)):
            if int(families_in_dir[i].split(".")[0].split("/")[-1]) in list(families):
                copyfile(
                    families_in_dir[i],
                    "families_absrel/" + families_in_dir[i].split("/", 1)[1],
                )


def aggregate_absrel(wildcards):
    checkpoint_output = checkpoints.move_absrel.get(**wildcards).output[0]
    return expand(
        "families_absrel/logs/{j}.aln.codon.ABSREL.log",
        j=glob_wildcards(os.path.join(checkpoint_output, "trees/{j}.tree")).j,
    )


rule final_stats:
    input:
        aggregate_absrel,
    output:
        "final_results/gene_species_assoc.csv",
        "final_results/go_mapping.csv",
        "final_results/all_ann.csv",
    conda:
        "envs/final_stats.yaml"
    script:
        "scripts/final_stats.py"
