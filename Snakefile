import pandas as pd
import numpy as np
from Bio import SeqIO
import glob
import os
import re
from Bio import codonalign
from Bio import AlignIO
from shutil import copyfile
import phyphy
from ete3 import Tree

G, = glob_wildcards("samples/{g}.faa")

localrules: final, proteinortho, make_families, clean_faa, clean_fna, mafft,
localrules: fasttree, fubar, move_fubar, final_stats, absrel_stats

rule final:
    input: "all_ann.csv"

rule proteinortho:
    input: expand("samples/{g}.faa",g=G)
    output:
        "proteinortho/protein_families.poff.tsv",
        "AA.faa",
        "NT.fna",
    conda: "envs/proteinortho.yaml"
    shell:
        "proteinortho -clean -synteny -singles -nograph -project=protein_families samples/*.faa && "
        "mv protein_families* proteinortho/ && cat samples/*.faa > AA.faa && "
        "cat samples/*.fna > NT.fna"

checkpoint make_families:
    input: "proteinortho/protein_families.poff.tsv"
    output: directory("families/faas")
    script:
        "scripts/pillars.py"

rule clean_fna:
    input:
        "families/fnas/{fam}.fna"
    output:
        temp("families/cleaned_fnas/{fam}.fna.cleaned")
    run:
        new_sequences = []
        for record in SeqIO.parse(input[0], "fasta"):
            if len(record.seq) % 3 != 0:
                record.seq = record.seq + "N"
                new_sequences.append(record)
            else:
                new_sequences.append(record)
        SeqIO.write(new_sequences, output[0], "fasta")

rule clean_faa:
    input:
        "families/faas/{fam}.faa"
    output:
        temp("families/cleaned_faas/{fam}.faa.cleaned")
    run:
        new_sequences = []
        for record in SeqIO.parse(input[0], "fasta"):
            if record.seq.find("Z")!=-1:
                position = record.seq.find("Z")
                record.seq = record.seq.tomutable()
                record.seq[position] = "X"
                new_sequences.append(record)
            elif record.seq.find("J")!=-1:
                position = record.seq.find("J")
                record.seq = record.seq.tomutable()
                record.seq[position] = "X"
                new_sequences.append(record)
            elif record.seq.find("U")!=-1:
                position = record.seq.find("U")
                record.seq = record.seq.tomutable()
                record.seq[position] = "X"
                new_sequences.append(record)
            else:
                new_sequences.append(record)
        SeqIO.write(new_sequences, output[0], "fasta")

rule mafft:
    input:
        "families/cleaned_faas/{fam}.faa.cleaned"
    output:
        temp("families/alns/{fam}.aln")
    conda:
        "envs/mafft.yaml"
    shell:
        "mafft --auto --thread 1 {input} > {output}"
rule fasttree:
    input:
        "families/alns/{fam}.aln"
    output:
        "families/trees/{fam}.tree"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -nosupport {input} > {output} || true"
rule codonaln:
    input:
        pro_align = "families/alns/{fam}.aln",
        nucl_seqs = "families/cleaned_fnas/{fam}.fna.cleaned"
    output:
        alignment = "families/codon_alns/{fam}.aln.codon"
    run:
        aa_aln = AlignIO.read(input.pro_align, "fasta")
        na_seq = SeqIO.to_dict(SeqIO.parse(input.nucl_seqs, "fasta"))

        align = codonalign.build(aa_aln, na_seq, max_score=20)

        for record in range(0, len(align)):
            align._records[record].description = ""  # removes description by looping through the records

        SeqIO.write(align, output.alignment, "fasta")

rule fubar:
    input:
        tree = "families/trees/{fam}.tree",
        align = "families/codon_alns/{fam}.aln.codon"
    output:
        json = temp("families/codon_alns/{fam}.aln.codon.FUBAR.json"),
        log = "families/logs/{fam}.aln.codon.FUBAR.log",
        cache = temp("families/codon_alns/{fam}.aln.codon.FUBAR.cache")
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy fubar --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log} ||"
        "touch {output.log} {output.json} {output.cache}"

def aggregate_fams(wildcards):
    checkpoint_output = checkpoints.make_families.get(**wildcards).output[0]
    return expand("families/logs/{fam}.aln.codon.FUBAR.log",
           fam=glob_wildcards(os.path.join(checkpoint_output, "{fam}.faa")).fam)

rule aggregate_fams:
    input:
        aggregate_fams
    output:
        "final_results/fams_fubar.txt",
    run:
        with open(output[0], "w") as out:
            for currentFile in input:
                with open(currentFile) as f:
                    for line in f:
                        if "## FUBAR" in line:
                            if "no" not in line:
                                result = re.search('inferred(.*)sites', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n") # family num_selected_sites

checkpoint move_fubar:
    input:
        "final_results/fams_fubar.txt"
    output:
        directory("families_fubar")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=None)
        families = fams.iloc[:,0]
        families_in_dir = glob.glob("families/**/*")
        if not os.path.exists('families_fubar/fnas'):
            os.makedirs('families_fubar/fnas')
        if not os.path.exists('families_fubar/faas'):
            os.makedirs('families_fubar/faas')
        if not os.path.exists('families_fubar/logs'):
            os.makedirs('families_fubar/logs')
        if not os.path.exists('families_fubar/trees'):
            os.makedirs('families_fubar/trees')
        if not os.path.exists('families_fubar/codon_alns'):
            os.makedirs('families_fubar/codon_alns')
        for i in range(0,len(families_in_dir)):
            if int(families_in_dir[i].split(".")[0].split("/")[-1]) in list(families):
                copyfile(families_in_dir[i], "families_fubar/"+families_in_dir[i].split("/",1)[1])

rule absrel:
    input:
        tree = "families_fubar/trees/{i}.tree",
        align = "families_fubar/codon_alns/{i}.aln.codon"
    output:
        json = "families_fubar/codon_alns/{i}.aln.codon.ABSREL.json",
        log = "families_fubar/logs/{i}.aln.codon.ABSREL.log"
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy absrel --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log}"


def aggregate_fubar(wildcards):
    checkpoint_output = checkpoints.move_fubar.get(**wildcards).output[0]
    return expand("families_fubar/logs/{i}.aln.codon.ABSREL.log",
           i=glob_wildcards(os.path.join(checkpoint_output, "trees/{i}.tree")).i)

rule absrel_stats:
    input:
        aggregate_fubar
    output:
        "final_results/fams_absrel.txt"
    run:
        with open(output[0], "w") as out:
            for currentFile in input:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0** branches" not in line:
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + "\n")

checkpoint move_absrel:
    input:
        "final_results/fams_absrel.txt"
    output:
        directory("families_absrel")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=None)
        families = fams[0]
        families_in_dir = glob.glob("families_fubar/**/*")
        if not os.path.exists('families_absrel/fnas'):
            os.makedirs('families_absrel/fnas')
        if not os.path.exists('families_absrel/faas'):
            os.makedirs('families_absrel/faas')
        if not os.path.exists('families_absrel/logs'):
            os.makedirs('families_absrel/logs')
        if not os.path.exists('families_absrel/trees'):
            os.makedirs('families_absrel/trees')
        if not os.path.exists('families_absrel/codon_alns'):
            os.makedirs('families_absrel/codon_alns')
        for i in range(0,len(families_in_dir)):
            if int(families_in_dir[i].split(".")[0].split("/")[-1]) in list(families):
                copyfile(families_in_dir[i], "families_absrel/"+families_in_dir[i].split("/",1)[1])

def aggregate_absrel(wildcards):
    checkpoint_output = checkpoints.move_absrel.get(**wildcards).output[0]
    return expand("families_absrel/logs/{j}.aln.codon.ABSREL.log",
           j=glob_wildcards(os.path.join(checkpoint_output, "trees/{j}.tree")).j)

rule final_stats:
    input:
        aggregate_absrel
    output:
        "gene_species_assoc.csv",
        "go_mapping",
        "all_ann.csv"
    run:
        absrel = glob.glob("families_absrel/logs/*.ABSREL.log")

        family_list = []
        branch_list = []
        pvalue_list = []

        for file in absrel:
            with open(file) as myfile:
                for line in myfile:
                    if re.search(r'^\* \w.+ p-value',line):
                        family, branch, pvalue = file[21:].split(".")[0], line[2:].split(",")[0],line.split()[-1]
                        family_list.append(family)
                        branch_list.append(branch)
                        pvalue_list.append(pvalue)

        ps = pd.DataFrame({"family": family_list,'branch': branch_list, 'p-value': pvalue_list})

        absrel_json = glob.glob("families_absrel/codon_alns/*.ABSREL.json")

        tree_list=[]
        family_list=[]

        for file in absrel_json:
            if os.stat(file).st_size > 0:
                tree,family = phyphy.Extractor(file).extract_input_tree(),file[27:].split(".")[0]
                tree_list.append(tree)
                family_list.append(family)  # the problem is that there are empty json files and the script files

        tdf = pd.DataFrame({"tree":tree_list,"family":family_list})

        ps = ps.merge(tdf)

        children_list=[]

        for i in range(0,len(ps)):
            t = Tree(ps.tree[i], format=1)
            node = t.search_nodes(name=ps.branch[i])[0]
            children_list.append(node.get_leaf_names())

        ps["children"] = children_list

        lst_col = 'children'

        r = pd.DataFrame({
              col:np.repeat(ps[col].values, ps[lst_col].str.len())
              for col in ps.columns.drop(lst_col)}
            ).assign(**{lst_col:np.concatenate(ps[lst_col].values)})[ps.columns]

        r_dd = r.drop_duplicates(subset=['family', 'children'])

        files = os.listdir("samples/")

        taxa=[]
        for file in files:
            if file.endswith(".faa"):
                taxa.append(file.split(".")[0])
        empty_dict={}
        for t in taxa:
            test = SeqIO.to_dict(SeqIO.parse("samples/"+t+".faa","fasta"))
            for key,value in test.items():
                test[key] = t
                empty_dict.update(test)

        spec_name= pd.DataFrame.from_dict(empty_dict,orient='index',
               columns=["species"])

        spec_name["children"] = spec_name.index

        spec_name.replace("\.","_",regex=True,inplace=True)

        spec_name.to_csv("gene_species_assoc.csv",index=None)

        proteome = glob.glob("samples/*.tsv")

        col_names = ["protein_accession","md5","length","analysis", "signature_accession",
                     "signature_description", "start","stop","score","status","date",
                     "ip_accesion","ip_description","go","pathway"]

        annotations = []

        for infile in proteome:
            data = pd.read_csv(infile, sep='\t',names=col_names)
            data["family"] = infile[8:].split(".")[0]
            annotations.append(data)

        all_annotations = pd.concat(annotations).reset_index()

        all_annotations["protein_accession"] = all_annotations["protein_accession"].replace("\.","_",regex=True)

        r_dd.to_csv("r_dd.csv")

        all_annotations["ps"] = np.where(all_annotations["protein_accession"].isin(r_dd["children"]), 1, 0)

        families = pd.read_csv("fam.txt","\s+",header=None)

        families = families[1].str.split(",",expand=True).fillna(value=np.nan, inplace=True)

        families = families.dropna(subset=[1]).melt().dropna()

        dups = families["value"]

        all_annotations["dups"] = np.where(all_annotations["protein_accession"].isin(dups), 1, 0)

        # Make GO associations file

        all_annotations["go"] = all_annotations["go"].str.replace("|",",")

        all_annotations[["protein_accession","go"]].dropna().to_csv("go_mapping.csv", sep="\t", index= False)

        # Write annotations to file

        all_annotations.to_csv("all_ann.csv")