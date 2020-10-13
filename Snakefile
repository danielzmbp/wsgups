import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import glob

from Bio import codonalign
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
from shutil import copyfile

import os

G, = glob_wildcards("samples/{g}.faa")

rule final:
    input: "final_results/fams_absrel.txt"

rule proteinortho:
    input: expand("samples/{g}.faa",g=G)
    output:
        "proteinortho/protein_families.proteinortho.tsv",
        "AA.faa",
        "NT.fna",
    conda: "envs/proteinortho.yaml"
    shell:
        "proteinortho -clean -project=protein_families samples/*.faa && "
        "mv protein_families* proteinortho/ && cat samples/*.faa > AA.faa && "
        "cat samples/*.fna > NT.fna"

checkpoint make_families:
    input: "proteinortho/protein_families.proteinortho.tsv"
    output: directory("families/faas")
    run:
        poff_tsv = input[0]  # replace with your file name

        pillars = pd.read_csv(poff_tsv, "\t")

        pillars = pillars.replace("*", np.nan)

        pillars = pillars.iloc[:, 3:]

        pillars["count"] = pillars.count(1)

        pillars["family"] = pillars.index

        melted = pillars.melt(["family", "count"])

        melted[melted["count"] > 0].dropna()[["family", "value"]].to_csv("fam.txt",     # change filter value to select cutoff
                                                                         "\t",          # for min number of family members
                                                                         index=False,
                                                                         header=False)


        # # Modify fam.txt to put duplications in the same family

        fam = pd.read_csv("fam.txt", "\t", header=None)

        # separate the second column by the comma.

        fam_separated = fam[1].str.split(",", expand=True)


        # merge the separated dataframe with the regular one so that I get
        # the family names

        fam[1] = np.nan

        merged = fam.merge(fam_separated, left_index=True, right_index=True).melt(
                "0_x")

        # replace none in dataframe

        merged.replace(to_replace=[None], value=np.nan, inplace=True)

        merged["value"].dropna().str.contains(",").value_counts()


        merged[["0_x", "value"]].dropna().to_csv("fam_dup.txt", "\t", index=False,
                                                 header=False)

        # # Make files for the nucleotide sequences

        node_file = "fam_dup.txt"
        na_file = "NT.fna"

        famDict = {}

        if not os.path.exists('families/fnas'):
            os.makedirs("families/fnas")
        if not os.path.exists('families/faas'):
            os.makedirs("families/faas")

        with open(node_file) as f:
            for line in f:
                row = line.split()
                if row[0] not in famDict:
                    famDict[row[0]] = [row[1]]
                else:
                    famDict[row[0]].append(row[1])

        naseqDict = SeqIO.to_dict(SeqIO.parse(na_file, "fasta"))

        for i in famDict.keys():
            file = "families/fnas/" + i + ".fna"
            with open(file, "w") as out:
                for j in famDict[i]:
                    out.write('>' + j + '\n')
                    out.write(str(naseqDict[j].seq) + '\n')

        # # Make files for the amino acid sequences

        aa_file = "AA.faa"

        aaseqDict = SeqIO.to_dict(SeqIO.parse(aa_file, "fasta"))

        for i in famDict.keys():
            file = "families/faas/" + i + ".faa"
            with open(file, "w") as out:
                for j in famDict[i]:
                    out.write('>' + j + '\n')
                    out.write(str(aaseqDict[j].seq) + '\n')

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
        aa_aln = AlignIO.read(input.pro_align, "fasta", alphabet=generic_protein)
        na_seq = SeqIO.to_dict(SeqIO.parse(input.nucl_seqs, "fasta", alphabet=generic_dna))

        codonalign.default_codon_table.forward_table.update(
            {'NTT': "X", 'NNN': "X", 'GNN': "X", 'RTT': "X", 'GGN': "G", 'GGY': "G", 'GGK': "G", 'GGR': "G", 'CCS': "P",
             'CCN': "P", 'CCY': "P", 'CCK': "P", 'CCM': "P", 'CTN': "L", 'CTR': "L", 'CTK': "L", 'CWC': "X", 'CWS': "X",
             'WAC': "X", 'CAN': "X", 'TRS': "X", 'TRK': "X", 'RCA': "X", 'CAM': "X", 'MRC': "X", 'NCC': "X", 'NNC': "X",
             'NNA': "X", 'ANN': "X", 'CNN': "X", 'CKC': "X", 'NNG': "X", 'AAN': "X", 'NAA': "X", 'NAT': "X", 'TTN': "X",
             'NTC': "X", 'NTA': "X", 'NTG': "X", 'NGC': "X", 'NGA': "X", 'NGT': "X", 'ATN': "X", 'ATY': "I", 'AAS': "X",
             'ASC': "X", 'ATS': "X", 'TAN': "X", 'TAY': "Y", 'GAN': "X", 'AYC': "X", 'CNC': "X", 'NCG': "X", 'NGG': "X",
             'NNT': "X", 'CYC': "X", 'NAG': "X", 'ACN': "T", 'ACK': "T", 'ACR': "T", 'ACY': "T", 'ACS': "T", 'KCG': "X",
             'RCT': "X", 'RCC': "X", 'RCG': "X", 'RGC': "X", 'NAC': "X", 'TCN': "S", 'TCR': "S", 'TCK': "S", 'TCM': "S",
             'TCS': "S", 'TCY': "S", 'AGN': "X", 'SGG': "X", 'KAT': "X", 'KAC': "X", 'CGN': "R", 'CGR': "R", 'TGN': "X",
             'TNN': "X", 'GCN': "A", 'GCY': "A", 'NCT': "X", 'NCN': "X", 'GTN': "V", 'GTY': "V", 'GTK': "V", 'GTM': "V",
             'GRT': "X", 'TMC': "X", 'TKG': "X", 'TKM': "X", 'NCA': "X", 'GNG': "X", 'WGC': "X", 'GNC': "X", 'ANA': "X",
             'GNT': "X", 'ANT': "X", 'RAC': "B", 'RAA': "X", 'NGN': "X", 'KGG': "X", 'RWY': "X", 'GNA': "X", 'CNG': "X",
             'CNA': "X", 'ANG': "X", 'CNT': "X", 'TNC': "X", 'TTK': "X", 'ANC': "X", 'TNA': "X", 'TNT': "X", 'MCC': "X",
             'TNG': "X", 'NTN': "X", 'YGT': "X", 'NAN': "X", 'GAM': "X", 'GAY': "D", 'GAR': "E", 'AMC': "X", 'YTC': "X",
             'YTG': "L", 'YTA': "L", 'YAC': "X", 'WCG': "X", 'GMA': "X", 'GMC': "X", 'ACM': "T", 'ACW': "T", 'CGY': "R",
             'GYC': "X", 'GYT': "X", 'WCC': "X", 'CSC': "X", 'SCT': "X", 'SCC': "X", 'CRT': "X", 'CRC': "X", 'YCC': "X",
             'GAK': "X", 'AKC': "X", 'AWC': "X", 'AYA': "X", 'AYG': "X", 'CAY': "H", 'CAR': "Q", 'AAR': "K", 'GCR': "A",
             'CGS': "R", 'AAW': "X", 'RAT': "B", 'GGM': "G", 'GGS': "G", 'CTM': "L", 'CTW': "L", 'MGG': "R", 'GCK': "A",
             'GTR': "V", 'GTS': "V", 'ATW': "I", 'MGA': "R", 'ARC': "X", 'CGW': "R", 'GCM': "A", 'CCR': "P", 'TTR': "L",
             'TTY': "F", 'AGK': "X", 'GKG': "X", 'GYA': "X", 'ATM': "I", 'AGR': "R", 'AAY': "N", 'AAK': "X", 'GCS': "A",
             'KCT': "X", 'RTG': "X", 'GKT': "X", 'ARA': "X", 'CAK': "X", 'MCA': "X", 'MGT': "X", 'CMC': "X", 'TGY': "C",
             'CGM': "R", 'TYC': "X", 'CTY': "L", 'GGW': "G", 'KCC': "X", 'TYT': "X", 'GTW': "V", 'WCT': "X", 'CAS': "X",
             'GYG': "X", 'GCW': "A", 'MCT': "X", 'GMG': "X", 'ASA': "X", 'GCB': "A", 'WCA': "X", 'AYT': "X", 'CGH': "R",
             'RTC': "X", 'MAG': "X", 'YCG': "X", 'RGG': "X", 'MAT': "X", 'GSA': "X", 'AKG': "X", 'CSA': "X", 'MAA': "X",
             'CYG': "X", 'CKA': "X", 'CRG': "X", 'YTT': "X", 'ARG': "X", 'RGA': "X", 'KTA': "X", 'CGK': "R", 'TYA': "X",
             'ACH': "T", 'RAG': "X", 'CYA': "X", 'GMT': "X", 'GYK': "X", 'YGC': "X", 'AST': "X", 'KTT': "X", 'YAT': "X",
             'STA': "X", 'GWT': "X", 'CKG': "X", 'CSG': "X", 'TRT': "X", 'AMT': "X", 'RRM': "X", 'RRA': "X", 'MTS': "X",
             'CMT': "X", 'YYT': "X", 'YCA': "X", 'ACB': "T", 'MGR': "R", 'AKT': "X", 'AWT': "X", 'RWA': "X", 'CKT': "X",
             'TTM': "X", 'YGG': "X", 'KTG': "X", 'AGM': "X", 'CRA': "X", 'TTW': "X", 'CMA': "X", 'WAT': "X", 'TWC': "X",
             'AMA': "X", 'RTA': "X", 'TRC': "X", 'GWA': "X", 'SGA': "X", 'TYG': "X", 'CWA': "X", 'MCG': "X", 'TKT': "X",
             'GSC': "X", 'AWA': "X", 'TGK': "X", 'YCT': "X", 'SGT': "X", 'TCW': "S", 'GRA': "X", 'CMG': "X", 'ARM': "X",
             'ACD': "T", 'RRG': "X", 'GAW': "X", 'KGT': "X", 'STT': "X", 'MAC': "X", 'SAG': "X", 'ATK': "X", 'KCA': "X",
             'CCH': "P", 'SAA': "X", 'AMG': "X", 'ATR': "X", 'RGT': "X", 'GST': "X", 'AGS': "X", 'GKK': "X", 'MTC': "X",
             'RMR': "X", 'MTT': "X", 'TWT': "X", 'GRY': "X", 'TKK': "X", 'RCY': "X", 'GWC': "X", 'WCR': "X", 'GTH': "V",
             'SCA': "X", 'ATB': "X", 'SWT': "X", 'YGA': "X", 'TTS': "X", 'ARS': "X", 'AKA': "X", 'MTA': "X", 'MTG': "X",
             'SYC': "X", 'RRC': "X", 'KAG': "X", 'MKG': "X", 'AGW': "X", 'TSC': "X", 'ART': "X", 'STG': "X", 'CYT': "X",
             'WTG': "X", 'WYC': "X", 'RTS': "X", 'TMT': "X", 'KTY': "X", 'KTC': "X", 'WTT': "X", 'AGY': "S", 'GSG': "X",
             'GAS': "X", 'AAM': "X", 'YKG': "X", 'CTS': "L", 'GWW': "X", 'GMR': "X", 'GRC': "X", 'CCW': "P", 'SAC': "X",
             'GKA': "X", 'STC': "X", 'WTC': "X", 'RGM': "X", 'KGC': "X", 'WTA': "X", 'MMS': "X", 'AYR': "X", 'GRS': "X",
             'RCR': "X", 'GKC': "X", 'MGY': "X", 'KGA': "X", 'YYC': "X", 'GCD': "A", 'GYY': "X", 'GTB': "V", 'TMA': "X",
             'MAR': "X", 'KYC': "X", 'AWG': "X", 'CST': "X", 'DCC': "X", 'TGS': "X", 'TYK': "X", 'GRR': "X", 'WCS': "X",
             'SRM': "X", 'ARY': "X", 'CWG': "X", 'SAT': "X", 'AYY': "X", 'GYW': "X", 'RWG': "X", 'YTY': "X", 'RSM': "X",
             'TST': "X", 'CYK': "X", 'SMS': "X", 'MRA': "X", 'CRW': "X", 'GBC': "X", 'RAR': "X", 'ARR': "X", 'GRW': "X",
             'RAY': "B", 'SAW': "X", 'AMR': "X", 'GWS': "X", 'YKT': "X", 'GWK': "X", 'KMT': "X", 'YCK': "X", 'MGS': "X",
             'SSS': "X", 'AWK': "X", 'ASM': "X", 'SSM': "X", 'RSG': "X", 'AMM': "X", 'WMG': "X", 'RSA': "X", 'AAH': "X",
             'CSW': "X", 'MSM': "X", 'RAK': "X", 'WRC': "X", 'RGW': "X", 'RKC': "X", 'RWC': "X", 'GWR': "X", 'RRY': "X",
             'GWY': "X", 'RWS': "X", 'RMC': "X", 'SRW': "X", 'TSG': "X", 'WCY': "X", 'SCG': "X", 'RRS': "X", 'GSM': "X",
             'RCS': "X", 'RGY': "X", 'RSC': "X", 'GRG': "X", 'CKS': "X", 'SWY': "X", 'RAS': "X", 'RYC': "X", 'MRG': "X",
             'TGW': "X", 'SGC': "X", 'RMG': "X", 'GAH': "X", 'CCB': "P", 'GWG': "X", 'YTK': "X", 'KCS': "X", 'RKA': "X",
             'CAW': "X", 'TKC': "X", 'KYT': "X", 'WMT': "X", 'MGC': "X", 'ASG': "X", 'WGT': "X", 'RMA': "X", 'MKC': "X",
             'SMG': "X", 'AKY': "X", 'WMA': "X", 'SMR': "X", 'SCW': "X", 'AAB': "X", 'KAA': "X", 'YKC': "X", 'ASY': "X",
             'TMG': "X", 'RKT': "X", 'ARW': "X", 'AYM': "X", 'TYS': "X", 'SYT': "X", 'RAM': "X", 'MYA': "X", 'WMC': "X",
             'WMM': "X", 'YRT': "X", 'CKK': "X", 'SRR': "X", 'SKM': "X", 'WKG': "X", 'SGS': "X", 'WYA': "X", 'AMK': "X",
             'CWT': "X", 'TSY': "X", 'WAG': "X", 'SKT': "X", 'AMW': "X", 'RKS': "X", 'GMW': "X", 'MAK': "X", 'YSS': "X",
             'TRY': "X", 'TYY': "X", 'SRG': "X", 'SCS': "X", 'SRX': "X", 'SMC': "X", 'YCY': "X", 'KWA': "X", 'AYW': "X",
             'KCK': "X", 'YMM': "X", 'MRM': "X", 'ARK': "X", 'AMS': "X", 'TAM': "X", 'KYG': "X", 'AKM': "X", 'SRC': "X",
             'MRS': "X", 'WYW': "X", 'DGG': "X", 'KKG': "X", 'GSR': "X", 'WWC': "X", 'GKM': "X", 'TTH': "X", 'GTD': "V",
             'CTH': "L", 'CYS': "X", 'MYG': "X", 'CYM': "X", 'CRR': "X", 'GGH': "X", 'GYR': "X", 'GYM': "X", 'RGR': "X",
             'MMA': "X", 'RYG': "X", 'GRK': "X", 'GKW': "X", 'RWT': "sX", 'SKC': "X", 'AWY': "X", 'GRM': "X", 'MWG': "X"})

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
        "hyphy fubar --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log} || touch {output.log} {output.json} {output.cache}"

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
        "hyphy absrel --alignment {input.align} --tree {input.tree} --output {output.json} > {output.log} || touch {output.log} {output.json}"


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
            for currentFile in input[0]:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0**" not in line:
                                result = re.search('found(.*)branches', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")


"""
checkpoint move_absrel:
    input:
        "final_results/fams_absrel.txt"
    output:
        directory("families_absrel")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=False)
        families = fams["family"]
        families_in_dir = glob.glob("families_fubar/**/*")
        if not os.path.exists('families_fubar/fnas'):
            os.makedirs('families_absrel/fnas')
        if not os.path.exists('families_fubar/faas'):
            os.makedirs('families_absrel/faas')
        if not os.path.exists('families_fubar/logs'):
            os.makedirs('families_absrel/logs')
        if not os.path.exists('families_absrel/trees'):
            os.makedirs('families_absrel/trees')
        if not os.path.exists('families_absrel/codon_alns'):
            os.makedirs('families_fubar/codon_alns')
        for i in range(0,len(families_in_dir)):
            if int(families_in_dir[i].split(".")[0].split("/")[-1]) in list(families):
                copyfile(families_in_dir[i], "families_absrel/"+families_in_dir[i].split("/",1)[1])

def aggregate_absrel(wildcards):
    checkpoint_output = checkpoints.move_fubar.get(**wildcards).output[0]
    return expand("families_absrel/logs/family_{j}.aln.codon.ABSREL.log",
           j=glob_wildcards(os.path.join(checkpoint_output, "family_{j}.tree")).j)

rule final_stats:
    input:
        aggregate_absrel
    output:
        "final_results/fams_absrel.txt"
    run:
        with open(output[0], "w") as out:
            for currentFile in input[0]:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0**" not in line:
                                result = re.search('found(.*)branches', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")
"""

