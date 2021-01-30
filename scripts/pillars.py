# coding: utf-8

# ## Homology groups OGOB study
# Use this groups to create protein families than can be alined and tested for
# positive selection later.

import os
import pandas as pd
import numpy as np
from Bio import SeqIO


poff_tsv = snakemake.input[0]  # replace with your file name

pillars = pd.read_csv(poff_tsv, "\t").replace("*", np.nan).iloc[:, 3:]

pillars["count"] = pillars.count(1)

pillars["family"] = pillars.index

melted = pillars.melt(["family", "count"])

# change filter value to select cutoff
# for min number of family members

melted[melted["count"] > 4].dropna()[["family", "value"]].to_csv("fam.txt",
                                                                 "\t",
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
