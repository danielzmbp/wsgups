import pandas as pd
import glob
import os
import re
import phyphy
from ete3 import Tree
import numpy as np

absrel = glob.glob("families_absrel/logs/*.ABSREL.log")

family_list = []
branch_list = []
pvalue_list = []

for file in absrel:
    with open(file) as myfile:
        for line in myfile:
            if re.search(r"^\* \w.+ p-value", line):
                family, branch, pvalue = (
                    file[21:].split(".")[0],
                    line[2:].split(",")[0],
                    line.split()[-1],
                )
                family_list.append(family)
                branch_list.append(branch)
                pvalue_list.append(pvalue)

ps = pd.DataFrame(
    {"family": family_list, "branch": branch_list, "p-value": pvalue_list}
)

# ps: dataframe that includes all families, branches and pvalues of the ABSREL analysis

absrel_json = glob.glob("families_absrel/codon_alns/*.ABSREL.json")

tree_list = []
family_list = []

for file in absrel_json:
    if os.stat(file).st_size > 0:
        tree, family = (
            phyphy.Extractor(file).extract_input_tree(),
            file[27:].split(".")[0],
        )
        tree_list.append(tree)
        family_list.append(family)

tdf = pd.DataFrame({"tree": tree_list, "family": family_list})

# tdf: dataframe that includes all trees from ABSREL positives and their family

ps = ps.merge(tdf)

children_list = []

for i in range(0, len(ps)):
    t = Tree(ps.tree[i], format=1)
    node = t.search_nodes(name=ps.branch[i])[0]
    children_list.append(node.get_leaf_names())

ps["children"] = children_list

# children list includes all children from the nodes under selection

lst_col = "children"

r = pd.DataFrame(
    {
        col: np.repeat(ps[col].values, ps[lst_col].str.len())
        for col in ps.columns.drop(lst_col)
    }
).assign(**{lst_col: np.concatenate(ps[lst_col].values)})[ps.columns]

r_dd = r.drop_duplicates(subset=["family", "children"])

files = os.listdir("samples/")

taxa = []
for file in files:
    if file.endswith(".faa"):
        taxa.append(file.split(".")[0])

spec_name_list = []

for t in taxa:
    with open("samples/" + t + ".faa", "r") as file:
        for line in file:
            if line.startswith(">"):
                spec_name_list.append([line[1:].strip(), t])

spec_name = pd.DataFrame(spec_name_list, columns=["index", "species"]).set_index(
    "index"
)

spec_name["children"] = spec_name.index

spec_name.replace(r"\.", "_", regex=True, inplace=True)

spec_name.replace(r"-", "_", regex=True, inplace=True)

spec_name["ps"] = spec_name["children"].isin(r_dd["children"])

spec_name.to_csv(snakemake.output[0], index=False)

proteome = glob.glob("samples/*.tsv")

col_names = [
    "protein_accession",
    "md5",
    "length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
    "ip_accession",
    "ip_description",
    "go",
    "pathway",
]

annotations = []

for infile in proteome:
    data = pd.read_csv(infile, sep="\t", names=col_names)
    data["family"] = infile[8:].split(".")[0]
    annotations.append(data)

all_annotations = pd.concat(annotations).reset_index()

all_annotations["protein_accession"].replace(r"\.", "_", regex=True, inplace=True)

all_annotations["protein_accession"].replace(r"-", "_", regex=True, inplace=True)

all_annotations["ps"] = np.where(
    all_annotations["protein_accession"].isin(r_dd["children"]), 1, 0
)

families = pd.read_csv("fam.txt", sep=r"\s+", header=None)

families = families[1].str.split(",", expand=True).fillna(value=np.nan)

families = families.dropna(subset=[1]).melt().dropna()

dups = families["value"].replace(r"\.", "_", regex=True)

all_annotations["dups"] = np.where(
    all_annotations["protein_accession"].isin(dups), 1, 0
)

# Make GO associations file

all_annotations["go"] = all_annotations["go"].str.replace("|", ",")

all_annotations[["protein_accession", "go"]].dropna().to_csv(
    snakemake.output[1], sep="\t", index=False
)

# Write annotations to file

all_annotations.to_csv(snakemake.output[2], index=False)
