import os
from shutil import copyfile
import pandas as pd

FAM, = glob_wildcards("families/family_{fam}.faa")

rule final:
    input:
#        "final_results/famsUnderSelection.txt"
        dynamic("families_fubar/{fam}.faa")
rule finalStatistics:
    input:
        log = expand("families/family_{fam}.aln.codon.FUBAR.log", fam=FAM)
    output:
        "final_results/famsUnderSelection.txt"
    run:
        with open(output[0], "w") as out:
            out.write("family numSitesUnderSelection\n")
            for currentFile in input.log:
                with open(currentFile) as f:
                    for line in f:
                        if "## FUBAR" in line:
                            if "no" not in line:
                                result = re.search('inferred(.*)sites', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")

rule move_files:
    input:
        "final_results/famsUnderSelection.txt"
    output:
        dynamic("families_fubar/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False)
        families = fams["family"]
        families_in_dir = os.listdir("families")


        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0] in list(families):
                copyfile("families/"+families_in_dir[i], "families_fubar/"+families_in_dir[i])
