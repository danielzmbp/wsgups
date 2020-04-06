import os
from shutil import move
import pandas as pd
import re
import glob

FAM, = glob_wildcards("families/family_{fam}.aln.codon.FUBAR.log")

rule final:
    input:
        dynamic("families_fubar/{fam}.faa")

rule finalStatistics:
    input:
        log = expand("families/family_{fam}.aln.codon.FUBAR.log", fam=FAM)
    output:
        "final_results/fams_selection.txt",
    run:
        with open(output[0], "w") as out:
            out.write("family num_SUS\n")
            for currentFile in input.log:
                with open(currentFile) as f:
                    for line in f:
                        if "## FUBAR" in line:
                            if "no" not in line:
                                result = re.search('inferred(.*)sites', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")
                            else:
                                for file in glob.glob(currentFile.split(".")[0] + ".*"):
                                    os.remove(file)

rule move_files:
    input:
        "final_results/fams_selection.txt"
    output:
        dynamic("families_fubar/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False)
        families = fams["family"]
        families_in_dir = os.listdir("families")

        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0] in list(families):
                move("families/"+families_in_dir[i], "families_fubar/"+families_in_dir[i])
