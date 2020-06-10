from shutil import move
import pandas as pd
import re
import glob
import os

FAM, = glob_wildcards("families/logs/{fam}.aln.codon.FUBAR.log")

rule final:
    input:
        dynamic("families_fubar/faas/{fam}.faa")

rule finalStatistics:
    input:
        log = expand("families/logs/{fam}.aln.codon.FUBAR.log", fam=FAM)
    output:
        "final_results/fams_selection.txt",
    run:
        with open(output[0], "w") as out:
            for currentFile in input.log:
                with open(currentFile) as f:
                    for line in f:
                        if "## FUBAR" in line:
                            if "no" not in line:
                                result = re.search('inferred(.*)sites', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n") # family num_selected_sites
                            else:
                                for file in glob.glob(currentFile.split(".")[0] + ".*"):
                                    os.remove(file)

rule move_files:
    input:
        "final_results/fams_selection.txt"
    output:
        dynamic("families_fubar/faas/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=None)
        families = fams.iloc[:,0]
        families_in_dir = os.listdir("families")

        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0] in list(families):
                for folder in glob.glob("families/*"):
                    move("families/" + folder + families_in_dir[i],
                         "families_fubar/" + folder + families_in_dir[i])
