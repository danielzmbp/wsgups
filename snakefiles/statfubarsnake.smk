from shutil import copyfile
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

rule move_files:
    input:
        "final_results/fams_selection.txt"
    output:
        dynamic("families_fubar/faas/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=None)
        families = fams.iloc[:,0]
        families_in_dir = glob.glob("families/**/*")

        os.mkdir("families_fubar")
        os.mkdir("families_fubar/codon_alns")
        os.mkdir("families_fubar/faas")
        os.mkdir("families_fubar/fnas")
        os.mkdir("families_fubar/logs")
        
        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0].split("/")[-1] in list(map(str,families)):
                copyfile(families_in_dir[i], "families_fubar/" + families_in_dir[i].split("/",1)[-1])
