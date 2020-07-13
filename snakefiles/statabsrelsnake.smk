from shutil import copyfile
import pandas as pd
import re
import glob
import os

FAM, = glob_wildcards("families_fubar/logs/{fam}.aln.codon.ABSREL.log")

localrules: final, finalStatistics, move_files

rule final:
    input:
        dynamic("families_absrel/faas/{fam}.faa")

rule finalStatistics:
    input:
        log = expand("families_fubar/logs/{fam}.aln.codon.ABSREL.log", fam=FAM)
    output:
        "final_results/absrel_selection.txt",
    run:
        with open(output[0], "w") as out:
            for currentFile in input.log:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0**" not in line:
                                result = re.search('found(.*)branches', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")
rule move_files:
    input:
        "final_results/absrel_selection.txt"
    output:
        dynamic("families_absrel/faas/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False,header=None)
        families = fams.iloc[:,0]
        families_in_dir = glob.glob("families_fubar/**/*")

        os.mkdir("families_absrel")
        os.mkdir("families_absrel/codon_alns")
        os.mkdir("families_absrel/faas")
        os.mkdir("families_absrel/fnas")
        os.mkdir("families_absrel/logs")
	os.mkdir("families_absrel/alns")
	os.mkdir("families_absrel/trees")

        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0].split("/")[-1] in list(map(str,families)):
                copyfile(families_in_dir[i], "families_absrel/" + families_in_dir[i].split("/",1)[-1])
