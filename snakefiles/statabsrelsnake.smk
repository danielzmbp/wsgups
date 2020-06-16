import os
from shutil import move
import pandas as pd

FAM, = glob_wildcards("families_absrel/family_{fam}.aln.codon.ABSREL.log")

rule final:
    input:
        dynamic("families_absrel/faas/{fam}.faa")

rule finalStatistics:
    input:
        json = dynamic(
            "families_fubar/codon_alns/{fam}.aln.codon.ABSREL.json"),
        log = dynamic(
            "families_fubar/logs/{fam}.aln.codon.ABSREL.log")
    output:
        "final_results/absrel_selection.txt"
    run:
        with open(output[0], "w") as out:
            out.write("family numSitesUnderSelection\n")
            for currentFile in input.log:
                with open(currentFile) as f:
                    for line in f:
                        if "Likelihood" in line:
                            if "**0**" not in line:
                                result = re.search('found(.*)branches', line)
                                out.write(currentFile.split(
                                    "/")[-1].split(".")[0] + " " + result.group(1) + "\n")
                            else:
                                for files in glob.glob(currentFile.split(".")[0] + ".*"):
                                    os.remove(file)
rule move_files:
    input:
        "final_results/absrel_selection.txt"
    output:
        dynamic("families_absrel/faas/{fam}.faa")
    run:
        fams = pd.read_csv(input[0],"\s+",index_col=False)
        families = fams["family"]
        families_in_dir = os.listdir("families/**/*")    # not sure if this would work
        os.mkdir("families_absrel/")
        folders = ["faas/","fnas/","codon_alns/", "alns/", "jsons/",
                   "logs/", "trees/", "cleaned_faas/", "cleaned_fnas/"]

        for i in range(0,len(families_in_dir)):
            if families_in_dir[i].split(".")[0] in list(map(str, families)):
                for folder in folders:
                    move("families_fubar/" + folder + families_in_dir[i],
                         "families_absrel/" + folder + families_in_dir[i])
