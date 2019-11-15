FAM, = glob_wildcards("families_absrel/family_{fam}.faa")


rule final:
    input:
        "final_results/famsUnderAbsrel.txt"

rule finalStatistics:
    input:
        json = dynamic(
            "families_fubar/family_{fam}.aln.codon.ABSREL.json"),
        log = dynamic(
            "families_fubar/family_{fam}.aln.codon.ABSREL.log")
    output:
        "final_results/famsUnderAbsrel.txt"
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

