FAM, = glob_wildcards("families/family_{fam}.faa")


rule final:
    input:
        "final_results/famsUnderSelection.txt"
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
