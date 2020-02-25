from Bio import SeqIO
from Bio import codonalign
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
import os

wd = os.getcwd()

FAM, = glob_wildcards("families/family_{fam}.fna")

rule final:
    input:
        expand("families/family_{fam}.aln.codon.FUBAR.json", fam= FAM)

rule clean_fna:
    input:
        "families/family_{fam}.fna"
    output:
        "families/family_{fam}.fna.cleaned"
    run:
        new_sequences = []
        for record in SeqIO.parse(input[0], "fasta"):
            if len(record.seq) % 3 != 0:
                record.seq = record.seq + "N"
                new_sequences.append(record)
            else:
                new_sequences.append(record)
        SeqIO.write(new_sequences, output[0], "fasta")

rule mafft:
    input:
        "families/family_{fam}.faa"
    output:
        "families/family_{fam}.aln"
    shell:
        "mafft --auto --thread 1 {input} > {output}"
rule fasttree:
    input:
        "families/family_{fam}.aln"
    output:
        "families/family_{fam}.tree"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -nosupport {input} > {output} || true"
rule codonaln:
    input:
        pro_align = "families/family_{fam}.aln",
        nucl_seqs = "families/family_{fam}.fna.cleaned"
    output:
        alignment = "families/family_{fam}.aln.codon"
    run:
        aa_aln = AlignIO.read(input.pro_align, "fasta", alphabet=generic_protein)
        na_seq = SeqIO.to_dict(SeqIO.parse(input.nucl_seqs, "fasta", alphabet=generic_dna))

        codonalign.default_codon_table.forward_table["NTT"] = "X"
        codonalign.default_codon_table.forward_table["NNN"] = "X"  # quick fix for Ns in the dataset
        codonalign.default_codon_table.forward_table["GNN"] = "X"
        codonalign.default_codon_table.forward_table["GGN"] = "G"
        codonalign.default_codon_table.forward_table["GGY"] = "G"
        codonalign.default_codon_table.forward_table["GGK"] = "G"
        codonalign.default_codon_table.forward_table["CCN"] = "P"
        codonalign.default_codon_table.forward_table["CCY"] = "P"
        codonalign.default_codon_table.forward_table["CCK"] = "P"
        codonalign.default_codon_table.forward_table["CCM"] = "P"
        codonalign.default_codon_table.forward_table["CTN"] = "L"
        codonalign.default_codon_table.forward_table["CTR"] = "L"
        codonalign.default_codon_table.forward_table["CTK"] = "L"
        codonalign.default_codon_table.forward_table["CWC"] = "X"
        codonalign.default_codon_table.forward_table["CWS"] = "X"
        codonalign.default_codon_table.forward_table["WAC"] = "X"
        codonalign.default_codon_table.forward_table["CAN"] = "X"
        codonalign.default_codon_table.forward_table["TRS"] = "X"
        codonalign.default_codon_table.forward_table["TRK"] = "X"
        codonalign.default_codon_table.forward_table["CAM"] = "X"
        codonalign.default_codon_table.forward_table["MGR"] = "X"
        codonalign.default_codon_table.forward_table["MRC"] = "X"
        codonalign.default_codon_table.forward_table["NCC"] = "X"
        codonalign.default_codon_table.forward_table["NNC"] = "X"
        codonalign.default_codon_table.forward_table["NNA"] = "X"
        codonalign.default_codon_table.forward_table["ANN"] = "X"
        codonalign.default_codon_table.forward_table["CNN"] = "X"
        codonalign.default_codon_table.forward_table["CKC"] = "X"
        codonalign.default_codon_table.forward_table["NNG"] = "X"
        codonalign.default_codon_table.forward_table["AAN"] = "X"
        codonalign.default_codon_table.forward_table["NAA"] = "X"
        codonalign.default_codon_table.forward_table["NAT"] = "X"
        codonalign.default_codon_table.forward_table["TTN"] = "X"
        codonalign.default_codon_table.forward_table["NTC"] = "X"
        codonalign.default_codon_table.forward_table["NTA"] = "X"
        codonalign.default_codon_table.forward_table["NTG"] = "X"
        codonalign.default_codon_table.forward_table["NGC"] = "X"
        codonalign.default_codon_table.forward_table["NGA"] = "X"
        codonalign.default_codon_table.forward_table["NGT"] = "X"
        codonalign.default_codon_table.forward_table["ATN"] = "X"
        codonalign.default_codon_table.forward_table["ATY"] = "I"
        codonalign.default_codon_table.forward_table["AAS"] = "X"
        codonalign.default_codon_table.forward_table["ASC"] = "X"
        codonalign.default_codon_table.forward_table["ATS"] = "X"
        codonalign.default_codon_table.forward_table["TAN"] = "X"
        codonalign.default_codon_table.forward_table["TAY"] = "Y"
        codonalign.default_codon_table.forward_table["GAN"] = "X"
        codonalign.default_codon_table.forward_table["AYC"] = "X"
        codonalign.default_codon_table.forward_table["CNC"] = "X"
        codonalign.default_codon_table.forward_table["NCG"] = "X"
        codonalign.default_codon_table.forward_table["NGG"] = "X"
        codonalign.default_codon_table.forward_table["NNT"] = "X"
        codonalign.default_codon_table.forward_table["CYC"] = "X"
        codonalign.default_codon_table.forward_table["NAG"] = "X"
        codonalign.default_codon_table.forward_table["ACN"] = "T"
        codonalign.default_codon_table.forward_table["ACK"] = "T"
        codonalign.default_codon_table.forward_table["ACR"] = "T"
        codonalign.default_codon_table.forward_table["ACS"] = "T"
        codonalign.default_codon_table.forward_table["KCG"] = "X"
        codonalign.default_codon_table.forward_table["RCT"] = "X"
        codonalign.default_codon_table.forward_table["RCC"] = "X"
        codonalign.default_codon_table.forward_table["RCG"] = "X"
        codonalign.default_codon_table.forward_table["RGC"] = "X"
        codonalign.default_codon_table.forward_table["NAC"] = "X"
        codonalign.default_codon_table.forward_table["TCN"] = "S"
        codonalign.default_codon_table.forward_table["TCR"] = "S"
        codonalign.default_codon_table.forward_table["TCK"] = "S"
        codonalign.default_codon_table.forward_table["TCM"] = "S"
        codonalign.default_codon_table.forward_table["TCS"] = "S"
        codonalign.default_codon_table.forward_table["TCY"] = "S"
        codonalign.default_codon_table.forward_table["AGN"] = "X"
        codonalign.default_codon_table.forward_table["SGG"] = "X"
        codonalign.default_codon_table.forward_table["KAT"] = "X"
        codonalign.default_codon_table.forward_table["KAC"] = "X"
        codonalign.default_codon_table.forward_table["CGN"] = "R"
        codonalign.default_codon_table.forward_table["CGR"] = "R"
        codonalign.default_codon_table.forward_table["TGN"] = "X"
        codonalign.default_codon_table.forward_table["TNN"] = "X"
        codonalign.default_codon_table.forward_table["GCN"] = "A"
        codonalign.default_codon_table.forward_table["GCY"] = "A"
        codonalign.default_codon_table.forward_table["NCT"] = "X"
        codonalign.default_codon_table.forward_table["NCN"] = "X"
        codonalign.default_codon_table.forward_table["GTN"] = "V"
        codonalign.default_codon_table.forward_table["GTY"] = "V"
        codonalign.default_codon_table.forward_table["GTK"] = "V"
        codonalign.default_codon_table.forward_table["GTM"] = "V"
        codonalign.default_codon_table.forward_table["GRT"] = "X"
        codonalign.default_codon_table.forward_table["TMC"] = "X"
        codonalign.default_codon_table.forward_table["TKG"] = "X"
        codonalign.default_codon_table.forward_table["TKM"] = "X"
        codonalign.default_codon_table.forward_table["NCA"] = "X"
        codonalign.default_codon_table.forward_table["GNG"] = "X"
        codonalign.default_codon_table.forward_table["GNC"] = "X"
        codonalign.default_codon_table.forward_table["ANA"] = "X"
        codonalign.default_codon_table.forward_table["GNT"] = "X"
        codonalign.default_codon_table.forward_table["ANT"] = "X"
        codonalign.default_codon_table.forward_table["RAC"] = "B"
        codonalign.default_codon_table.forward_table["RAA"] = "X"
        codonalign.default_codon_table.forward_table["NGN"] = "X"
        codonalign.default_codon_table.forward_table["KGG"] = "X"
        codonalign.default_codon_table.forward_table["RWY"] = "X"
        codonalign.default_codon_table.forward_table["GNA"] = "X"
        codonalign.default_codon_table.forward_table["CNG"] = "X"
        codonalign.default_codon_table.forward_table["CNA"] = "X"
        codonalign.default_codon_table.forward_table["ANG"] = "X"
        codonalign.default_codon_table.forward_table["CNT"] = "X"
        codonalign.default_codon_table.forward_table["TNC"] = "X"
        codonalign.default_codon_table.forward_table["ANC"] = "X"
        codonalign.default_codon_table.forward_table["TNA"] = "X"
        codonalign.default_codon_table.forward_table["TNT"] = "X"
        codonalign.default_codon_table.forward_table["MCC"] = "X"
        codonalign.default_codon_table.forward_table["TNG"] = "X"
        codonalign.default_codon_table.forward_table["NTN"] = "X"
        codonalign.default_codon_table.forward_table["YGT"] = "X"
        codonalign.default_codon_table.forward_table["NAN"] = "X"
        codonalign.default_codon_table.forward_table["GAM"] = "X"
        codonalign.default_codon_table.forward_table["GAY"] = "D"
        codonalign.default_codon_table.forward_table["AMC"] = "X"
        codonalign.default_codon_table.forward_table["YTC"] = "X"
        codonalign.default_codon_table.forward_table["YTG"] = "L"
        codonalign.default_codon_table.forward_table["YTA"] = "L"
        codonalign.default_codon_table.forward_table["YAC"] = "X"
        codonalign.default_codon_table.forward_table["WCG"] = "X"
        codonalign.default_codon_table.forward_table["GMC"] = "X"
        codonalign.default_codon_table.forward_table["ACM"] = "T"
        codonalign.default_codon_table.forward_table["CGY"] = "R"
        codonalign.default_codon_table.forward_table["GYC"] = "X"
        codonalign.default_codon_table.forward_table["GYT"] = "X"
        codonalign.default_codon_table.forward_table["WCC"] = "X"
        codonalign.default_codon_table.forward_table["CSC"] = "X"
        codonalign.default_codon_table.forward_table["SCT"] = "X"
        codonalign.default_codon_table.forward_table["SCC"] = "X"
        codonalign.default_codon_table.forward_table["CRT"] = "X"
        codonalign.default_codon_table.forward_table["CRC"] = "X"
        codonalign.default_codon_table.forward_table["YCC"] = "X"
        codonalign.default_codon_table.forward_table["GAK"] = "X"
        codonalign.default_codon_table.forward_table["AKC"] = "X"
        codonalign.default_codon_table.forward_table["AWC"] = "X"
        codonalign.default_codon_table.forward_table["AYA"] = "X"
        codonalign.default_codon_table.forward_table["AYG"] = "X"
        codonalign.default_codon_table.forward_table["CAY"] = "H"


        align = codonalign.build(aa_aln, na_seq, max_score=20)

        for record in range(0, len(align)):
            align._records[record].description = ""  # removes description by looping through the records

        SeqIO.write(align, output.alignment, "fasta")

rule hyphy:
    input:
        tree = "families/family_{fam}.tree",
        align = "families/family_{fam}.aln.codon"
    output:
        json = "families/family_{fam}.aln.codon.FUBAR.json",
        log = "families/family_{fam}.aln.codon.FUBAR.log"
    conda:
        "envs/hyphy.yaml"
    shell:
        "tmpThing=$(find " + wd +
        "/.. -name FUBAR.bf|tail -n 1);(echo 1;echo " +
        wd + "/{input.align}; echo " + wd +
        "/{input.tree}; echo 20;echo 5;echo 3;echo 0.5 )|hyphy $tmpThing > {output.log} || touch {output.log} {output.json} "

rule finalStatistics:
    input:
        json = dynamic("families/family_{fam}.aln.codon.FUBAR.json"),
        log = dynamic("families/family_{fam}.aln.codon.FUBAR.log")
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
