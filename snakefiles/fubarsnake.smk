from Bio import SeqIO
from Bio import codonalign
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
import os

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

rule clean_faa:
    input:
        "families/family_{fam}.faa"
    output:
        "families/family_{fam}.faa.cleaned"
    run:
        new_sequences = []
        for record in SeqIO.parse(input[0], "fasta"):
            if record.seq.find("Z")!=-1:
                position = record.seq.find("Z")
                record.seq = record.seq.tomutable()
                record.seq[position] = "X"
                new_sequences.append(record)
            elif record.seq.find("J")!=-1:
                position = record.seq.find("J")
                record.seq = record.seq.tomutable()
                record.seq[position] = "X"
                new_sequences.append(record)
            else:
                new_sequences.append(record)
        SeqIO.write(new_sequences, output[0], "fasta")

rule mafft:
    input:
        "families/family_{fam}.faa.cleaned"
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

        codonalign.default_codon_table.forward_table.update({"NTT":"X","NNN":"X",
        "GNN":"X","RTT":"X","GGN":"G","GGY":"G","GGK":"G","GGR":"G","CCS":"P","CCN":"P",
        "CCY":"P","CCK":"P","CCM":"P","CTN":"L","CTR":"L","CTK":"L","CWC":"X","CWS":"X",
        "WAC":"X","CAN":"X","TRS":"X","TRK":"X","RCA":"X","CAM":"X","MGR":"X","MRC":"X",
        "NCC":"X","NNC":"X","NNA":"X","ANN":"X","CNN":"X","CKC":"X","NNG":"X","AAN":"X",
        "NAA":"X","NAT":"X","TTN":"X","NTC":"X","NTA":"X","NTG":"X","NGC":"X","NGA":"X",
        "NGT":"X","ATN":"X","ATY":"I","AAS":"X","ASC":"X","ATS":"X","TAN":"X","TAY":"Y",
        "GAN":"X","AYC":"X","CNC":"X","NCG":"X","NGG":"X","NNT":"X","CYC":"X","NAG":"X",
        "ACN":"T","ACK":"T","ACR":"T","ACY":"T","ACS":"T","KCG":"X","RCT":"X","RCC":"X",
        "RCG":"X","RGC":"X","NAC":"X","TCN":"S","TCR":"S","TCK":"S","TCM":"S","TCS":"S",
        "TCY":"S","AGN":"X","SGG":"X","KAT":"X","KAC":"X","CGN":"R","CGR":"R","TGN":"X",
        "TNN":"X","GCN":"A","GCY":"A","NCT":"X","NCN":"X","GTN":"V","GTY":"V","GTK":"V",
        "GTM":"V","GRT":"X","TMC":"X","TKG":"X","TKM":"X","NCA":"X","GNG":"X","WGC":"X",
        "GNC":"X","ANA":"X","GNT":"X","ANT":"X","RAC":"B","RAA":"X","NGN":"X","KGG":"X",
        "RWY":"X","GNA":"X","CNG":"X","CNA":"X","ANG":"X","CNT":"X","TNC":"X","TTK":"X",
        "ANC":"X","TNA":"X","TNT":"X","MCC":"X","TNG":"X","NTN":"X","YGT":"X","NAN":"X",
        "GAM":"X","GAY":"D","GAR":"E","AMC":"X","YTC":"X","YTG":"L","YTA":"L","YAC":"X",
        "WCG":"X","GMA":"X","GMC":"X","ACM":"T","ACW":"T","CGY":"R","GYC":"X","GYT":"X",
        "WCC":"X","CSC":"X","SCT":"X","SCC":"X","CRT":"X","CRC":"X","YCC":"X","GAK":"X",
        "AKC":"X","AWC":"X","AYA":"X","AYG":"X","CAY":"H","CAR":"Q","AAR":"K","GCR":"A",
        "CGS":"R","AAW":"X","RAT":"B","GGM":"G","GGS":"G","CTM":"L","CTW":"L","MGG":"R",
        "GCK":"A","GTR":"V","GTS":"V","ATW":"I","CTY":"C","MGA":"R","CAR":"Q","ARC":"X",
        "CGW":"R","GCM":"A","CCS":"P","CCR":"P","TTR":"L","TTY":"F","AGK":"X","GKG":"X",
        "GYA":"X","ATM":"I","AGR":"R","AAY":"N","AAK":"X","GCS":"A","KCT":"X","TTY":"F",
        "RTG":"X","GKT":"X","ARA":"X","CAK":"X","MCA":"X","MGT":"X","GTS":"V","CMC":"X",
        "TGY":"C","CGM":"R","GTR":"V","TYC":"X","CTY":"L","GGW":"G","KCC":"X","TYT":"X",
        "GTW":"V","WCT":"X","GGS":"G","CAS":"X","GYG":"X","GCW":"A","MCT":"X","GMG":"X",
        "ASA":"X","GCB":"A","WCA":"X","AYT":"X","CGH":"R","RTC":"X","MAG":"X","YCG":"X",
        "RGG":"X","MAT":"X","GSA":"X","AKG":"X","CSA":"X","MAA":"X","CYG":"X","CKA":"X",
        "CRG":"X","YTT":"X","ARG":"X","RGA":"X","KTA":"X","CGK":"R","MGA":"R","TYA":"X",
        "ACH":"T","RAG":"X","CYA":"X","GMT":"X","GYK":"X","YGC":"X","AST":"X","KTT":"X",
        "YAT":"X","STA":"X","GWT":"X","CKG":"X","CSG":"X","TRT":"X","AMT":"X","RRM":"X",
        "RRA":"X","MTS":"X","CMT":"X","YYT":"X","YCA":"X","ACB":"T","MGR":"R","AKT":"X",
        "AWT":"X","RWA":"X","CKT":"X","TTM":"X","YGG":"X","KTG":"X","AGM":"X","CRA":"X",
        "TTW":"X","STA":"X","RTG":"X","CMA":"X","WAT":"X","TWC":"X","CGK":"R","AMA":"X",
        "RTA":"X","TRC":"X","GWA":"X","SGA":"X","TYG":"X","CWA":"X","MCG":"X","TKT":"X",
        "GSC":"X","AWA":"X","TGK":"X","YCT":"X","SGT":"X","TCW":"S","GRA":"X","CMG":"X",
        "ARM":"X","ACD":"T","RRG":"X","GAW":"X","KGT":"X","STT":"X","ATN":"X","MAC":"X",
        "SAG":"X","ATK":"X","KCA":"X","CCH":"P","SAA":"X","AMG":"X","ATR":"X","RGT":"X",
        "GST":"X","AGS":"X","GKK":"X","ATK":"X","MTC":"X","RMR":"X","MTT":"X","TWT":"X",
        "GRY":"X","TKK":"X","RCY":"X","GWC":"X","WCR":"X","GTH":"V","SCA":"X","ATB":"X",
        "SWT":"X","YGA":"X","TTS":"X","ARS":"X","AKA":"X","MTA":"X","MTG":"X","SYC":"X",
        "RRC":"X","KAG":"X","MKG":"X","AGW":"X","TSC":"X","ART":"X","STG":"X","CYT":"X",
        "WTG":"X","WYC":"X","RTS":"X","TMT":"X","KTY":"X","KTC":"X","WTT":"X","AGY":"S",
        "GSG":"X","GAS":"X","AAM":"X","YKG":"X","CTS":"L","AGY":"S","GWW":"X","GMR":"X",
        "GRC":"X","CCW":"P","SAC":"X","GKA":"X","STC":"X","AYR":"X","WTC":"X","RGM":"X",
        "KGC":"X","WTA":"X","MMS":"X","AYR":"X","GRS":"X","RCR":"X","GKC":"X","MGY":"X",
        "KGA":"X","YYC":"X","GCD":"A","GYY":"X","DCC":"X","GTB":"V","TMA":"X","WCS":"X",
        "MAR":"X","GRR":"X","KYC":"X","AWG":"X","CST":"X","DCC":"X","TGS":"X","GCD":"A",
        "GYY":"X","TYK":"X","TMA":"X","GRR":"X","GTB":"V","WCS":"X","SRM":"X","ARY":"X",
        "CWG":"X","SAT":"X","AYY":"X","GYW":"X","RWG":"X","YTY":"X","RSM":"X","TST":"X",
        "CYK":"X","SMS":"X","MRA":"X","CRW":"X","GBC":"X","RAR":"X","ARR":"X","GRW":"X",
        "RAY":"B","SAW":"X","AMR":"X","GWS":"X","YKT":"X","GWK":"X","KMT":"X","YCK":"X",
        "MGS":"X","SSS":"X","AWK":"X","ASM":"X","SSM":"X","RSG":"X","AMM":"X","WMG":"X",
        "RSA":"X","AAH":"X","CSW":"X","MSM":"X","RAK":"X","WRC":"X","RGW":"X","RKC":"X",
        "RWC":"X","GWR":"X","RRY":"X","GWY":"X","RWS":"X","RMC":"X","SRW":"X","TSG":"X",
        "WCY":"X","SCG":"X","RRS":"X","GSM":"X","RCS":"X","RGY":"X","RSC":"X","GRG":"X",
        "CKS":"X","SWY":"X","RAS":"X","RYC":"X","MRG":"X","TGW":"X","SGC":"X","RMG":"X",
    	"GAH":"X","CCB":"P","GWG":"X","YTK":"X","KCS":"X","RKA":"X","CAW":"X","TKC":"X",
	    "KYT":"X","WMT":"X","MGC":"X","ASG":"X","WGT":"X","RMA":"X","MKC":"X","SMG":"X",
    	"AKY":"X","WMA":"X","SMR":"X","SCW":"X","AAB":"X","KAA":"X","YKC":"X","ASY":"X",
	    "TMG":"X","RKT":"X","ARW":"X","AYM":"X","TYS":"X","SYT":"X","RAM":"X","MYA":"X",
    	"WMC":"X","WMM":"X","YRT":"X","CKK":"X","SRR":"X","SKM":"X","WKG":"X","SGS":"X",
    	"WYA":"X","AMK":"X","CWT":"X","TSY":"X","WAG":"X","SKT":"X","AMW":"X","RKS":"X",
    	"GMW":"X","MAK":"X","YSS":"X","TRY":"X","TYY":"X","SRG":"X","SCS":"X","SRX":"X",
    	"SMC":"X","YCY":"X","KWA":"X","AYW":"X","KCK":"X","YMM":"X","MRM":"X","ARK":"X",
    	"AMS":"X","TAM":"X","KYG":"X","AKM":"X","SRC":"X","MRS":"X","WYW":"X","DGG":"X",
	    "KKG":"X","GSR":"X","MRS":"X","WWC":"X","GKM":"X","TTH":"X","GTD":"V","CTH":"L",
    	"CYS":"X","MYG":"X","CYM":"X","CRR":"X","GGH":"X"})



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
        """
        tmpThing=$(find /.. -name FUBAR.bf|tail -n 1);(echo 1;echo {input.align}; echo {input.tree};
        echo 20;echo 5;echo 3;echo 0.5 )|hyphy $tmpThing > {output.log} || touch {output.log} {output.json}
        """

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
