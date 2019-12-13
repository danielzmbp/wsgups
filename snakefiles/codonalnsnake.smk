from Bio import SeqIO
from Bio import codonalign
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein


FAM, = glob_wildcards("families/family_{fam}.faa")


rule final:
    input:
        expand("families/family_{fam}.aln.codon", fam= FAM)


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
        nucl_seqs = "families/family_{fam}.fna"
    output:
        alignment = "families/family_{fam}.aln.codon"
    run:
        aa_aln = AlignIO.read(input.pro_align, "fasta",alphabet= generic_protein)
        na_seq = SeqIO.to_dict(SeqIO.parse(input.nucl_seqs, "fasta", alphabet= generic_dna))


        codonalign.default_codon_table.forward_table["NTT"] = "X"
        codonalign.default_codon_table.forward_table["NNN"] = "X" #quick fix for Ns in the dataset
        codonalign.default_codon_table.forward_table["GNN"] = "X"
        codonalign.default_codon_table.forward_table["GGN"] = "G"
        codonalign.default_codon_table.forward_table["GGY"] = "G"
        codonalign.default_codon_table.forward_table["CCN"] = "P"
        codonalign.default_codon_table.forward_table["CTN"] = "L"
        codonalign.default_codon_table.forward_table["CTR"] = "L"
        codonalign.default_codon_table.forward_table["CAN"] = "X"
        codonalign.default_codon_table.forward_table["NCC"] = "X"
        codonalign.default_codon_table.forward_table["NNC"] = "X"
        codonalign.default_codon_table.forward_table["NNA"] = "X"
        codonalign.default_codon_table.forward_table["ANN"] = "X"
        codonalign.default_codon_table.forward_table["CNN"] = "X"
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
        codonalign.default_codon_table.forward_table["TAN"] = "X"
        codonalign.default_codon_table.forward_table["GAN"] = "X"
        codonalign.default_codon_table.forward_table["CNC"] = "X"
        codonalign.default_codon_table.forward_table["NCG"] = "X"
        codonalign.default_codon_table.forward_table["NGG"] = "X"
        codonalign.default_codon_table.forward_table["NNT"] = "X"
        codonalign.default_codon_table.forward_table["NAG"] = "X"
        codonalign.default_codon_table.forward_table["ACN"] = "T"
        codonalign.default_codon_table.forward_table["NAC"] = "X"
        codonalign.default_codon_table.forward_table["TCN"] = "S"
        codonalign.default_codon_table.forward_table["AGN"] = "X"
        codonalign.default_codon_table.forward_table["CGN"] = "R"
        codonalign.default_codon_table.forward_table["TGN"] = "X"
        codonalign.default_codon_table.forward_table["TNN"] = "X"
        codonalign.default_codon_table.forward_table["GCN"] = "A"
        codonalign.default_codon_table.forward_table["NCT"] = "X"
        codonalign.default_codon_table.forward_table["GTN"] = "V"
        codonalign.default_codon_table.forward_table["GTY"] = "V"
        codonalign.default_codon_table.forward_table["NCA"] = "X"
        codonalign.default_codon_table.forward_table["GNG"] = "X"
        codonalign.default_codon_table.forward_table["GNC"] = "X"
        codonalign.default_codon_table.forward_table["ANA"] = "X"
        codonalign.default_codon_table.forward_table["GNT"] = "X"
        codonalign.default_codon_table.forward_table["CNA"] = "X"
        codonalign.default_codon_table.forward_table["CNG"] = "X"
        codonalign.default_codon_table.forward_table["GNA"] = "X"
        codonalign.default_codon_table.forward_table["TNC"] = "X"
        codonalign.default_codon_table.forward_table["ANG"] = "X"
        codonalign.default_codon_table.forward_table["NGN"] = "X"
        codonalign.default_codon_table.forward_table["ANC"] = "X"
        codonalign.default_codon_table.forward_table["NAN"] = "X"
        codonalign.default_codon_table.forward_table["TNA"] = "X"
        codonalign.default_codon_table.forward_table["TNT"] = "X"
        codonalign.default_codon_table.forward_table["ANT"] = "X"
        codonalign.default_codon_table.forward_table["TNG"] = "X"
        codonalign.default_codon_table.forward_table["NTN"] = "X"
        codonalign.default_codon_table.forward_table["CNT"] = "X"
        codonalign.default_codon_table.forward_table["NCN"] = "X"
        codonalign.default_codon_table.forward_table["YGT"] = "X"
        codonalign.default_codon_table.forward_table["RWY"] = "X"
        codonalign.default_codon_table.forward_table["ACK"] = "T"
        codonalign.default_codon_table.forward_table["TAY"] = "Y"
        codonalign.default_codon_table.forward_table["GAM"] = "X"
        codonalign.default_codon_table.forward_table["ATS"] = "X"
        codonalign.default_codon_table.forward_table["YTC"] = "X"
        codonalign.default_codon_table.forward_table["WCG"] = "X"
        codonalign.default_codon_table.forward_table["CKC"] = "X"


        align = codonalign.build(aa_aln, na_seq, max_score= 20)


        for record in range(0,len(align)):
            align._records[record].description = "" # removes description by looping through the records


        SeqIO.write(align, output.alignment, "fasta")

