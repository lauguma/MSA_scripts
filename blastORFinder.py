#!/usr/local/bin/python
# -*- coding: utf-8 -*-


__author__ = "Laura G. Macias"
__email__ = "laugmacias@gmail.com"


import csv
from argparse import ArgumentParser
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# USAGE python get_seqs_blast.py contigs.fna file.blastn

'''
The aim of this script is to extract nucleotide coding sequences from a contigs file
according to a tabular file obtained from blast searches.
e. g. if you want to get a gene sequence from a contigs file, you use a gene sequence from
a close organism to get the location of this gene in your contigs file and this script uses
this information to automatically extract ORF corresponding to that gene sequence.
It might create two type of results:
- complete CDS if the ORF found is complete: it starts by ATG and ends by a stop codon
- partial CDS: if the ORF it is not complete.
'''


# Blast results class
class BlastResult:
    query = ""
    contig = ""
    length = int()
    identity = float()
    seq_start = int()
    seq_end = int()
    strand = ""
    def return_contig(self):
        return self.contig
    # Strand attribute: direct strand or reverse strand
    def which_strand(self):
        if self.seq_start < self.seq_end:
            self.strand = "+"
        elif self.seq_start > self.seq_end:
            self.strand = "-"
        return self.strand
    def contig_range(self):
        return (self.seq_start, self.seq_end)


# Function to get dna seq if gene is in reverse complement strand
def SeqStrand(sequence, strand):
    if strand == "-":
        sequence = sequence.reverse_complement()
    else:
        sequence = sequence

    return sequence
    

# Function to extract one sequence from multifasta file by ID
def extractSeq(fasta_file, seqID):
    fasta = SeqIO.parse(fasta_file, "fasta")
    for seq in fasta:
        if seq.id == seqID:
            return seq.seq

    fasta.close()

# Translate the three open-reading frames to amino-acid sequences
def ReadingFrames(seq):
    readingframes = [Seq.translate(seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]
    results = []
    for frame in readingframes:
        for peptide in frame.split('*'): # Split translation over stopcodons
            if len(peptide) > 30:
                results.append(peptide)

    orf = max(results, key=len) # most likely ORF, max length

    if orf.startswith("M"):   # the ORF chosen will be the largest starting by ATG - Met
        return orf 
    else:
        return False


def main():
    
    parser = ArgumentParser()
    parser.add_argument("-fasta", dest="contigs", help="multi fasta file from you want to extract coding sequences found with blast", type=str)
    parser.add_argument("-blast", dest="outfmt", help="blast results file (outfmt 6 format)", type=str)
    args = parser.parse_args()
    
    fasta = args.contigs
    met = "ATG"
    with open(args.outfmt,"r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            blast = BlastResult()
            blast.query = row[0]
            blast.contig = row[1]
            blast.identity = row[2]
            blast.length = int(row[3])
            blast.seq_start =  int(row[8])
            blast.seq_end = int(row[9])
            blast.strand = blast.which_strand()
            contig_name = blast.return_contig() # contig name where orf will be searched
            gene_range = blast.contig_range()  # base range
            gene_range = sorted(gene_range)
            
            contig_seq = extractSeq(fasta, contig_name)
            
            if blast.length > 300: # quality control (arbitrary, can be changed)
                ini = gene_range[0] - 1
                fin = gene_range[1] - 1
                
                info = [blast.query,blast.identity,str(blast.length),contig_name]
                orf_final = contig_seq[ini:fin]   # complete initial sequence, if orf are not found I will kept this sequence (partial CDS result)
                header = "_".join(info)
                
                gene = contig_seq[ini:fin]
                seguir = True
                while seguir == True:
                    gene = SeqStrand(gene, blast.strand)
                    orf = ReadingFrames(gene)
                    if orf != False:
                        orf_final = gene  # if I found a complete ORF I kept this result
                        header = "_".join(info)
                        seguir = False

                    else:
                        ini = ini - 1
                        fin = fin + 1
                        if ini < 0 or fin > len(contig_seq):
                            header = header+"_PARTIAL" # if ORF found is not complete --> partial ORF result will be returned
                            seguir = False
                        else:
                            gene = contig_seq[ini:fin]

                sequence_object = Seq(str(orf_final), IUPAC.unambiguous_dna)
                record = SeqRecord(sequence_object,id=header, description="")
                
                
                if not header.endswith("_PARTIAL"):
                    with open("completeCDS_"+contig_name+".fasta","a") as out_file:
                        SeqIO.write(record,out_file,"fasta")

                else:
                    with open("partialCDS_"+contig_name+".fasta","a") as out_file2:
                        SeqIO.write(record,out_file2,"fasta")
                    
                

        

if __name__ == "__main__" :
    main()
