#!/usr/bin/python



# USAGE: python backTranslate.py

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import os

def main():

    parser = ArgumentParser()
    parser.add_argument("-nt", dest="suffix_dna", help="extension of nucleotide sequences", type=str)
    parser.add_argument("-ali",dest="suffix_aa", help="extension of aminoacid multiple alignment files", type=str)
    parser.add_argument("-out",dest="suffix_out", help="extension of nucleotide multiple alignment output files", type=str)

    args = parser.parse_args()

    path = os.getcwd() # working directory
    multi_files = [f for f in os.listdir(path) if f.endswith(args.suffix_dna)]
    ali_files = [f for f in os.listdir(path) if f.endswith(args.suffix_aa)]

    for x in ali_files:
        gene_name = x.replace(args.suffix_aa,"")
        ali = SeqIO.parse(x,"fasta")
        out_file = open(gene_name+"."+args.suffix_out,"w")
        
        for s in ali:
            aa_seq = str(s.seq)
            dna = SeqIO.parse(gene_name+args.suffix_dna,"fasta")
            for cds in dna:
                if cds.id == s.id:
                    #print cds.id
                    nt_seq = str(cds.seq)
                    c1 = 0 # codon start
                    c2 = 3 # codon end
                    chain = ""
                    for aa in aa_seq:
                        if aa != "-":
                            codon = nt_seq[c1:c2]
                            c1 += 3
                            c2 += 3
                            chain += codon
                        else:
                            chain += aa*3
                            


                    sequence_object = Seq(chain, IUPAC.unambiguous_dna)
                    record = SeqRecord(sequence_object,id=cds.id, description="")
                    SeqIO.write(record,out_file,"fasta")
                            
            dna.close()  
        ali.close()
        out_file.close()
        


if __name__ == "__main__" :
    main()
