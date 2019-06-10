#!/usr/bin/python
# -*- coding: utf-8 -*-



# USAGE: python cat_MSA.py -msa < .cd.mafft > -out < cat_msa.fasta > 

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import os

def main():

    parser = ArgumentParser()
    parser.add_argument("-msa", dest="suffix_msa", help="suffix of the multiple sequence alignment files", type=str)
    parser.add_argument("-out",dest="out_name", help="name of the concatenated alignment file", type=str)

    args = parser.parse_args()

    path = os.getcwd() # working directory
    multi_files = [f for f in os.listdir(path) if f.endswith(args.suffix_msa)]

    
    


    # create dictionary with all sequences to be concatenated
    # Genome name format in fasta file written before '_' 
    # Example : SCERS288C_15G2550
    

    f_in = SeqIO.parse(multi_files[1],"fasta")
    spp = dict()
    for line in f_in:
        i = line.id.split("_")[0]
        spp["%s" %i] = ""
        
    f_in.close()


    for f in multi_files:
        fasta = SeqIO.parse(f, "fasta")
        for s in fasta:
            tag = s.id.split("_")[0]
            spp[tag] += str(s.seq)
                    

        fasta.close()

    cat_file = open(args.out_name,"w")
    for k, v in spp.items():
        sequence_object = Seq(str(v), IUPAC.unambiguous_dna)
        record = SeqRecord(sequence_object,id=k, description="")
        SeqIO.write(record,cat_file,"fasta")
    cat_file.close()

                    



if __name__ == "__main__" :
    main()
