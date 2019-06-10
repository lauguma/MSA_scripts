#!/usr/bin/python



from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import re
import os



# USAGE: python multifasta2ffn.py -fasta < fasta extension of the CDS files >  -out < output fasta extension of the single CDS fasta files >


def main():

    parser = ArgumentParser()
    parser.add_argument("-fasta", dest="fasta_suffix", help="multifasta files extension containing CDS", type=str)
    parser.add_argument("-out", dest="fasta_out", help="gene fasta files output extension", type=str)

    args = parser.parse_args()

    path = os.getcwd() # working directory
    suffix = str(args.fasta_suffix)  # suffix of multifasta files containing all genes
    multi_files = [f for f in os.listdir(path) if f.endswith(suffix)]
    n_list = range(0, len(multi_files))

    # Dictionary with key as strain name and value as a list of all genes in each strain
    d = dict()
    for i in n_list:
        d["lista_%s" %i] = []

    # Set one list with the first list just to go throw them later and compare
    common = list()
    f = os.path.join(path,multi_files[0])
    in_file = SeqIO.parse(f,"fasta")
    for line in in_file:
        gen = line.id.split("_")[0]
        common.append(gen)

    for x in range(len(multi_files)):
        #print multi_files[x]
        f = os.path.join(path,multi_files[x])
        in_file = SeqIO.parse(f,"fasta")
        strain = multi_files[x].split("_")[0]
        for line in in_file:
            gen = line.id.split("_")[0]
            d["lista_%s" %x].append(gen)
        common = list(set(common).intersection(d["lista_%s" %x]))

        in_file.close()

    print "Nr. of common genes among multifasta files provided: ", len(common)
    print len(common), " fasta files will be created in the working directory"

    # Read again multifasta files, if gen in common gene list write gene in a new file containing all strain sequences for a same gene
    for x in range(len(multi_files)):
        f = os.path.join(path,multi_files[x])
        in_file = SeqIO.parse(f,"fasta")
        strain = multi_files[x].split("_")[0]    
        for line in in_file:
            gen = line.id.split("_")[0]  ##
            if gen in common:
                sequence_object = Seq(str(line.seq), IUPAC.unambiguous_dna)
                record = SeqRecord(sequence_object,id=strain, description="")
                out_name = gen + str(args.fasta_out)
                with open(out_name, "a") as out_file:
                    SeqIO.write(record,out_file,"fasta")
                
            
        in_file.close()

        

    


if __name__ == "__main__" :
    main()

