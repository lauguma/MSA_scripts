#!/usr/bin/python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio import AlignIO


# USAGE python msa2phy.py -msa < MSA_file.fasta > -out < outFile.phy >


def main():

    parser = ArgumentParser()
    parser.add_argument("-msa", dest="msa_input", \
                        help="fasta file containing MSA in fasta format", type=str)
    parser.add_argument("-out",dest="phy_out",  \
                        help="File output name where msa will be written in Phylip sequential format", type=str)

    args = parser.parse_args()


    input_handle = open(args.msa_input, "rU")
    output_handle = open(args.phy_out, "w")


    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip-relaxed")

    output_handle.close()
    input_handle.close()


    

if __name__ == "__main__" :
    main()
