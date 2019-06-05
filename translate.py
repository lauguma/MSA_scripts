#!/usr/bin/python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import *
from Bio.SeqRecord import SeqRecord
import sys
import copy


# USAGE python translate.py -h


def main():

    parser = ArgumentParser()
    parser.add_argument("-nt", dest="fasta_file", help="fasta file containing nucleotide sequences", type=str)
    parser.add_argument("-aa",dest="file_out", help="fasta file output containing amino acid sequences", type=str)

    args = parser.parse_args()


    fasta = SeqIO.parse(args.fasta_file,"fasta")
    prot = open(args.file_out, "w")

    for s in fasta:
        aas = s.seq.translate(stop_symbol="")
        sequence_object = Seq(str(aas), ExtendedIUPACProtein)
        record = SeqRecord(sequence_object,id=s.id, description="")
        SeqIO.write(record, prot, "fasta")




    fasta.close()
    prot.close()

        


if __name__ == "__main__" :
    main()
