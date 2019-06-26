#!/usr/bin/python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio import AlignIO

# Function to convert MSA in fasta format into Phylip
def fasta2phy(msa_input, phy_out):
    input_handle = open(msa_input, "rU")
    output_handle = open(phy_out, "w")
    seqs = []
    headers = []
    alignments = SeqIO.parse(input_handle, "fasta")
    for a in alignments:
        headers.append(str(a.id.split("_")[0]))
        seqs.append(str(a.seq))
    input_handle.close()

    output_handle.write("  " + str(len(headers)) + "  " + str(len(seqs[0])) + "  " + "\n")

    for x in range(0, len(headers)):
        output_handle.write(headers[x] + "  " + seqs[x] + "\n")
    output_handle.close()


# USAGE python msa2phy.py -msa < MSA_file.fasta > -out < outFile.phy >


def main():

    parser = ArgumentParser()
    parser.add_argument("-msa", dest="msa_input", \
                        help="fasta file containing MSA in fasta format", type=str)
    parser.add_argument("-out",dest="phy_out",  \
                        help="File output name where msa will be written in Phylip format", type=str)
    parser.add_argument("-sequential", dest="phylip_seq", action="store_true", \
                        help="Set this option to convert fasta file into Phylip sequential format otherwise it will \
                             be written in Phylip interleaved format")

    args = parser.parse_args()


    input_handle = open(args.msa_input, "rU")
    output_handle = open(args.phy_out, "w")


    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip-relaxed")

    output_handle.close()
    input_handle.close()


    if args.phylip_seq == True:
        fasta2phy(args.msa_input, args.phy_out)


    

if __name__ == "__main__" :
    main()
