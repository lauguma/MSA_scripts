#!/usr/bin/python

import argparse
from Bio import SeqIO
import sys
import numpy as np
import os


def read_seq(fasta,header):
    """Extracts sequence from multifasta file by giving the seq.ID"""
    fna = SeqIO.parse(fasta,"fasta")
    for f in fna:
        if f.id == header:
            return f.seq

def pairs_comparisons(fasta):
    """Reads fasta MSA and calculates all pairwise comparisons."""
    headers = []
    fna = SeqIO.parse(fasta,"fasta")
    for f in fna:
        headers.append(f.id)
    tuples = [(x,y) for x in headers for y in headers if x != y]
    for entry in tuples:
        if (entry[1], entry[0]) in tuples:
            tuples.remove((entry[1],entry[0]))
    
    return tuples

def poisson_dist(seq1,seq2):
    """Calculates distances between 2 sequences"""
    i = 0
    dist = 0
    gap = 0
    gap_nulo = 0
    while i <= len(seq1)-1:
        aa1 = seq1[i]
        aa2 = seq2[i]
        if aa1 != aa2 and aa1 != "-" and aa2 != "-":
            dist += 1
        elif aa1 == "-" or aa2 == "-":
            gap +=0
        i += 1
    # Apply a correction by sequence length
    pdist = poisson(dist, (len(seq1)-gap))
    return pdist

def poisson(distance, length):
    """Corrects distance using log and sequence length"""
    calc = 1-(float(distance)/float(length))
    if distance == 0 :
        poisson_d = 0
    else:
        poisson_d = - np.log(calc)
    return poisson_d


def main():
    """Pairwise distances calculator

    This scripts takes as input a set of multiple sequence alignments
    with the same extension in fasta format.
    
    The output consits of a set tab-separated files, one per each
    pairwise genome comparison. Each tab-file contains a gene name
    per row and its distance value.

    This script requires that `Biopython` be installed within the Python(3)
    environment you are running this script in.
    """

    parser = argparse.ArgumentParser(description=main.__doc__, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-suffix", dest="suffix_msa", \
                        help="extension of nucleotide sequences", type=str)
    args = parser.parse_args()

    path = os.getcwd()

    aln_files = [f for f in os.listdir(path) if f.endswith(args.suffix_msa)]
    
    fasta = aln_files[0]
    pairs = pairs_comparisons(fasta)
    print("Calculating distances between the following sequences...")
    for p in pairs:
        print(p)

    for f in aln_files:
        gene = f.replace(args.suffix_msa,"")
        for p in pairs:
            comparison = str(p[0])+"_"+str(p[1])
            f_out = open(comparison +".txt","a")
            sec1 = read_seq(f,p[0])
            sec2 = read_seq(f,p[1])
            distance = poisson_dist(sec1,sec2)
            f_out.write(gene+"\t"+str(distance)+"\n")

if __name__ == "__main__" :
    main()
