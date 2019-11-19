# Scripts

- translate.py

Python script for translating DNA sequences into aminoacid sequences

```
python translate.py -nt dna_sequence.fna -aa protein_sequence.faa
```


- backTranslate.py

Python script for backtranslating protein sequences in fasta format (e.g. multiple alignment sequences) into aligned dna sequences.
It must be run from the directory in which the original DNA sequences and the MSA files are located.


```
python backTranslate.py -nt .ffn -ali .faa.aligned -out .cd.aligned
```

- cat_MSA.py

Python script for concatenating multiple sequence alignments in a single multifasta file alignment.
Extension of the single multiple sequences alignment files must be provided (e.g. .cd.mafft). 
Species name in multifasta files is identified as the string written before the underscore.
E.g. >SCERS288C_15G2550 --> species name is 'SCERS288C' and '15G2550' is the locus tag.


```
python cat_MSA.py -msa .cd.mafft -out cat.cd.mafft
```



- multifasta2ffn.py

This script takes multifasta files containing annotated CDS coming from different genomes as input.
It takes the gene ID and calculate the number of common genes among genomes to create
new fasta files: one fasta file per gene containing all input genomes. 

```
python multifasta2ffn.py -fasta .cds.fna -out .ffn
```

- msa2phy.py

Python script to format conversion for Multiple Sequence Alignments:
from Fasta format to Phylip sequential/interleaved formats (if sequential: set -sequential option). 

```
python msa2phy.py -msa MSA_file.fasta -out outFile.phy -sequential < optional >
```

- blastORFinder.py

Python script to find a coding sequence in a fasta file (e. g. a contigs file) after using blast to figure out its location

Use example:
Imagine that you want to extract a gene sequence after a genome assembly. You are only interested in one gene sequence or in a low number of
gene sequences so it is not necessary a genome annotation for the moment. If you have a gene sequence from a close relative of the organism you
can figure out if the gene is present in your genome and in which position and contig by using blastn. 
e. g.

```
makeblastdb -in contigs.fasta -dbtype nucl -parse_seqids
blastn -db contigs.fa -query query.fa -outfmt 6 -out blast_result.blastn
```

In the file blast_result.blastn you can find the location of that gene if the gene has been found.
You can use this python script to extract the Open Reading Frame corresponding to that location.
If the ORF is found and is complete (It starts by ATG and ends by a stop codon) the sequence will be 
saved in a file called completeCDS_xx.fna. If the ORF is not complete, you will obtain a fasta file named partialCDS_xx.fna

```
python blastORFinder.py -fasta contigs.fasta -blast blast_results.blastn 
```




