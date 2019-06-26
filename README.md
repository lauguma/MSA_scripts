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
from Fasta format to Phylip sequential format

```
python msa2phy.py -msa MSA_file.fasta -out outFile.phy
```




