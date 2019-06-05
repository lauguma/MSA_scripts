# Scripts

- translate.py

Python scripts for translating DNA sequences into aminoacid sequences

```
python translate.py -nt dna_sequence.fna -aa protein_sequence.faa
```


- backTranslate.py

Python scripts for backtranslating protein sequences in fasta format (e.g. multiple alignment sequences) into aligned dna sequences.
It must be run from the directory in which the original DNA sequences and the MSA files are located.


```
python backTranslate.py -nt .ffn -ali .faa.aligned -out .cd.aligned
```


