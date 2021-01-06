#!/bin/bash
python main.py --refseq refseq_addEnsemble.sorted.txt --variant $1 --reference Homo_sapiens_assembly19.fasta --spliceSize 2
