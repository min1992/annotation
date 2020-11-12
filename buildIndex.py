#!/usr/bin/python
# -*- coding=utf8 -*-
##this script is for build index

import sys
import re
refseq=sys.argv[1]
gene=sys.argv[2]

def buildIndex(infile,col,outfile):
    output=open(outfile,"w")
    total=0
    info={}
    chromList=[]
    for line in open(infile,"r"):
        lineinfo=line.strip("\n").split("\t")
        total+=len(line)
        try:
            info[lineinfo[col]]+=len(line)
        except:
            info[lineinfo[col]]=len(line)
        if lineinfo[col] in chromList:
            continue
        else:
            chromList.append(lineinfo[col])
    for chrom in chromList:
        output.write(chrom+"\t"+str(info[chrom])+"\n")
    output.close()

refseqindex=refseq+".fai"
buildIndex(refseq,2,refseqindex)
geneindex=gene+".fai"
buildIndex(gene,1,geneindex)

