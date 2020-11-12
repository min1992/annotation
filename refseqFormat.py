#!/usr/bin/python
# -*- coding=utf8 -*-
###this script is used to produce formated refseq file and index
###python version 2.7

import os
import re
import sys
from argparse import ArgumentParser
import gzip

###增加基因+染色体位置+转录本index,先通过区域锁定基因，然后得到索引，根据索引获得转录本及注释
###index    transcript  gene  chromsome   start   end   strand startbytes  endbytes genestart   geneend
def buildIndex(refseq,outfile):
    output=open(outfile,"w")
    startbyte=0
    for line in open(refseq):
        lineinfo=line.strip("\n").split("\t")
        index=lineinfo[0]
        chrom=lineinfo[2]
        gene=lineinfo[12]
        trans=lineinfo[1]
        start=lineinfo[4]
        end=lineinfo[5]
        strand=lineinfo[3]
        endbyte=len(line)+startbyte
        output.write(index+"\t"+trans+"\t"+gene+"\t"+chrom+"\t"+start+"\t"+end+"\t"+strand+"\t"+str(startbyte)+"\t"+str(endbyte)+"\n")
        startbyte+=endbyte
    output.close()
    
        

##从gencode中获取ensemble转录本信息,1-based,no used
def readGFF3(gff3,geneinfo):
    output=open(geneinfo,"w")
    info={} 
    readWay="r"
    if re.search("gz$",gff3):
        readWay="rt"
    with gzip.open(gff3,readWay) as f:
        for line in f:
            lineinfo=line.strip("\n").split("\t")
            if re.match("#",line):
                continue
            genetype=lineinfo[2] 
            if genetype=="gene":
                chrom=lineinfo[0]
                start=lineinfo[3]
                end=lineinfo[4]
                g=re.findall(r"gene_name=(.*?);",lineinfo[8])[0] 
                #info[g]=start+"\t"+end
                output.write(g+"\t"+chrom+"\t"+start+"\t"+end+"\n")
    output.close()


###获取转录本间的对应关系
def getRelations(relation):
    relaDict={}
    readWay="r"
    if "," in relation:
        for infile in re.split(",",relation):
            if re.search("gz$",infile):
                readWay="rt"
                with gzip.open(infile,readWay) as f:
                    for line in f:
                        lineinfo=line.strip("\n").split("\t")
                        ensemble=lineinfo[0]
                        ensembleNonum=re.sub("\.\d","",ensemble)
                        for i in range(1,len(lineinfo)):
                            relaDict[lineinfo[i]]=ensemble
                            temp=re.sub("\.\d","",lineinfo[i])
                            relaDict[temp]=ensembleNonum
    return relaDict


##原refseq文件增加ensemble转录本信息
def readRefseq(refseq,relaDict,outfile):
    output=open(outfile,"w")
    for line in open(refseq):
        lineinfo=line.strip("\n").split("\t")
        if re.match("#",line):
            output.write("\t".join(lineinfo)+"\tENSname\n")
            continue
        tranID=lineinfo[1]
        if tranID in relaDict:
            output.write("\t".join(lineinfo)+"\t"+relaDict[tranID]+"\n")
        else:
            temp=re.sub("\.\d","",tranID)
            if temp in relaDict:
                output.write("\t".join(lineinfo)+"\t"+relaDict[temp]+"\n")
            else:
                print(line)
                output.write("\t".join(lineinfo)+"\tNone\n")
    output.close()

def obtainPars():
    parser=ArgumentParser()
    parser.add_argument("--refseq",dest="refseq",action="store",help="refseq file from UCSC",required=True)
    parser.add_argument("--relation",dest="relation",action="store",help="relation for NCBI transcript and ENSEMBLE transcript",required=True)
    parser.add_argument("--nrefseq",dest="nrefseq",action="store",help="output formated refseq file",required=True)
    parser.add_argument("--index",dest="index",action="store",help="output index  for refseq file",required=True)
    parser.add_argument("--gff3",dest="gff3",action="store",help="NCBI gff3 file",required=True)
    parser.add_argument("--geneinfo",dest="geneinfo",action="store",help="gene info",required=True)
    pars=parser.parse_args()
    return pars

def main():
    pars=obtainPars()
    refseq,relation,nrefseq,index,gff3,geneinfo=pars.refseq,pars.relation,pars.nrefseq,pars.index,pars.gff3,pars.geneinfo
    relaDict=getRelations(relation)
    readRefseq(refseq,relaDict,nrefseq)
    readGFF3(gff3,geneinfo)
    #buildIndex(nrefseq,geneInfo,index)
    


if __name__=="__main__":
    main()




