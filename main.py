#!/usr/bin/python
# -*- coding=utf8 -*-
##this script is used to annotate

import os
import sys
import re
#import getSequence
from argparse import ArgumentParser
import annotation
from baseFunctions import *

def getEnd(pos,ref,alt):
    start=int(pos)
    end=int(pos)
    variantType=""
    if ref=="-":###对于输入的插入，在位置之前插入;
        start=pos-1
        end=pos
    elif alt=="-":
        start=pos
        end=start+len(ref)-1
    else:
        start=pos
        end=pos+len(ref)-1
    return start,end

def variantPreProcess(ref,alt,start):
    length=min(len(ref),len(alt))
    i=0
    for i in range(length):
        if alt[i]!=ref[i]:
            return ref[i:],alt[i:],start+i
    return ref[i:],alt[i:],start+i
     

def obtainPars():
    parser=ArgumentParser()
    parser.add_argument("--refseq",dest="refseq",action="store",help="refseq file from UCSC",required=True)
    parser.add_argument("--variant",dest="variant",action="store",help="variant format chr:start:ref:alt",required=True)
    parser.add_argument("--spliceSize",dest="spliceSize",action="store",help="splicing size, default 2",required=True,default=2)
    parser.add_argument("--reference",dest="reference",action="store",help="reference",required=True)
    parser.add_argument("--upstream",dest="upstream",action="store",help="upstream distance,default 2Kb",default=2000)
    parser.add_argument("--downstream",dest="downstream",action="store",help="downstream distance,default 500bp",default=500)
    pars=parser.parse_args()
    return pars

def main():
    pars=obtainPars()
    refseq,variant,spliceSize,reference,upstream,downstream=pars.refseq,pars.variant,pars.spliceSize,pars.reference,pars.upstream,pars.downstream
    variantInfo=variant.split(":")
    if len(variantInfo) != 4:
        raiseError("VariantFormat:")
    ###get variant info
    chrom,pos,ref,alt=variantInfo 
    start,end=getEnd(int(pos),ref,alt)
    ###preprocess variant
    ref,alt,start=variantPreProcess(ref,alt,start)
   # seq,left5bp,right5bp,leftdupList,rightdupList,leftseqList,rightseqList,left5bpList,right5bpList,faiInfo=getSequence.getVariantInfo(chrom,start,end,reference,ref,alt)
    ###move dup along transcript direction
    
    ###annotation
    run=annotation.anno(reference,refseq,chrom,start,end,ref,alt,upstream,downstream,spliceSize)
#    run=annotation.anno(reference,refseq,chrom,start,end,ref,alt,seq,left5bp,right5bp,leftdupList,rightdupList,leftseqList,rightseqList,left5bpList,right5bpList,upstream,downstream,faiInfo,spliceSize)
#    print(variant)
    annoList,vairiantTypeList=run.main()
    print(variant+"\t"+"\t".join(annoList))

if __name__ == "__main__":
    main()

    
    




    
