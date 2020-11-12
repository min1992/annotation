#!/usr/bin/python
# -*- coding=utf8 -*-
##this script is for annotation

import os
import sys
import re
import math
from baseFunctions import *
import getSequence

STARTCODON="ATG"
ENDCONDON=["TAA","TAG","TGA"]
CONDON2AA={"TT[T|C]":"Phe","TT[A|G]":"Leu","CT.":"Leu","AT[T|C|A]":"Ile","ATG":"Met","GT.":"Val","TC.":"Ser","CC.":"Pro","AC.":"Thr","GC.":"Ala",
"TA[T|C]":"Tyr","CA[T|C]":"His","CA[A|G]":"Gln","AA[T|C]":"Asn","AA[A|G]":"Lys","GA[T|C]":"Asp","GA[A|G]":"Glu","TG[T|C]":"Cys",
"TGG":"Trp","CG.":"Arg","AG[T|C]":"Ser","AG[A|G]":"Arg","GG.":"Gly","TAA":"Ter","TAG":"Ter","TGA":"Ter"}
abbCONDON2AA={"Phe":"F","Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
"Leu":"L","Lys":"K","Met":"M","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}
####main annotation part
class anno(object):
    def __init__(self,fasta,refseq,chrom,start,end,ref,alt,upstream,downstream,splicingSize):
        self.refseq=refseq
        self.fasta=fasta
        self.chrom=chrom
        self.start=start
        self.end=end
        self.ref=ref
        self.alt=alt
        self.upstream=upstream
        self.downstream=downstream
        self.splicingSize=splicingSize

    ###get index info
    def geneIndex(self,gene):
        end=0
        chromInfo={}
        fai=gene+".fai"
        for line in open(fai,"r"):
            lineinfo=line.strip("\n").split("\t")
            chrom=lineinfo[0]
            bytesPos=int(lineinfo[1])
            chromInfo[chrom]={"start":end,"end":end+bytesPos}
            end+=bytesPos
        return chromInfo

    ###get bucket according to index
    def readIn(self,refseq,startByte,endByte):
        readInfo=open(refseq,"r")
        readInfo.seek(startByte,0)
        info=readInfo.read(endByte-startByte)
        readInfo.close()
        return info

    ###judge if two sections are intersectant
    def judgeIntersect(self,start1,start2,end1,end2):
        maxStart=max(int(start1),int(start2))
        minEnd=min(int(end1),int(end2)) 
        if maxStart<=minEnd:
            return True,"middle"
        elif int(start2) > int(start1):
            return False,"right"
        elif int(start2) < int(start1):
            return False,"left"

    ###find neadby multiple transcript accornding one index
    def findnearby(self,lines,index,start,end,indexList,label):
        ###find other index
        if label=="left":
            i=index-1
        elif label=="right":
            i=index+1
        while i >=0:
            line=lines[i]
            lineinfo=line.split("\t")
            stemp=int(lineinfo[4])
            etemp=int(lineinfo[5])
            lastend=int(lines[i-1].split("\t")[5]) if i>0 else 0
            status,orientation=judgeIntersect(stemp,start,etemp,end,lastend)
            if status:
                indexList.append(i)
            else:
                return indexList
            if label=="left":
                i-=1
            elif label=="right":
                i+=1

    ###get cDNA position
    def getCDNA(self,start1,start2,end1,end2):
        status,pos=self.judgeIntersect(start1,start2,end1,end2)
        return status

    ###strive for the remainder
    def striveRemainder(self,num):
        res=int(math.ceil(float(num)/3))
        return res
   
    ##parse aa
    def parseaa(self,temp): 
        condons=list(dict.keys(CONDON2AA))
        indexInfo=[len(re.findall(condon,temp)) for condon in condons]	
        aa=CONDON2AA[condons[indexInfo.index(1)]]
        return aa

    ##get amino acid
    def getaas(self,seq2aa):
        if len(seq2aa)%3 !=0:
            raiseError("SeqError: the sequence needed to translate into aa is wrong!")
        i=0
        aas=""
        while i < len(seq2aa):
            temp=seq2aa[i:i+3]
            aas+=self.parseaa(temp)
            i+=3
        return aas

    ###compare refaa and altaa, exclude same aa in header or tail
    def compareAA(self,refaa,altaa):
        i,leftexcnum,rightexcnum=0,0,0
        length=min(len(refaa),len(altaa))
        while i <length:
            reft=refaa[i:i+3]
            altt=altaa[i:i+3]
            if reft!=altt:
                break
            i+=3
        refaa=refaa[i:]
        altaa=altaa[i:]
        leftexcnum=i
        i=0
        while i < length:
            reft=refaa[-i-3:-i] if i>0 else refaa[-i-3:]
            altt=altaa[-i-3:-i] if i>0 else altaa[-i-3:]
            #print(reft,altt)
            if reft!=altt:
                break
            i+=3
        rightexcnum=i
        refaa=refaa[0:-i] if i>0 else refaa[0:]
        altaa=altaa[0:-i] if i>0 else altaa[0:] 
        return refaa,altaa,leftexcnum,rightexcnum

    ###search stop condon for fs
    def searchSC(self,altCDNA):
        i=0
        length=int(len(altCDNA)/3)*3
        label="?"
        while i < length:
            temp=altCDNA[i:i+3]
            aa=self.parseaa(temp)
            if temp in ENDCONDON:
                label=str(i/3+1)
                break
            i+=3
        return "fs*"+label

    ###parse frameshift
    def parseFS(self,seq2aaref,seq2aaalt,aasref,leftproteinNo,fixedleftCDNA,cdnaseq,downUTR):
        temp="fs"
        pchange,vairiantType="",""
        cdnaseq+=downUTR
        print(fixedleftCDNA,cdnaseq)
        refCDNA=cdnaseq[fixedleftCDNA-1:]
        #print(refCDNA,fixedleftCDNA)
        print(seq2aaref,seq2aaalt,refCDNA)
        altCDNA=re.sub("^"+seq2aaref,seq2aaalt,refCDNA)
        print(refCDNA,altCDNA)
        #print(seq2aaref,seq2aaalt)
        i=0
        while i < min(int(len(refCDNA)/3)*3,int(len(altCDNA)/3)*3):
            temprefaa=self.getaas(refCDNA[i:i+3])
            tempaltaa=self.getaas(altCDNA[i:i+3])
            if temprefaa!=tempaltaa:
                break
            i+=3
        leftproteinNo+=i/3
        refCDNA=refCDNA[i:]
        altCDNA=altCDNA[i:]
        print(i,refCDNA,altCDNA)
        reffirstaa=self.getaas(refCDNA[0:3])
        altfirstaa=self.getaas(altCDNA[0:3])
        if altfirstaa =="Ter":
            pchange="p."+reffirstaa+str(leftproteinNo)+altfirstaa
            vairiantType="nonsense"
        else:
            pchange="p."+reffirstaa+str(leftproteinNo)+altfirstaa
            pchange+=self.searchSC(altCDNA)
            vairiantType="frameshift"
        return pchange,vairiantType
    ###移动蛋白质
    def moveAA(self,aasref,upprotein,downprotein,leftproteinNo,rightproteinNo):
        #print(leftproteinNo,rightproteinNo,aasref,downprotein)
        length=min(len(aasref),len(downprotein))
        i=0
        while i <length:
            aa=aasref[i:i+3]
            daa=downprotein[i:i+3]
            if aa!=daa:
                break
            i+=3
        leftproteinNo+=i/3
        rightproteinNo+=i/3
        upprotein+=aasref[:i]
        aasref=aasref[i:]+downprotein[:i]
        ###dup
        downprotein=downprotein[i:]
        i=0
        while i < len(downprotein):
            daa=downprotein[i:i+len(aasref)]
            if aasref!=daa:
                break
            i+=len(aasref)
        leftproteinNo+=i/3
        rightproteinNo+=i/3
        upprotein+=downprotein[:i]
        downprotein=downprotein[i:]
        ###
        i=0
        while i <length:
            aa=aasref[i:i+3]
            daa=downprotein[i:i+3]
            if aa!=daa:
                break
            i+=3
        leftproteinNo+=i/3
        rightproteinNo+=i/3
        newaasref=aasref[i:]+downprotein[:i]
        upprotein+=aasref[:i]
        downprotein=downprotein[i:]
        #print(leftproteinNo,rightproteinNo)
        return newaasref,leftproteinNo,rightproteinNo,upprotein,downprotein

    ###parse upprotein for nonframeshift protein
    def parseupdownaa(self,upprotein,downprotein,leftproteinNo,rightproteinNo,leftexcnum,rightexcnum,rawaasref,aasref,aasalt):
        #print(leftexcnum,rightexcnum)
        upprotein+=rawaasref[0:leftexcnum]
        downprotein=rawaasref[-rightexcnum:]+downprotein if rightexcnum!=0 else downprotein
        #print(upprotein,downprotein)
        leftproteinNo+=leftexcnum/3
        rightproteinNo-=rightexcnum/3
        #print(aasref,aasalt,leftproteinNo,rightproteinNo)
        if leftproteinNo==rightproteinNo:
            if len(aasalt)==3 and len(aasref)==3:
                temp=""
            elif len(aasalt)>3:
                temp="delins"
            else:
                temp="del"
        #        newdownprotein=rawaasref[-rightexcnum:]+downprotein
#                dup=newdownprotein[0:len(aasref)]
#                if dup==aasref:
#                    leftproteinNo+=len(aasref)/3
                aasref,leftproteinNo,rightproteinNo,upprotein,downprotein=self.moveAA(aasref,upprotein,downprotein,leftproteinNo,rightproteinNo)
            pchange="p."+aasref+str(leftproteinNo)+temp+aasalt
        elif leftproteinNo < rightproteinNo and len(aasalt)>0:
            temp="delins"
#            pchange="p."+upprotein[-3:]+str(leftproteinNo)+"_"+downprotein[0:3]+str(rightproteinNo)+temp+aasalt
            pchange="p."+aasref[0:3]+str(leftproteinNo)+"_"+aasref[-3:]+str(rightproteinNo)+temp+aasalt
        elif rightproteinNo<leftproteinNo:
            temp="ins"
            aasref,rightproteinNo,leftproteinNo,upprotein,downprotein=self.moveAA(aasalt,upprotein,downprotein,rightproteinNo,leftproteinNo)
            dup=upprotein[-len(aasalt):]
            if dup==aasalt:
                temp="dup"
                leftproteinNo=rightproteinNo
                rightproteinNo-=len(aasalt)/3-1
                pchange="p."+aasalt[0:3]+str(rightproteinNo)+"_"+aasalt[-3:]+str(leftproteinNo)+temp if len(aasalt)>3 else "p."+aasalt+str(rightproteinNo)+temp
            else:
                pchange="p."+upprotein[-3:]+str(rightproteinNo)+"_"+downprotein[0:3]+str(leftproteinNo)+temp+aasalt
        else:
            temp="del"
#            pchange="p."+upprotein[-3:]+str(leftproteinNo)+"_"+downprotein[0:3]+temp+aasalt
            ###氨基酸后移
#            newdownprotein=rawaasref[-rightexcnum:]+downprotein
            #print(downprotein)
            #print(aasref)
            aasref,leftproteinNo,rightproteinNo,upprotein,downprotein=self.moveAA(aasref,upprotein,downprotein,leftproteinNo,rightproteinNo)
            pchange="p."+aasref[0:3]+str(leftproteinNo)+"_"+aasref[-3:]+str(rightproteinNo)+temp+aasalt
        return pchange

        

    ##parse c.,if variant locate in coding region,need to judge dup and get c.dup，lef5bp=left2bp+leftlength+lef
    def getcpincoding(self,leftCDNA,rightCDNA,ref,alt,left5bp,right5bp,nearbySize,cdnaseq,seq,downUTR): ###此处左右指上下游，-链left5bp和right5bp要做反向互补处理
        #print(leftCDNA,rightCDNA,seq)
        cchange=""
        pchange=""
        vairiantType=""
        fixedleftCDNA=leftCDNA
        upstart=0
        if leftCDNA-nearbySize-1>=0:
            left5bp=cdnaseq[leftCDNA-nearbySize-1:leftCDNA-1]
        else:
            left5bp=cdnaseq[0:leftCDNA-1]
        if rightCDNA+nearbySize<=len(cdnaseq):
            right5bp=cdnaseq[rightCDNA:rightCDNA+nearbySize]
        else:
            right5bp=cdnaseq[rightCDNA:]
        complement={1:["right",2],2:["both",1],0:["left",2]} ##record how to find nearby bases to form one condon
        leftresidual=int(leftCDNA)%3
        rightresidual=int(rightCDNA)%3
        leftproteinNo=self.striveRemainder(leftCDNA)#int(math.ceil(float(leftCDNA)/3))###start condon encode protein and recode 1,while exclude in protein modification
        rightproteinNo=self.striveRemainder(rightCDNA)#int(math.ceil(float(rightCDNA)/3))
        leftcomType=complement[leftresidual]
        rightcomType=complement[rightresidual]
        #print(leftCDNA,rightCDNA,ref,alt)
        leftbases,rightbases,upprotein,downprotein="","","",""
        if leftcomType[0]!="right":
            leftbases+=left5bp[-leftcomType[1]:]
        #    leftbases=cdnaseq[leftCDNA-leftcomType[1]-1:leftCDNA-1]
            fixedleftCDNA-=leftcomType[1]
            upstart-=leftcomType[1]
            t=int((len(left5bp)-leftcomType[1])/3)*3
#            upprotein=self.getaas(left5bp[-leftcomType[1]-t:-leftcomType[1]])
            upprotein=self.getaas(cdnaseq[0:leftCDNA-1-leftcomType[1]])
            print(fixedleftCDNA,leftCDNA,rightCDNA,upprotein)
        else:
#            t=int(len(left5bp)/3)*3
#            upprotein=self.getaas(left5bp[-t:])
            upprotein=self.getaas(cdnaseq[0:leftCDNA-1])
            print(leftCDNA,rightCDNA,upprotein)
        if rightcomType[0]!="left":
            rightbases+=right5bp[0:rightcomType[1]]
         #   rightbases=cdnaseq[rightCDNA:rightCDNA+rightcomType[1]]
        #    fixedrightCDNA+=rightcomType[1]
            t=int((len(cdnaseq)-rightCDNA-rightcomType[1])/3)*3
#            downprotein=self.getaas(right5bp[rightcomType[1]:rightcomType[1]+t])
            downprotein=self.getaas(cdnaseq[rightCDNA+rightcomType[1]:rightCDNA+rightcomType[1]+t])
            print(leftCDNA,rightCDNA,downprotein)
        else:
            t=int((len(cdnaseq)-rightCDNA)/3)*3
#            downprotein=self.getaas(right5bp[0:t])
            downprotein=self.getaas(cdnaseq[rightCDNA:rightCDNA+t])
            print(leftCDNA,rightCDNA,downprotein)
        seq2aaref=leftbases+ref+rightbases
        seq2aaalt=leftbases+alt+rightbases
        print(seq2aaref,seq2aaalt)
        if ref=="-": ###insertion
            temp="ins"
            #print(seq)
            left5bp+=seq[0]
            right5bp=seq[1]+right5bp
            leftbases+=seq[0]
            rightbases=seq[1]+rightbases
            dup=left5bp[-len(alt):]
            cchange="c."+str(leftCDNA)+"_"+str(rightCDNA)+temp+alt
            ###judge if upstream of ins is same with ins
            if dup==alt:
                temp="dup"
                cchange="c."+str(leftCDNA-len(alt)+1)+"_"+str(leftCDNA)+temp if len(alt)>1 else "c."+str(leftCDNA)+temp
#            elif alt==self.transRevAndComple(ref):
#                temp="inv"
#                cchange="c."+str(leftCDNA)+"_"+str(rightCDNA)+temp+ref if len(alt)>1 else "c."+str(leftCDNA)+temp+ref
            ##p.需要增加上游序列判断dup:nonframe shift
            seq2aaref=leftbases+rightbases
            seq2aaalt=leftbases+alt+rightbases
            #print(seq2aaref,seq2aaalt,leftproteinNo,rightproteinNo)
            aasref=self.getaas(seq2aaref) 
            if len(seq2aaalt)%3==0:#non frameshift
                temp="ins"
                rawaasref=aasref
                proteinNo=leftproteinNo    
                aasalt=self.getaas(seq2aaalt)
                #print(aasref,aasalt)
                aasref,aasalt,leftexcnum,rightexcnum=self.compareAA(aasref,aasalt)
                #print(aasref,aasalt,leftexcnum,rightexcnum)
                pchange=self.parseupdownaa(upprotein,downprotein,leftproteinNo,rightproteinNo,leftexcnum,rightexcnum,rawaasref,aasref,aasalt)
                vairiantType="nonframeshift"
            else: ##frameshift
                pchange,vairiantType=self.parseFS(seq2aaref,seq2aaalt,aasref,leftproteinNo,fixedleftCDNA,cdnaseq,downUTR)
        elif alt=="-" or (len(ref)>1 or len(alt)>1): ###deletion，fs需要在c.和p.之后识别终止密码子位置
            temp,middle="",""
            if alt=="-":
                temp="del"
                cchange="c."+str(leftCDNA)+temp if leftCDNA==rightCDNA else "c."+str(leftCDNA)+"_"+str(rightCDNA)+temp
            else:
                temp="delins"
                cchange="c."+str(leftCDNA)+"_"+str(rightCDNA)+temp+alt if len(ref)>1 else "c."+str(leftCDNA)+temp+alt
                if len(ref)==len(alt):
                    dup=left5bp[-len(alt):]
                    #if dup==alt:
                    #    temp="dup"
                    #    cchange="c."+str(leftCDNA)+"_"+str(rightCDNA)+temp+alt if len(ref)>1 else "c."+str(leftCDNA)+temp+alt
                    if alt==self.transRevAndComple(ref):
                        temp="inv"
                        cchange="c."+str(leftCDNA)+"_"+str(rightCDNA)+temp+ref if len(ref)>1 else "c."+str(leftCDNA)+temp+ref
                middle=alt
            seq2aaalt=leftbases+middle+rightbases
            #print(seq2aaalt)
            aasref=self.getaas(seq2aaref)
            #print(seq2aaref,seq2aaalt)
            if len(seq2aaalt)%3==0:###nonframeshit,remaining bases%3=0
                aasalt=self.getaas(seq2aaalt) 
                rawaasref=aasref
                #print(aasref,aasalt)
                aasref,aasalt,leftexcnum,rightexcnum=self.compareAA(aasref,aasalt)
                #print(leftproteinNo,rightproteinNo,leftexcnum,rightexcnum,aasref,aasalt,downprotein)
                pchange=self.parseupdownaa(upprotein,downprotein,leftproteinNo,rightproteinNo,leftexcnum,rightexcnum,rawaasref,aasref,aasalt)
                vairiantType="nonframe shift"
            else: ###frameshift, del non integral aas
                pchange,vairiantType=self.parseFS(seq2aaref,seq2aaalt,aasref,leftproteinNo,fixedleftCDNA,cdnaseq,downUTR)
        elif len(ref)==1 and len(alt)==1: ##SNV
            temp=">"
            cchange="c."+str(leftCDNA)+ref+temp+alt
            refaa=self.getaas(seq2aaref)
            altaa=self.getaas(seq2aaalt)
            if refaa==altaa: ###sysnoymous
                pchange="p."+refaa+str(leftproteinNo)+"="
                vairiantType="sysnoymous"
            elif refaa != altaa and altaa in ENDCONDON:###nonsense
                pchange="p."+refaa+str(leftproteinNo)+altaa
                vairiantType="nonsense"
            else: ##missense
                pchange="p."+refaa+str(leftproteinNo)+altaa
                vairiantType="missense"
        return cchange,pchange,vairiantType

    ##getcp in intron,增加dup的判断
    def getcpinintron(self,leftCDNA,rightCDNA,s,e,ref,alt,lr,left5bp,right5bp):
        cchange,pchange,vairiantType="","","intronic"
        corr={"left":leftCDNA,"right":rightCDNA}
        updown={"left":"+","right":"-"} ###upstream and downstream
        if (lr[0][1]<=8 and updown[lr[0][0]]=="+") or (lr[1][1]<=8 and updown[lr[1][0]]=="+"):
            vairiantType="splice_region_variant"
        elif (lr[0][1]<=8 and updown[lr[0][0]]=="-") or (lr[1][1]<=8 and updown[lr[1][0]]=="-"):
            vairiantType="splice_region_variant"
        elif (lr[0][1]<=2 and updown[lr[0][0]]=="+") or (lr[1][1]<=2 and updown[lr[1][0]]=="+"):
            vairiantType="splice_donor_variant"
        elif (lr[0][1]<=2 and updown[lr[0][0]]=="-") or (lr[1][1]<=2 and updown[lr[1][0]]=="-"):
                vairiantType="splice_acceptor_variant"
        if ref=="-": ##insersion,
            temp="ins"
            dup=left5bp[-len(alt):]
            cchange="c."+str(corr[lr[0][0]])+updown[lr[0][0]]+str(lr[0][1])+"_"+str(corr[lr[1][0]])+updown[lr[1][0]]+str(lr[1][1])+temp+alt
        elif alt=="-":
            temp="del"
            cchange="c."+str(corr[lr[0][0]])+updown[lr[0][0]]+str(lr[0][1])+"_"+str(corr[lr[1][0]])+updown[lr[1][0]]+str(lr[1][1])+temp if len(ref)>1 else "c."+str(corr[lr[0][0]])+updown[lr[0][0]]+str(lr[0][1])+temp
        elif (len(ref)>1 and len(alt)>1):
            temp="delins"
            cchange="c."+str(corr[lr[0][0]])+updown[lr[0][0]]+str(lr[0][1])+"_"+str(corr[lr[1][0]])+updown[lr[1][0]]+str(lr[1][1])+temp+alt
        else:
            temp=">" 
            cchange="c."+str(corr[lr[0][0]])+str(updown[lr[0][0]])+str(lr[0][1])+ref+temp+alt
        return cchange,pchange,vairiantType

    ##get distance for UTR
    def getDis(self,s,e,strand,position):
        leftborder,rightborder,length,label=position
        UTRlength=abs(leftborder-rightborder)+1 
        if strand=="+":
            if label=="5' UTR":
                leftdis=str(s-rightborder-1)
                rightdis=str(e-rightborder-1)
                vairiantType="5_prime_UTR_variant"
                if abs(s-rightborder-1)>UTRlength:
                    print("variant locate UTR and intergenic region")
                    leftdis="-"+str(UTRlength)
                return leftdis,rightdis,vairiantType
            else:
                leftdis="*"+str(s-leftborder+1)
                rightdis="*"+str(e-leftborder+1)
                vairiantType="3_prime_UTR_variant"
                if abs(e-leftborder+1)>UTRlength:
                    print("variant locate UTR and intergenic region")
                    rightdis="*"+str(UTRlength)
                return leftdis,rightdis,vairiantType
        else:
            if label=="5' UTR":
                leftdis=leftborder-e-1
                rightdis=leftborder-s-1
                vairiantType="5_prime_UTR_variant"
                if abs(leftdis)>UTRlength:
                    print("variant locate UTR and intergenic region")
                    leftdis="-"+str(UTRlength)
                return leftdis,rightdis,vairiantType
            else:
                leftdis="*"+str(rightborder-e+1)
                rightdis="*"+str(rightborder-s+1)
                vairiantType="3_prime_UTR_variant"
                if abs(rightborder-e+1)>UTRlength:
                    print("variant locate UTR and intergenic region")
                    rightdis="*"+str(UTRlength)
                return leftdis,rightdis,vairiantType
                  
    ##get cp in UTR
    def getcpinUTR(self,s,e,ref,alt,strand,position):
        cchange,pchange="",""
        leftdis,rightdis,vairiantType=self.getDis(s,e,strand,position)
        if ref=="-":
            temp="ins"
            cchange="c."+str(leftdis)+"_"+str(rightdis)+temp+alt
        elif alt=="-":
            temp="del"
            cchange="c."+str(leftdis)+"_"+str(rightdis)+temp if leftdis!=rightdis else "c."+str(leftdis)+temp
        elif (len(ref)>1 and len(alt)>1):
            temp="delins"
            cchange="c."+str(leftdis)+"_"+str(rightdis)+temp+alt
        else:
            temp=">"
            cchange="c."+str(leftdis)+ref+temp+alt
        return cchange,pchange,vairiantType

    ##get leftCDNA and rightCDNA so on
    def getexonInfo(self,s,e,strand,position):
        leftborder,rightborder,length,label=position 
        #print(leftborder,rightborder,length,label)
        vairiantType=""
        stemp,etemp=s,e
        upborder,downborder=leftborder,rightborder
        if strand=="+":
            leftCDNA=s-length
            rightCDNA=e-length 
        else:
            leftCDNA=length-e
            rightCDNA=length-s
            stemp,etemp=e,s
            upborder,downborder=rightborder,leftborder
        if abs(stemp-upborder)<=3:
            vairiantType="splice_region_variant"
        elif abs(etemp-rightborder)<=3:
            vairiantType="splice_region_variant"
        elif abs(stemp-upborder) <=2:
            vairiantType="splice_accpetor_variant"
        elif abs(etemp-downborder)<=2:
            vairiantType="splice_donor_variant"
        return leftCDNA,rightCDNA,vairiantType
    ###merge s and e
    def merge(self,lr):
        if lr[0][0]=="both" and lr[1][0]=="both":
            lr[0][0]="left"
            lr[1][0]="left"
        elif lr[1][0]=="both":
            lr[1][0]=lr[0][0]
        elif lr[0][0]=="both":
            lr[1][0]=lr[0][0]
        else:
            pass
        return lr
    ###judge distance
    def judgeD(self,s,e,leftborder,rightborder):
        lr=[]
        for i in [s,e]:
            leftdis=abs(i-leftborder)
            rightdis=abs(i-rightborder)
            if leftdis<rightdis:
                lr.append(["left",leftdis])
            elif leftdis>rightdis:
                lr.append(["right",rightdis])
            else:
                lr.append(["both",leftdis])
        return self.merge(lr)

    ##get leftCDNA and rightCDNA for intron
    def getintronInfo(self,s,e,strand,position):
        leftborder,rightborder,length,label=position
        if strand=="+":
            leftborder-=1
            rightborder+=1
            intronL=(rightborder-leftborder-1)
            leftCDNA=leftborder-length+intronL
            rightCDNA=rightborder-length
            lr=self.judgeD(s,e,leftborder,rightborder)
            return leftCDNA,rightCDNA,lr
        else:
            leftborder-=1
            rightborder+=1
            intronL=(rightborder-leftborder-1)
            leftCDNA=length+intronL-rightborder
            rightCDNA=length-leftborder
            lr=self.judgeD(e,s,rightborder,leftborder)
            return leftCDNA,rightCDNA,lr
            #leftdis=e-rightborder ##upstream distance
            #rigthdis=s-leftborder ##downstram distance
    
    ###get cp cross regions
    def getcpcrossregions(self,positions,s,e,ref,alt,strand):
        temp="del"
        if alt!="-":
            temp="delins"+alt
        symbol=""
        cchange,pchange,vairiantType="","",""
        leftborder1,rightborder1,length1,label1=positions[0]
        leftborder2,rightborder2,length2,label2=positions[1]
        label=label1+"-"+label2
        updis,downdis=0,0
        if strand=="+":
            updis=abs(s-rightborder1)
            downdis=abs(e-leftborder2)
            if "intron" in label1 and "exon" in label2:
                CDNA=e-length2
                ctemp=CDNA-downdis
                stemp="-"
                dis=updis+1
                if abs(s-leftborder1) < abs(s-rightborder1):
                    ctemp-=1
                    stemp="+"
                    dis=abs(s-leftborder1)+1
                cchange="c."+str(ctemp)+stemp+str(dis)+"_"+str(CDNA)+temp
                symbol="-"
            elif "exon" in label1 and "intron" in label2:
                CDNA=s-length1
                ctemp=CDNA+updis
                stemp="+"
                dis=downdis+1
                if abs(e-leftborder2)>abs(e-rightborder2):
                    ctemp+=1
                    stemp="-"
                    dis=abs(e-rightborder2)+1
                cchange="c."+str(CDNA)+"_"+str(ctemp)+str(stemp)+str(dis)+temp
                symbol="+"
            elif ("5' UTR" in label1 or "upstream" in label1) and "exon" in label2:
                CDNA=e-length2
                ctemp=CDNA-downdis
                stemp="-"
                dis=updis+1
                cchange="c."+str(ctemp)+stemp+str(dis)+"_"+str(CDNA)+temp
            elif ("3' UTR" in label2 or "downstream" in label2) and "exon" in label1:
                CDNA=s-length1
                ctemp=CDNA+updis
                stemp="+"
                dis=downdis+1
                cchange="c."+str(CDNA)+"_"+str(ctemp)+stemp+str(dis)+temp
            else:
                pass
        else:
            updis=abs(e-leftborder1)
            downdis=abs(s-rightborder2)
            if "intron" in label1 and "exon" in label2:
                CDNA=length2-s
                ctemp=CDNA-downdis
                stemp="-"
                dis=updis+1
                if abs(e-leftborder1)>abs(e-rightborder1):
                    ctemp-=1
                    stemp="+"
                    dis=abs(e-rightborder1)+1
                cchange="c."+str(ctemp)+stemp+str(dis)+"_"+str(CDNA)+temp
                symbol="-"
            elif "exon" in label1 and "intron" in label2:
                CDNA=length1-e
                ctemp=CDNA+updis
                stemp="+"
                dis=downdis+1
                if abs(s-leftborder2)<abs(s-rightborder2):
                    ctemp+=1
                    stemp="-"
                    dis=abs(s-leftborder2)+1
                cchange="c."+str(CDNA)+"_"+str(ctemp)+stemp+str(dis)+temp
                symbol="+"
            elif ("5' UTR" in label1 or "upstream" in label1) and "exon" in label2:
                CDNA=length2-s
                ctemp=CDNA-downdis
                stemp="-"
                dis=updis+1
                cchange="c."+str(ctemp)+stemp+str(dis)+"_"+str(CDNA)+temp
            elif ("3' UTR" in label2 or "downstream" in label2) and "exon" in label1:
                CDNA=length1-e
                ctemp=CDNA+updis
                stemp="+"
                dis=downdis+1
                cchange="c."+str(CDNA)+"_"+str(ctemp)+stemp+str(dis)+temp
            else:
                pass
        if updis<=2 or downdis<=2 and symbol=="+":
            vairiantType="splice_donor_variant"
        elif updis<=2 or downdis<=2 and symbol=="-":
            vairiantType="splice_acceptor_variant"
        else:
            vairiantType="splice_region_variant"
        return cchange,pchange,vairiantType,label

    ###parse cp for each region
    def getCDNAchange(self,positions,strand,s,e,ref,alt,left5bp,right5bp,nearbySize,cdnaseq,seq,downUTR,splicingSize=2):
        ###judge variant locate upstream or downstream of exon,splicing
        cchange,pchange,vairiantType,label="","","","" 
        if len(positions)<1:
            print("variant is in intergenic!")
            return "","","",""
        elif len(positions)==1:
            label=positions[0][-1] 
            leftCDNA,rightCDNA,rawvairiantType= self.getexonInfo(s,e,strand,positions[0])
            print(s,e,leftCDNA,rightCDNA,ref,alt)
            #print(leftCDNA,rightCDNA)
            if "exon" in label:      
                cchange,pchange,vairiantType=self.getcpincoding(leftCDNA,rightCDNA,ref,alt,left5bp,right5bp,nearbySize,cdnaseq,seq,downUTR)
                if rawvairiantType!="":
                    vairiantType=rawvairiantType
            elif "intron" in label:
                leftCDNA,rightCDNA,lr=self.getintronInfo(s,e,strand,positions[0]) 
                cchange,pchange,vairiantType=self.getcpinintron(leftCDNA,rightCDNA,s,e,ref,alt,lr,left5bp,right5bp)
            elif re.search("UTR",label):
                cchange,pchange,vairiantType=self.getcpinUTR(s,e,ref,alt,strand,positions[0])
            elif re.search("upstream",label) or re.search("downstream",label):
                cchange,pchange,vairiantType="","",label
            else:
                print("variant locate intergenic")
        elif len(positions)==2:
            cchange,pchange,vairiantType,label=self.getcpcrossregions(positions,s,e,ref,alt,strand)
        return cchange,pchange,vairiantType,label
            
    ###trans reverse and complementary
    def transRevAndComple(self,seq): 
        basecomple={"A":"T","G":"C","C":"G","T":"A","-":"-"}
        revcompleSeq=[basecomple[base] for base in seq[::-1]]
        return "".join(revcompleSeq)
    ###locate CDSstart and CDSend
    def locateCDS(self,CDSstart,CDSend,exonstart,exonend):
        #print(exonstart)
        i,j=0,len(exonstart)-1
        while i < len(exonstart):
            status,ori=self.judgeIntersect(CDSstart,exonstart[i],CDSend,exonend[i])
            if status:
                break
            i+=1
        while j >=0:
            status,ori=self.judgeIntersect(CDSstart,exonstart[j],CDSend,exonend[j])
            if status:
                break
            j-=1
        return exonstart[i:j+1],exonend[i:j+1],i,len(exonstart)-1-j
                

    ###split and define function region
    def defineFunction(self,lineinfo,fasta,strand,c,s,e,ref,alt,upstream,downstream,splicingSize):####0-based in refeq file 
        regionInfo={}
        strand=lineinfo[3]
        transtart=int(lineinfo[4])+1
        tranend=int(lineinfo[5])
        cdsstart=int(lineinfo[6])+1
        cdsend=int(lineinfo[7])
        exonstart=map(lambda x:int(x),re.split(",",lineinfo[9])[:-1])
        exonend=map(lambda x: int(x),re.split(",",lineinfo[10])[:-1]) 
        rawtotalExons=int(lineinfo[8])
        exonstart,exonend,leftexons,rightexons=self.locateCDS(cdsstart,cdsend,exonstart,exonend)
        totalExons=len(exonstart)
        if cdsstart==cdsend or re.match("NR",lineinfo[1]):
            print("no coding transcript")
            return "","","","",""
        positions=[]
        leftborder=0 ###left exon border for CDNA position
        rightborder=0 ###right exon border for CDNA position
        length=0 ###length for calculating CDNA change
        cdnaseq=""
        downUTR=""
       # print(s,e,ref,alt)
        s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(c,s,e,fasta,ref,alt,strand)
       # print(s,e,ref,alt)
        if strand=="+":
            if transtart!=cdsstart: 
                status=self.getCDNA(transtart,s,cdsstart-1,e)
                leftborder=transtart
                rightborder=cdsstart-1
                length=0
                if status:
                    positions.append([leftborder,rightborder,length,"5' UTR"])
            if tranend!=cdsend:
                status=self.getCDNA(cdsend+1,s,tranend,e)
                leftborder=cdsend+1
                rightborder=tranend
                length=0
                downUTR=getSeq(fasta,c,cdsend+1,tranend,faiInfo)[0]
                if status:
                    positions.append([leftborder,rightborder,length,"3' UTR"])
            status=self.getCDNA(transtart-2000,s,transtart-1,e)
            leftborder=transtart-2000
            rightborder=transtart-1
            if status:
                positions.append([leftborder,rightborder,length,"upstream"])
            status=self.getCDNA(tranend-1,s,tranend+500,e)
            leftborder=tranend-1
            rightborder=tranend+500
            if status:
                positions.append([leftborder,rightborder,length,"downstream"])
            start=0
            for i in range(totalExons):
                if i==0:
                    leftborder=cdsstart ###left border
                    rightborder=exonend[i]
                    if totalExons==1:
                        rightborder=cdsend
                    start=exonend[i]+1
                    status=self.getCDNA(leftborder,s,rightborder,e)
                    length+=cdsstart-1
                    if status:
                        positions.append([leftborder,rightborder,length,'exon'+str(i+1+leftexons)])
                    cdnaseq+=getSeq(fasta,c,cdsstart,exonend[i],faiInfo)[0]
                else:
                    status=self.getCDNA(start,s,exonstart[i],e)
                    leftborder=start
                    rightborder=exonstart[i]
                    length+=exonstart[i]-start+1
                    if status:
                        positions.append([leftborder,rightborder,length,'intron'+str(i+leftexons)])
                    if i<totalExons-1:
                        status=self.getCDNA(exonstart[i]+1,s,exonend[i],e)
                        leftborder=exonstart[i]+1
                        rightborder=exonend[i]
                        #print(status,leftborder,rightborder,'exon'+str(i+1),s,e)
                        if status:
                            positions.append([leftborder,rightborder,length,'exon'+str(i+1+leftexons)])
                        cdnaseq+=getSeq(fasta,c,exonstart[i]+1,exonend[i],faiInfo)[0]
                    else:
                       status=self.getCDNA(exonstart[i]+1,s,cdsend,e)
                       leftborder=exonstart[i]+1
                       rightborder=cdsend
                       cdnaseq+=getSeq(fasta,c,exonstart[i]+1,cdsend,faiInfo)[0]
                       if status:
                           positions.append([leftborder,rightborder,length,'exon'+str(i+1+leftexons)])
                    start=exonend[i]+1
#                i+=1
        else:
            ref=self.transRevAndComple(ref)
            alt=self.transRevAndComple(alt)
            left5bp=self.transRevAndComple(left5bp)
            right5bp=self.transRevAndComple(right5bp)
            transleft5bp=right5bp
            transright5bp=left5bp
            left5bp=transleft5bp
            right5bp=transright5bp
            seq=self.transRevAndComple(seq)
            if transtart!=cdsstart:
                status=self.getCDNA(transtart,s,cdsstart-1,e)
                leftborder=transtart
                rightborder=cdsstart-1
                length=0
                downUTR=self.transRevAndComple(getSeq(fasta,c,transtart,cdsstart-1,faiInfo)[0])
                if status:
                    positions.append([leftborder,rightborder,length,"3' UTR"])
            if tranend!=cdsend:
                status=self.getCDNA(cdsend+1,s,tranend,e)
                leftborder=cdsend+1
                rightborder=tranend
                length=0
                if status:
                    positions.append([leftborder,rightborder,length,"5' UTR"])
            status=self.getCDNA(tranend+1,s,tranend+2000,e)
            leftborder=tranend+1
            rightborder=tranend+2000
            if status:
                positions.append([leftborder,rightborder,length,"upstream"])
            status=self.getCDNA(transtart-500,s,transtart-1,e)
            leftborder=transtart-500
            rightborder=transtart-1
            if status:
                positions.append([leftborder,rightborder,length,"downstream"])
            end=0
            for i in range(totalExons-1,-1,-1): 
                if i==totalExons-1:
                    leftborder=exonstart[i]+1
                    rightborder=cdsend
                    if totalExons==1:
                        leftborder=cdsstart
                    status=self.getCDNA(leftborder,s,rightborder,e)
                    end=exonstart[i]
                    seqt=getSeq(fasta,c,leftborder,rightborder,faiInfo)[0] 
                    cdnaseq+=self.transRevAndComple(seqt)
                    length=cdsend+1 ###leftCDNA=length-e 
                    #print(leftborder,rightborder,length,"exon"+str(totalExons-i))
                    if status:
                        positions.append([leftborder,rightborder,length,"exon"+str(totalExons-i+rightexons)])
                else:
                    length-=(end-exonend[i])
                    status=self.getCDNA(exonend[i]+1,s,end,e)
                    leftborder=exonend[i]+1
                    rightborder=end
                    if status:
                        positions.append([leftborder,rightborder,length,"intron"+str(totalExons-i-1+rightexons)])
                    if i>0:
                        status=self.getCDNA(exonstart[i]+1,s,exonend[i],e)
                        seqt=getSeq(fasta,c,exonstart[i]+1,exonend[i],faiInfo)[0]
                        cdnaseq+=self.transRevAndComple(seqt)
                        leftborder=exonstart[i]+1
                        rightborder=exonend[i]
                        if status:
                            positions.append([leftborder,rightborder,length,"exon"+str(totalExons-i+rightexons)])
                    else:
                        status=self.getCDNA(cdsstart,s,exonend[i],e) 
                        seqt=getSeq(fasta,c,cdsstart,exonend[i],faiInfo)[0]
                        cdnaseq+=self.transRevAndComple(seqt)
                        leftborder=cdsstart
                        rightborder=exonend[i]
                        if status:
                            positions.append([leftborder,rightborder,length,"exon"+str(totalExons-i+rightexons)])
                    end=exonstart[i] 
        #print(positions)
        cchange,pchange,vairiantType,label=self.getCDNAchange(positions,strand,s,e,ref,alt,left5bp,right5bp,nearbySize,cdnaseq,seq,downUTR,splicingSize)
        return cchange,pchange,vairiantType,label,rawtotalExons


    def main(self):
        chromInfo=self.geneIndex(self.refseq)
        ###readin refseq
        if self.chrom not in dict.keys(chromInfo):
            self.chrom=re.sub("chr","",self.chrom) if re.match("chr",self.chrom) else "chr"+self.chrom
        startByte=chromInfo[self.chrom]['start']
        endByte=chromInfo[self.chrom]['end']
        info=self.readIn(self.refseq,startByte,endByte) 
        ###binary search for transcript index
        lines=info.split("\n")
        indexList=[]
        index=0
        #print(self.chrom,0,len(lines),self.start,self.end,index)
        index=binarySearch(lines,self.chrom,0,len(lines),self.start,self.end,index)
        if index=="none":
            print("no transcript or locate in intergenetic!")
            return [],[]
        indexList.append(index)
        ###find nearby multiple transcripts
        indexList=self.findnearby(lines,index,self.start,self.end,indexList,"left")
        indexList=self.findnearby(lines,index,self.start,self.end,indexList,"right")	
        ###annotation
        annoList=[]
        strandList=[]
        vairiantTypeList=[]
        for i in indexList:
            lineinfo=lines[i].split("\t")
            NCBINo=lineinfo[1]
            EnsemnleNo=lineinfo[-1]
            strand=lineinfo[3]
            cchange,pchange,vairiantType,label,totalExons=self.defineFunction(lineinfo,self.fasta,strand,self.chrom,self.start,self.end,self.ref,self.alt,self.upstream,self.downstream,self.splicingSize)
            hgvs=NCBINo+":"+EnsemnleNo+":"+strand+":"+label+"|"+str(totalExons)+":"+cchange+":"+pchange
            annoList.append(hgvs)
            vairiantTypeList.append(vairiantType)
        return annoList,vairiantTypeList

        ###find upstream and downstream of gene
        #updownIndexList=[]
        #binarySearch(lines,0,len(lines),self.start,self.end,updownIndexList,)
