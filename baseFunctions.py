#!/usr/bin/python
# -*- coding=utf8 -*-

import os
import re
import sys

def raiseError(string):
    sys.exit(string)

def judgeIntersect(start1,start2,end1,end2,lastend):
    maxStart=max(int(start1),int(start2))
    minEnd=min(int(end1),int(end2))
    if maxStart<=minEnd:
        return True,"middle"
    elif int(start2) > int(end1) and int(start2)>= int(lastend):
        return False,"right"
    elif int(start2) > int(end1):
        return False,"left"
    elif int(end2) < int(start1):
        return False,"left"
def mergeOri(oriList):
    left=oriList[0:5]
    right=oriList[-5:]
    leftori=max(left,key=left.count)
    rightori=max(right,key=right.count)
    if leftori==rightori:
        return leftori
    else:
        return False


def binarySearch(lines,chrom,startpos,endpos,start,end,index,upstream=2000,downstream=500): 
    if startpos==endpos or startpos<0:
        return index
    mid=int((startpos+endpos)/2)
#    if(mid==startpos) or (mid==endpos):
#        print("no transcript or locate in intergenetic!")
#        return "none"
    index=mid
#    midline=lines[mid]
#    midlineinfo=midline.split("\t")
    ###由于区间，无法准确通过二分法查找，介入前后五行
#    try:
#        lastlineinfo=lines[mid-1].split("\t")
#        lastend=int(lastlineinfo[5])
#    except:
#        lastend=0
    #print(lastend)
#    if chrom!=midlineinfo[2]:
 #       raiseError("Error: chromsome is wrong!")
    ori=[]
    for i in range(mid-5,mid+5+1):
        if i<1 and i>endpos-1:
            continue
        midline=lines[i]
        midlineinfo=midline.split("\t")
        #print(midlineinfo)
        stemp=int(midlineinfo[4])-upstream
        etemp=int(midlineinfo[5])+downstream
        lastlineinfo=lines[i-1].split("\t")
        lastend=int(lastlineinfo[5])
        nextlineinfo=lines[i-1].split("\t")
        nextend=int(nextlineinfo[5])
#        print(stemp,start,etemp,end,lastend) 
        status,orientation=judgeIntersect(stemp,start,etemp,end,lastend)
#        print(status,orientation)
        if orientation=="middle":
            return i
        if ori!=orientation:
            ori.append(orientation)
    orientation=mergeOri(ori)
    if orientation=="right":
        return binarySearch(lines,chrom,mid,endpos,start,end,index,upstream,downstream)
    elif orientation=="left":
        return binarySearch(lines,chrom,startpos,mid,start,end,index,upstream,downstream)
    else:
        print("no transcript or locate in intergenetic!")
        return "none"

###get start position in fasta
def getBtyes(start,length=60):
    startByte=start+int(start/length)
    return startByte

###get sequence
def getSeq(ref,chr,start,end,faiInfo,nearbySize=2):###1-based
    startBytes=getBtyes(start-nearbySize-1) ###byte positon for variant start
    endBytes=getBtyes(end+nearbySize-1) ###byte position for variant end
    if chr not in dict.keys(faiInfo):
        chr=re.sub("chr","",chr) if re.match("chr",chr) else "chr"+chr
    posBytes=faiInfo[chr] ###byte position for chr
    length=end-start+1  ###variant length
    ###get sequence
    fa=open(ref,"r")
    fa.seek(posBytes+startBytes,0) ###left 2bp
    sequence=fa.read(endBytes-startBytes+1) ###right 2bp
    sequence=re.sub("\n","",sequence) 
    fa.close() 
    return sequence[nearbySize:-(nearbySize)],sequence[:nearbySize],sequence[-(nearbySize):]

