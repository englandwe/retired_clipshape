#!/usr/bin/env python
#ENSMUST00000140248	2	90898388	NULL	0
#arguments: clipdata, min stops, desired length, output file
import sys
from collections import Counter
from operator import itemgetter

def getSites(cliplist,min_stops,range_len):
    site_list=[]
    blacklist=[]
    idx=1
    #grab every position with sufficient stops, sort by stops with the largest # first
    better_cliplist=[[x[0],x[1],x[2],int(x[3]),x[4],int(x[5])] for x in cliplist]
    enough_stops=sorted([x for x in better_cliplist if x[5] >= min_stops],key=itemgetter(5),reverse=True)
    for entry in enough_stops:
       #check that we didn't grab it already
       center_idx=better_cliplist.index(entry)
       ismatched=len([x for x in blacklist if x[0:2] == entry[0:2] and center_idx in range(x[2],x[3]+1)])
       if ismatched == 0:
           #grab the entry and +/- 50bp
           #start = int(entry[3]) - range_len
           #stop = int(entry[3]) + range_len
           #no; let's be smart about this and get the +/- 50 list entries, not numerical positions
           #first, find the index of the entry in better_cliplist, which is ordered
           start_idx=center_idx - range_len
           stop_idx=center_idx + range_len

           #add all those positions to the main list,with a ranked index and a local position
           #tmpsites=[[idx]+x for x in better_cliplist if x[0:2] == entry[0:2] and x[3] in range(start,stop+1)]
           
           tmpsites=[[idx]+x for x in better_cliplist[start_idx:stop_idx+1] if x[0] == entry[0]]
           if len(tmpsites) == 2*range_len+1:
               posidx=-range_len        
               for item in tmpsites:
                   site_list.append([posidx]+[str(item[0])+'_'+item[1]]+item)
                   posidx+=1
               #add the region to a blacklist so it doesn't get grabbed again
               blacklist.append([entry[0],entry[1],start_idx,stop_idx])
               idx+=1
    return site_list

def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

################################
clip_data=[]
with open(sys.argv[1]) as infile:
#with open('testclip2.out') as infile:
    for line in infile:
        clip_data.append(line.strip().split('\t'))

min_stops=int(sys.argv[2])
range_len=int(sys.argv[3])

clipout=getSites(clip_data,min_stops,range_len)

outlist=flattenList(clipout)
#print outlist
f=open(sys.argv[4],'w')
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tStrand\tOrigPos\tShape\tClip\n')
f.write(outlist+'\n')
f.close()
