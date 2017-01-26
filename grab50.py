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
    better_cliplist=[[x[0],x[1],int(x[2]),x[3],int(x[4])] for x in cliplist]
    enough_stops=sorted([x for x in better_cliplist if int(x[4]) >= min_stops],key=itemgetter(4),reverse=True)
    print "stopcheck complete: %s sites found" % (len(enough_stops))
    for entry in enough_stops:
       #check that we didn't grab it already
       ismatched=len([x for x in blacklist if x[0:2] == entry[0:2] and int(entry[2]) in range(x[2],x[3]+1)])
       print entry
       print ismatched
       if ismatched == 0:
           #grab the entry and +/- 50bp
           start = int(entry[2]) - range_len
           stop = int(entry[2]) + range_len
           #add all those positions to the main list,with a ranked index and a local position
           tmpsites=[[idx]+x for x in better_cliplist if x[0:2] == entry[0:2] and int(x[2]) in range(start,stop+1)]
           for item in tmpsites:
               posidx=int(item[3])-int(entry[2])
               site_list.append([posidx]+[str(item[0])+'_'+item[1]]+item)
           #add the region to a blacklist so it doesn't get grabbed again
           blacklist.append([entry[0],entry[1],start,stop])
           idx+=1
       else:
           print "that site has already been extracted"
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
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tOrigPos\tShape\tClip\n')
f.write(outlist)
f.close()
