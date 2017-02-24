#truerandom.py clipmap_output desired_sites range_len zero_or_any outfile
import sys
import random
import copy

def getRandomSites(cliplist,shuffled_cliplist,zeros_only,site_count,range_len,multiplier):
    site_list=[]
    blacklist=[]
    idx=1
    #cliplists are large, and we are lazy.  just grab 10x site_count's worth of entries, which chould be more than enough.  if we need more, we can get them.
    #if not, we can always grab more.
#    sys.stderr.write('shuffling'+ '\n')
#    random.shuffle(cliplist)
#    sys.stderr.write('shuffling complete.  subsampling cliplist'+ '\n')
    this_many=multiplier*site_count
    #better_cliplist=[[x[0],x[1],x[2],int(x[3]),x[4],int(x[5])] for x in cliplist]
    #grab every position with zero stops, or keep all of them
    if zeros_only=='zero':
        proper_stops=[]
        stopidx=0
        while len(proper_stops) < this_many:
            if int(shuffled_cliplist[stopidx][5]) == 0:
                proper_stops.append(shuffled_cliplist[stopidx])
            stopidx+=1
#        proper_stops=[x for x in better_cliplist if x[5] == 0]
    elif zeros_only=='any':
        proper_stops=shuffled_cliplist[0:this_many+1]
    else:
        failmsg="%s is not a valid option; please use 'zero' or 'any'" % (zeros_only)
        return failmsg
    sys.stderr.write('looping over entries' + '\n')
    for entry in proper_stops:
       #check that we didn't grab it already
       center_idx=cliplist.index(entry)
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
           
           tmpsites=[[idx]+x for x in cliplist[start_idx:stop_idx+1] if x[0] == entry[0]]
           if len(tmpsites) == 2*range_len+1:
               posidx=-range_len        
               for item in tmpsites:
                   site_list.append([posidx]+[str(item[0])+'_'+item[1]]+item)
                   posidx+=1
               #add the region to a blacklist so it doesn't get grabbed again
               blacklist.append([entry[0],entry[1],start_idx,stop_idx])
               idx+=1
               if idx > site_count:
                   break
    return site_list


def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

######
#import clipmap
clipmap_positions=[]
with open(sys.argv[1]) as inclip:
#with open('../debug_pen/clipmap-invitro-v2.out') as inclip:
    for line in inclip:
        clipmap_positions.append(line.strip().split('\t'))
clipmap_shuf=copy.deepcopy(clipmap_positions)
random.shuffle(clipmap_shuf)

#other vars
site_count=int(sys.argv[2])
range_len=int(sys.argv[3])
zeros_only=str(sys.argv[4])
#site_count=int(140)
#range_len=int(50)
#zeros_only=str('zero')
multiplier=10

randos=getRandomSites(clipmap_positions,clipmap_shuf,zeros_only,site_count,range_len,multiplier)
if type(randos) == str:
    print randos
else:
    print len(randos)
    while len(randos) < site_count:
        multiplier+=5
        randos=getRandomSites(clipmap_positions,zeros_only,site_count,range_len,multiplier)
    outlist=flattenList(randos)
#    f=open('testout','w')
    f=open(sys.argv[5],'w')
    f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tStrand\tOrigPos\tShape\tClip\n')
    f.write(outlist+'\n')
    f.close()

