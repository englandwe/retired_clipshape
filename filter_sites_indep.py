##Usage: filter_sites_indep.py clipmap_output rpkm_file newsites_file outfile samp_count
import sys
import os
import random
from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
#only needed for timing
from datetime import datetime

#########Functions###########

def addDataFast(newsites,clipmap_outputs,sampling_depth):
    #ALTERNATIVE IDEA: randomly subsample X newsites, find their central clippiness, if it's 0, keep it.  if not, toss.
    #continue sampling until you have enough true-zero newsites
    #in this test case, there were 141, so let's try for that
    biglist=[]
    outlist=[]
    sites_ret=0
    #shuffle list for randomness
    random.shuffle(newsites)
    for site in newsites:
        #grab clip/shape data for all positions in the site and add a local position
        clipshape=[[int(x[2])-int(site[4])] + x for x in clipmap_outputs if x[1] == site[1] and int(x[2]) in range(int(site[2]),int(site[3])+1)]
        #now it needs a rank, and the combo label,to match grab50.py
        #rank based on central clippiness, so save central clippiness
        #second part of this if statement will become unnecessary when transcripts are fixed
        #third part filters out any sites with null shape data anywhere in the range
        if len(clipshape) == 101 and 0 in [x[0] for x in clipshape] and 'NULL' not in [x[4] for x in clipshape]:
            central_clippiness=[x[-1] for x in clipshape if x[0] == 0][0]
            if float(central_clippiness) == 0:
                clipcombo=[central_clippiness,clipshape]
                biglist.append(clipcombo)
        if len(biglist) >= sampling_depth:
            break
    #now sort, rank and label
    #this is a bit silly with all zeroes, but i'll keep it in case we ever want higher vals in future
    #plus they do need some sort of identifier to append to the transcript name in case of multiple motifs per transcript
    ranked_biglist=sorted(biglist, key=lambda x: int(x[0]),reverse=True)
    for i in range(0,len(ranked_biglist)):
        for entry in ranked_biglist[i][1]:
            newentry=[entry[0]] + [str(i) + '_' + entry[1]] + [i] + entry[1:]
            outlist.append(newentry)
#            sys.stdout.write(newentry + os.linesep)
    return outlist


def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

#################################3
sys.stderr.write('starting the action: ' + str(datetime.now()) + os.linesep)

#and import the clipmap output
clipmap_positions=[]
#with open('testclip2.out') as inclip:
with open(sys.argv[1]) as inclip:
    for line in inclip:
        clipmap_positions.append(line.strip().split('\t'))

sys.stderr.write('clipmap imported: ' + str(datetime.now()) + os.linesep)

###########################
#motif_len=int(sys.argv[1])

#####################################
#Because picking out the data takes forever, filter by rpkm before adding data
#eventually do this even earlier to speed the process up

rpkmdict={}
#with open('SRR1534954.rpkm') as infile:
with open(sys.argv[2]) as infile:
    for line in infile:
        inline=line.strip().split('\t')
        if '#' not in inline[0]:
            rpkmdict[inline[0]] = float(inline[4])
#read in new sites
newlist=[]
#with open('intermed_newsites') as infile:
with open(sys.argv[3]) as infile:
    for line in infile:
        newlist.append(line.strip().split('\t'))

sys.stderr.write('rpkm and newsites read: ' + str(datetime.now()) + os.linesep)

samp_count=int(sys.argv[5])

#filter for RPKM >10, and output the list
filt_newsites=[newsite for newsite in newlist if rpkmdict[newsite[0]] > 10]

sys.stderr.write('newsites filtered: ' + str(datetime.now()) + os.linesep)

#filt_newsites_flat=flattenList(filt_newsites)
#newout=open('intermed_newsites_filtered','w')
#newout.write(filt_newsites_flat)
#newout.close()

#add the data
sites_with_data=addDataFast(filt_newsites,clipmap_positions,samp_count)
sys.stderr.write('metadata added: ' + str(datetime.now()) + os.linesep)

outlist=flattenList(sites_with_data)
sys.stderr.write('done! writing final output: ' + str(datetime.now()) + os.linesep)

f=open(sys.argv[4],'w')
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tOrigPos\tShape\tClip\n')
f.write(outlist)
f.close()

sys.stderr.write('all\'s well that ends well: ' + str(datetime.now()) + os.linesep)
