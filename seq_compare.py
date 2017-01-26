import sys
import os

sys.stderr.write('this is the beginning' + os.linesep)

from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
#only needed for timing
from datetime import datetime

sys.stderr.write('stuff imported: ' + str(datetime.now()) + os.linesep)

def getRange(list_of_positions):
    #get all unique entries by label (skip header)
    entry_list=list(set([x[1] for x in list_of_positions[1:]]))
    rangelist=[]
    for entry in entry_list:
        #take their min/max position and chromosome
        entry_chr=[x[4] for x in list_of_positions if x[1] == entry][0]
        entry_pos=[int(x[5]) for x in list_of_positions if x[1] == entry]
        entry_centerpoint=[int(x[5]) for x in list_of_positions if x[1] == entry and int(x[0]) == 0][0] 
        range_start=min(entry_pos)
        range_stop=max(entry_pos)
        left_dist=entry_centerpoint-range_start
        right_dist=range_stop-entry_centerpoint
        rangelist.append([entry,entry_chr,range_start,range_stop,entry_centerpoint,left_dist,right_dist])
    return rangelist

def getSeq(list_of_ranges,fasta_data):
    new_rangelist=list(list_of_ranges)
    for range in new_rangelist:
        #find that chromosome(it'll be the dict ID in fasta dict) and get the sequence for those positions
        range_seq=fasta_data[range[1]].seq[range[2]-1:range[3]]
        range.append(range_seq)
    return new_rangelist

def getMotif(new_rangelist,motif_len):
    #seqs=[x[5] for x in new_rangelist]
    #cut that out of each sequence
    motif_list=[]
    for entry in new_rangelist:
        #find the middle using distance from centerpoint to left edge and extend motif_len in each direction
        motif_start=entry[5]-int(motif_len)
        motif_stop=entry[5]+int(motif_len)
        #if  we don't have space to grab the motif (i.e. it's right at the end of a sequence), skip it
        #otherwise grab the motif
        if motif_start > 0 and motif_stop < len(entry[7]):
            motif_seq=entry[7][motif_start:motif_stop+1]
            motif_list.append(motif_seq)
    #make motif.  make lc uppercase
    clip_motif=motifs.create([x.upper() for x in motif_list])
    #set pseudocounts
    clip_motif.pseudocounts=0.5
    return clip_motif


def clipMapRanges(clipmap_output):
    trans_list=[]
    transcripts=list(set([x[0] for x in clipmap_output]))
    for trans in transcripts:
        trans_chr=[x[1] for x in clipmap_output if x[0] == trans][0]
        trans_pos=[int(x[2]) for x in clipmap_output if x[0] == trans]
        trans_start=min(trans_pos)
        trans_stop=max(trans_pos)
        trans_list.append([trans,trans_chr,trans_start,trans_stop])
    return trans_list


def findNewSites(clip_motif,fasta_data,clipmap_ranges):
    #use the pssm to search the fasta
    #but not all of the fasta - we only want transcripts, i.e. things with shape/clip data
    #so use the clipmap output to define ranges to search in
    clip_pssm=clip_motif.pssm
    newsites=[]
    for transcript in clipmap_ranges:
        tmpsites=[]
        tx_seq=fasta_data[transcript[1]].seq[transcript[2]-1:transcript[3]]
        for position,score in clip_pssm.search(tx_seq,threshold=3.0):
            tmpsites.append([position,score])
        #regardless of direction, the motif will be found from pos:pos+len(motif).  also translate back to original positions 
        for site in tmpsites:
            #position in the seq object
            subpos=[site[0],site[0]+len(clip_motif)]
            #position in reality - easy if positive
            if subpos[0] >= 0:
                realpos=[subpos[0]+transcript[2],subpos[1]+transcript[2]]
            else:
                #is negative, swap to positive indices
                realpos=[subpos[0]+len(tx_seq)+transcript[2], subpos[1]+len(tx_seq)+transcript[2]]
            #grab the transcript,chr, new positions (+-50 from motif center)
            #should always be odd
            centerpos=int((realpos[0] + realpos[1])/2)
            startpos_50=centerpos-50
            endpos_50=centerpos+50
            newsites.append([transcript[0],transcript[1],startpos_50,endpos_50,centerpos,score])
    return newsites

def addData(newsites,clipmap_outputs):
    biglist=[]
    outlist=[]
    for site in newsites:
        #add clip, shape, local position
        #some will be empty because of THING I NEED TO FIX: defined transcripts as min to max, forgetting about exons.  
        #so some of the ranges won't be in clipshape data)
        clipshape=[[int(x[2])-site[4]] + x for x in clipmap_outputs if x[1] == site[1] and int(x[2]) in range(site[2],site[3]+1)]
        #print len(clipshape)
        #now it needs a rank, and the combo label,to match grab50.py
        #rank based on central clippiness, so save central clippiness
        if len(clipshape) == 101:
            central_clippiness=[x[-1] for x in clipshape if x[0] == 0][0]
            biglist.append([central_clippiness,clipshape])
    #now sort, rank and label
    ranked_biglist=sorted(biglist, key=lambda x: int(x[0]),reverse=True)
    for i in range(0,len(ranked_biglist)):
        for entry in ranked_biglist[i][1]:
            newentry=[[entry[0]] + [str(i) + '_' + entry[1]] + [i] + entry[1:]]
            outlist.append(newentry)
    return outlist

def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

#################################3
#import fasta data
sys.stderr.write('begin data import: ' + str(datetime.now()) + os.linesep)

fasta_dict=SeqIO.index("Mus_musculus.GRCm38.dna_sm.primary_assembly.fa", "fasta", alphabet=IUPAC.unambiguous_dna)

sys.stderr.write('fasta import complete: ' + str(datetime.now()) + os.linesep)

#import the seqs of interest - output of grab50.py
grabbed_positions=[]
with open('testclip2.100.50') as ingrab:
    for line in ingrab:
        grabbed_positions.append(line.strip().split('\t'))

sys.stderr.write('grabbed positions imported: ' + str(datetime.now()) + os.linesep)

#and import the clipmap output
clipmap_positions=[]
with open('testclip2.out') as inclip:
    for line in inclip:
        clipmap_positions.append(line.strip().split('\t'))

sys.stderr.write('clipmap imported: ' + str(datetime.now()) + os.linesep)

###########################
motif_len=int(sys.argv[1])

sys.stderr.write('starting the action: ' + str(datetime.now()) + os.linesep)

ranges=getRange(grabbed_positions)
sys.stderr.write('ranges gotten: ' + str(datetime.now()) + os.linesep)

ranges_withseq=getSeq(ranges,fasta_dict)
sys.stderr.write('sequences gotten: ' + str(datetime.now()) + os.linesep)

clip_motif=getMotif(ranges_withseq,motif_len)
sys.stderr.write('motif computed: ' + str(datetime.now()) + os.linesep)

clipmap_ranges=clipMapRanges(clipmap_positions)
sys.stderr.write('clipmap ranges computed: ' + str(datetime.now()) + os.linesep)

new_sites=findNewSites(clip_motif,fasta_dict,clipmap_ranges)
sys.stderr.write('new sites identififed: ' + str(datetime.now()) + os.linesep)

sites_with_data=addData(new_sites,clipmap_positions)
sys.stderr.write('metadata added: ' + str(datetime.now()) + os.linesep)

outlist=flattenList(sites_with_data)
sys.stderr.write('done! writing final output: ' + str(datetime.now()) + os.linesep)

f=open('test_motif.out','w')
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tOrigPos\tShape\tClip\n')
f.write(outlist)
f.close()

sys.stderr.write('all\'s well that ends well: ' + str(datetime.now()) + os.linesep)
