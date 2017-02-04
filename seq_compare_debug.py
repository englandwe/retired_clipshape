import sys
import os
import random
from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from itertools import groupby
from operator import itemgetter

#only needed for timing
from datetime import datetime

sys.stderr.write('stuff imported: ' + str(datetime.now()) + os.linesep)

#grab the ranges of high-clip positions (grab50 output)
def getRange(list_of_positions):
    #get all unique entries by label (skip header)
    entry_list=list(set([x[1] for x in list_of_positions[1:]]))
    rangelist=[]
    for entry in entry_list:
        #get list of positions w/ chromosome, centerpoint, left/right distances
        entry_sublist=[x for x in list_of_positions if x[1] == entry]
	entry_pos_wchr=[[x[4],x[5],int(x[6])] for x in entry_sublist]
        entry_centerpoint=[int(x[6]) for x in entry_sublist if int(x[0]) == 0][0]
        left_dist=min([int(x[0]) for x in entry_sublist])
        right_dist=max([int(x[0]) for x in entry_sublist])
        #returns a list of list of lists; each sublist is a series of continous positions, each sub-sublist is the chromosome and position
        subrange_list=[]
        for k, g in groupby(enumerate(entry_pos_wchr), lambda (i, x): i-x[2]):
	    subrange_list.append(map(itemgetter(1), g))
        #if there are multiple, the range was split across an exon boundary. need to keep a separate chromosome, range start, and range stop for each in a sublist 
        #if it's a single range, the sublist will just have a single entry
        subtmp=[]
        for subrange in subrange_list:
            entry_chr=subrange[0][0]
            entry_strand=subrange[0][1]
            range_start=min([x[2] for x in subrange])
            range_stop=max([x[2] for x in subrange])
            subtmp.append([entry_chr,entry_strand,range_start,range_stop])
        rangelist.append([entry,subtmp,entry_centerpoint,left_dist,right_dist])
    return rangelist

#get the sequences for high-clip ranges
def getSeq(list_of_ranges,fasta_data):
    new_rangelist=list(list_of_ranges)
    for newrange in new_rangelist:
        tmpseq=''
        for segment in newrange[1]:
            #is segment + or -?
            if segment[1] == '+':
                #find that chromosome(it'll be the dict ID in fasta dict) and get the sequence for those positions
                range_seq=fasta_data[segment[0]].seq[segment[2]-1:segment[3]]
            elif segment[1] == '-':
                range_seq=fasta_data[segment[0]].seq[segment[2]-1:segment[3]].reverse_complement()
            #sequence will be concatentated range by range
            tmpseq=tmpseq + range_seq
        newrange.append(tmpseq)
    return new_rangelist

#find the motif for high-clip ranges
def getMotif(new_rangelist,motif_len):
    motif_list=[]
    for entry in new_rangelist:
        #start at the centerpoint and extend motif_len in each direction
        #to determine where the centerpoint is in the sequence, use the left edge (i.e. if it's -50, go in 51 positions aka seq[50])
        motif_start=abs(entry[3])-int(motif_len)
        motif_stop=abs(entry[3])+int(motif_len)
        #if  we don't have space to grab the motif (i.e. it's right at the end of a sequence), skip it
        #otherwise grab the motif
        if motif_start >= 0 and motif_stop < len(entry[5]):
            motif_seq=entry[5][motif_start:motif_stop+1]
            motif_list.append(motif_seq)
    #make motif.  make lc uppercase
    clip_motif=motifs.create([x.upper() for x in motif_list])
    #set pseudocounts
    clip_motif.pseudocounts=0.5
    return clip_motif

#get ranges for all transcripts in the clipmap output (because it's by position)
def clipMapRanges(clipmap_output):
    trans_list=[]
    transcripts=list(set([x[0] for x in clipmap_output]))
    for trans in transcripts:
        #trans_chr=set([x[1] for x in clipmap_output if x[0] == trans])
        #print trans_chr
        trans_pos_wchr=[[x[1],x[2],int(x[3])] for x in clipmap_output if x[0] == trans]
        #now this is where we need to be careful of nonadjacent ranges
        subrange_list=[]
        for k, g in groupby(enumerate(trans_pos_wchr), lambda (i, x): i-x[2]):
            subrange_list.append(map(itemgetter(1), g))
        subtmp=[]
        for subrange in subrange_list:
            trans_chr=subrange[0][0]
            trans_strand=subrange[0][1]
            trans_start=min([x[2] for x in subrange])
            trans_stop=max([x[2] for x in subrange])
            subtmp.append([trans_chr,trans_strand,trans_start,trans_stop])
        trans_list.append([trans,subtmp])
    return trans_list

#find sites matching the motif
def findNewSites(clip_motif,fasta_data,clipmap_ranges):
    #use the pssm to search the fasta
    #but not all of the fasta - we only want transcripts, i.e. things with shape/clip data
    #so use the clipmap output to define ranges to search in
    clip_pssm=clip_motif.pssm
    newsites=[]
    for transcript in clipmap_ranges:
        #tx_seq=fasta_data[transcript[1]].seq[transcript[2]-1:transcript[3]]
        tx_seq=''
        for segment in transcript[1]:
            #is segment + or -?
            if segment[1] == '+':
                #find that chromosome(it'll be the dict ID in fasta dict) and get the sequence for those positions
                range_seq=fasta_data[segment[0]].seq[segment[2]-1:segment[3]]
            elif segment[1] == '-':
                range_seq=fasta_data[segment[0]].seq[segment[2]-1:segment[3]].reverse_complement()
            #sequence will be concatentated range by range
            tx_seq=tx_seq + range_seq
        tmpsites=[]
        for position,score in clip_pssm.search(tx_seq,threshold=3.0):
            tmpsites.append([position,score])
        #build a list (dict?) of seq positions to real positions:
        #if transcript is +, realpos is direct.  just take the base #s from start to end.
        #ex:  exon1 50-100 exon2 200-300 exon3 500-1000
        #you'd get seq1 = 50, seq2 = 51...seqx=100 seqy=200...seqz=300, seqi=500

        #if transcript is negative, you have to start counting from the RIGHT SIDE of the exon
        #ex:  exon1 50-100 exon2 200-300 exon3 500-1000
        #you'd get seq1 = 100, seq2 = 99...seqx=50 seqy=300...seqz=200, seqi=1000 .. seqj=500

        #first, make the ranges.  this is the same regardless of strand 
        #(the first x positions are always from exon1, even if they're in reverse)
        startval=1
        tx_ranges=[]
        for idx,subset in enumerate(transcript[1]):
            range_seq_start=startval
            range_seq_stop=startval+(subset[3]-subset[2])
            startval+=range_seq_stop
            tx_ranges.append([idx,subset[0],subset[1],range_seq_start,range_seq_stop])
        #ok, now we can look inside the appropriate
        transcript_positions={}
        for i in range(1,len(tx_seq)+1):
            #find the appropriate transcript range
            proper_range=[x for x in tx_ranges if i in range(x[3],x[4]+1)][0]
            #if positive
            if proper_range[2] == '+':
                #find your number's position in the proper range
                range_idx=[x[0] for x in list(enumerate(range(proper_range[3],proper_range[4]+1))) if x[1] == i][0]
                #get the real position at that index, grab chr & pos.  strandedness no longer matters at this point.
                realpos=[proper_range[1],range(transcript[1][proper_range[0]][2],transcript[1][proper_range[0]][3]+1)[range_idx]]
                transcript_positions[i] = realpos
            elif proper_range[2] == '-':
                #find your number's position in the proper range; but now everything's reversed.  the range therefore counts down
                range_idx=[x[0] for x in list(enumerate(range(proper_range[4],proper_range[3]-1,-1))) if x[1] == i][0]
                #get the real position at that index, grab chr & pos.  strandedness no longer matters at this point.
                realpos=[proper_range[1],range(transcript[1][proper_range[0]][2],transcript[1][proper_range[0]][3]+1)[range_idx]]
                transcript_positions[i] = realpos
        #general useful note: regardless of direction, the motif will be found from pos:pos+len(motif).  also translate back to original positions 
        #we don't actually want matches on the other strand, because these are (presumably single-stranded) transcripts.
        #neg-strand transcripts have already been RC'd so they will be handled correctly
        tmpsites_pos=[x for x in tmpsites if x[0] >= 1]
        for site in tmpsites_pos:
            #position in the seq object
            subpos=[site[0],site[0]+len(clip_motif)]
                #so this bit is tricky
                #I have the positions in the sequence above
                # if the sequence is contiguous, this would be easy; just add the sequence pos to the tx 
                #but they aren't; need to account for cross-splice motifs
                #end result: need txid, chr, centerpos of motif,positions 50bp up, psoition 50bp down
                #previously it was ranges; but it's probably easier to build a positions list now, then change add_data_fast
            #step one: find centerpos, start and stop pos based on the sequence
            centerpos=int((subpos[0] + subpos[1])/2)
            startpos_50=centerpos-50
            endpos_50=centerpos+50
            #we want to make sure we're not too close to the end of the tx, so startpos <1 or endpos > the numer of positions are right out
            if startpos_50 >= 1 and endpos_50 <=len(tx_seq):
                #loop from start to endpos, finding [trans,chr,realpos]
                #strand is no longer important
                posout=[transcript_positions[pos] for pos in range(startpos_50, endpos_50+1)]
                real_centerpos=transcript_positions[centerpos]
                newsites.append([transcript[0],posout,real_centerpos,score])
    return newsites

#grab new matching sites with low clip values
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
        positions_with_clipshape=[]
        #let's check for a full length region (101bp in this case)
            if len(site[1]) == 101:
                #before we even bother with a site, let's see if it has a central clippiness of zero
                central_chrom,central_position=site[2]
                central_clippiness=[x[4] for x in clipmap_outputs if x[0] == site[0] and x[1] == central_chrom and x[2] == central_position]
                if len(central_clippiness) > 0 and int(central_clippiness[0]) == 0:
                #grab clip/shape data for all positions in the site
                    for chrom, position in site[1]:
                        #get the clip/shape data
                        pos_clipshape = [[chrom,position,x[3],x[4]] for x in clipmap_outputs if x[0] == site[0] and x[1] == chrom and x[2] == position]
                        if len(pos_clipshape) > 0:
                            positions_with_clipshape.append(pos_clipshape[0])
                    #check for shape/clip data for all positions.  only continue if no empty clipshape, and no NULL in teh shape data
                    if len([x for x in positions_with_clipshape if len(x) == 0]) == 0 and len([x for x in positions_with_clipshape if 'NULL' in x]) == 0:
                        #all basic checks passed.  add local position, 'rank', and label
                        local_start=int(-50)
                        mock_rank=1
                        for item in positions_with_clipshape:
                            fmt_site=[local_start,str(mock_rank)+'_'+site[0],mock_rank,site[0],item[0],item[1],item[2],item[3]]
                            biglist.append(fmt_site)
        if len(biglist) >= sampling_depth:
            break
    return biglist

#what it says on the tin
def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

#################################3
#import fasta data
sys.stderr.write('begin data import: ' + str(datetime.now()) + os.linesep)

fasta_dict=SeqIO.index("../Mus_musculus.GRCm38.dna_sm.primary_assembly.fa", "fasta", alphabet=IUPAC.unambiguous_dna)

sys.stderr.write('fasta import complete: ' + str(datetime.now()) + os.linesep)

#import the seqs of interest - output of grab50.py
grabbed_positions=[]
with open('../fake_testclip2.100.50') as ingrab:
    for line in ingrab:
        grabbed_positions.append(line.strip().split('\t'))

sys.stderr.write('grabbed positions imported: ' + str(datetime.now()) + os.linesep)

#and import the clipmap output
clipmap_positions=[]
with open('../testclip2.out.cut') as inclip:
    for line in inclip:
        clipmap_positions.append(line.strip().split('\t'))

sys.stderr.write('clipmap imported: ' + str(datetime.now()) + os.linesep)

rpkmdict={}
with open('SRR1534954.rpkm') as infile:
    for line in infile:
        inline=line.strip().split('\t')
        if '#' not in inline[0]:
            rpkmdict[inline[0]] = float(inline[4])


###########################
motif_len=int(5)

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

###intermediate output for new sites###
tmpsites=flattenList(new_sites)
tmpout=open('intermed_newsites','w')
tmpout.write(tmpsites)
tmpout.close()
####

samp_count=int(100)

#filter for RPKM >10, and output the list
filt_newsites=[newsite for newsite in new_sites if rpkmdict[newsite[0]] > 10]

sites_with_data=addDataFast(filt_newsites,clipmap_positions,samp_count)
sys.stderr.write('metadata added: ' + str(datetime.now()) + os.linesep)

outlist=flattenList(sites_with_data)
sys.stderr.write('done! writing final output: ' + str(datetime.now()) + os.linesep)

f=open('test_motif.out','w')
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tOrigPos\tShape\tClip\n')
f.write(outlist)
f.close()

sys.stderr.write('all\'s well that ends well: ' + str(datetime.now()) + os.linesep)
