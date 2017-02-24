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
        print trans_chr
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

def importClipRanges(cliprange_list):
    trans_list=[]
    transcripts=list(set([x[0] for x in cliprange_list]))
    for trans in transcripts:
        subtmp=[x[1:] for x in cliprange_list if x[0] == trans]
        trans_list.append([trans,subtmp])
    return trans_list

def findNewSitesFast(transcript,clip_motif,clip_pssm,fasta_data,clipmap_outputs,rpkmdict,mock_rank):
    newsites=[] 
    #we don't want none unless you got 10RPKM hon                          
    if rpkmdict[transcript[0]] > 10:
        #tx_seq=fasta_data[transcript[1]].seq[transcript[2]-1:transcript[3]]
        print 'getting tx data'
        tx_seq=''
        for segment in transcript[1]:
            #is segment + or -?
            if segment[1] == '+':
                #find that chromosome(it'll be the dict ID in fasta dict) and get the sequence for those positions
                range_seq=fasta_data[segment[0]].seq[int(segment[2])-1:int(segment[3])]
            elif segment[1] == '-':
                range_seq=fasta_data[segment[0]].seq[int(segment[2])-1:int(segment[3])].reverse_complement()
            else:
                print "WHAT???? %s" % segment 
            #sequence will be concatentated range by range
            tx_seq=tx_seq + range_seq
        #see if there are any hits
        print 'finding hits'
        tmpsites=[]
        for position,score in clip_pssm.search(tx_seq,threshold=3.0):
            tmpsites.append([position,score])
        #pointless to go further if there are no appropriate hits
        tmpsites_pos=[x for x in tmpsites if x[0] >= 1]
        if len(tmpsites_pos) > 0:
            print 'got this many hits: %s' % (len(tmpsites_pos))
            #build a list (dict?) of seq positions to real positions:
            #if transcript is +, realpos is direct.  just take the base #s from start to end.
            #ex:  exon1 50-100 exon2 200-300 exon3 500-1000
            #you'd get seq1 = 50, seq2 = 51...seqx=100 seqy=200...seqz=300, seqi=500
    
            #if transcript is negative, you have to start counting from the RIGHT SIDE of the exon
            #ex:  exon1 50-100 exon2 200-300 exon3 500-1000
            #you'd get seq1 = 100, seq2 = 99...seqx=50 seqy=300...seqz=200, seqi=1000 .. seqj=500
    
            #first, make the ranges.  this is the same regardless of strand 
            #(the first x positions are always from exon1, even if they're in reverse)
            print 'getting txranges'
            startval=1
            tx_ranges=[]
            for idx,subset in enumerate(transcript[1]):
                range_seq_start=startval
                range_seq_stop=startval+(int(subset[3])-int(subset[2]))
                startval=range_seq_stop+1
                tx_ranges.append([idx,subset[0],subset[1],range_seq_start,range_seq_stop])           
            #ok, now we can look inside the appropriate
            print 'making transcript pos dict'
            transcript_positions={}
            for i in range(1,len(tx_seq)+1):
                #find the appropriate transcript range
                proper_range=[x for x in tx_ranges if i in range(x[3],x[4]+1)][0]
                #if positive
                if proper_range[2] == '+':
                    #find your number's position in the proper range
                    range_idx=[x[0] for x in list(enumerate(range(proper_range[3],proper_range[4]+1))) if x[1] == i][0]
                    #get the real position at that index, grab chr & pos & strandedness. 
                    realpos=[proper_range[1],proper_range[2],range(int(transcript[1][proper_range[0]][2]),int(transcript[1][proper_range[0]][3])+1)[range_idx]]
                    transcript_positions[i] = realpos
                elif proper_range[2] == '-':
                    #find your number's position in the proper range; but now everything's reversed.  the range therefore counts down
                    range_idx=[x[0] for x in list(enumerate(range(proper_range[4],proper_range[3]-1,-1))) if x[1] == i][0]
                    #get the real position at that index, grab chr & pos & strandedness.
                    realpos=[proper_range[1],proper_range[2],range(int(transcript[1][proper_range[0]][2]),int(transcript[1][proper_range[0]][3])+1)[range_idx]]
                    transcript_positions[i] = realpos
            #general useful note: regardless of direction, the motif will be found from pos:pos+len(motif).  also translate back to original positions 
            #we don't actually want matches on the other strand, because these are (presumably single-stranded) transcripts.
            #neg-strand transcripts have already been RC'd so they will be handled correctly
            print 'cycling through potential sites'
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
                    print 'ok, length check passed'
                    #loop from start to endpos, finding [trans,chr,realpos]
                    #strand is no longer important
                    posout=[transcript_positions[pos] for pos in range(startpos_50, endpos_50+1)]
                    real_centerpos=transcript_positions[centerpos]
                    potential_newsite=[transcript[0],posout,real_centerpos,score]
##################integrate add_data here
                    #let's check for a full length region (101bp in this case)
                    if len(potential_newsite[1]) == 101:
                        print 'calculating central clippiness'
                    #before we even bother with a site, let's see if it has a central clippiness of zero
                        central_chrom,central_strand,central_position=potential_newsite[2]
                        #central_clippiness=[x[4] for x in clipmap_outputs if x[0] == potential_newsite[0] and x[1] == central_chrom and x[3] == central_position]
                        #with a generator
                        clipgen=(x for x in clipmap_outputs if x[0] == potential_newsite[0] and x[1] == central_chrom and int(x[3]) == central_position)
                        clipgen_search=next(clipgen,None)
                        if clipgen_search is not None:
                            central_clippiness=int(clipgen_search[5])
                        else:
                           central_clippiness=None
                        print central_clippiness
                        if central_clippiness == 0:
                            #grab clip/shape data for all positions in the site
                            #base_posgen=(x for x in clipmap_outputs if x[0] == potential_newsite[0])
                            base_list=[x for x in clipmap_outputs if x[0] == potential_newsite[0]]
                            positions_with_clipshape=[]
                            for chrom, strand, position in potential_newsite[1]:
                                #get the clip/shape data
                                #generator here too
                                #pos_clipshape = [[chrom,position,x[4],x[5]] for x in clipmap_outputs if x[0] == potential_newsite[0] and x[1] == chrom and x[3] == position]
                                print 'getting clip/shape data'
                                for entry in base_list:
                                    if entry[1] == chrom and int(entry[3]) == position:
                                        positions_with_clipshape.append([chrom,position,entry[4],entry[5]])
                                #posgen=([chrom,position,x[4],x[5]] for x in clipmap_outputs if x[0] == potential_newsite[0] and x[1] == chrom and int(x[3]) == position)
                                #posgen_search=next(posgen,None)
                                #if posgen_search is not None:
                                #    positions_with_clipshape.append(posgen_search)
                            #check for shape/clip data for all positions.  only continue if no empty clipshape, and no NULL in teh shape data
                            print 'checking we have complete data'    
                            if len([x for x in positions_with_clipshape if len(x) == 0]) == 0 and len([x for x in positions_with_clipshape if 'NULL' in x]) == 0:
                                #all basic checks passed.  add local position, 'rank', and label
                                print 'indeed we do. adding to list' 
                                local_start=int(-50)
                                for item in positions_with_clipshape:
                                    fmt_site=[local_start,str(mock_rank)+'_'+potential_newsite[0],mock_rank,potential_newsite[0],item[0],central_strand,item[1],item[2],item[3]]
                                    newsites.append(fmt_site)
                                    local_start+=1
                                #mock_rank+=1
                                #only one site allowed per transcript, so let's break out
                                return newsites
            return('No sites passed checks (length,clip/shape data present,clip=0)')
        else:
            return('No Sufficient Hits')
    else:
        return('Low RPKM')

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

fasta_dict=SeqIO.index(sys.argv[1], "fasta", alphabet=IUPAC.unambiguous_dna)

sys.stderr.write('fasta import complete: ' + str(datetime.now()) + os.linesep)

#import the seqs of interest - output of grab50.py
grabbed_positions=[]
with open(sys.argv[3]) as ingrab:
    for line in ingrab:
        grabbed_positions.append(line.strip().split('\t'))

sys.stderr.write('grabbed positions imported: ' + str(datetime.now()) + os.linesep)

#and import the clipmap output
clipmap_positions=[]
with open(sys.argv[2]) as inclip:
    for line in inclip:
        clipmap_positions.append(line.strip().split('\t'))
#trying this out as a set for speed, since there should be no duplicates

sys.stderr.write('clipmap imported: ' + str(datetime.now()) + os.linesep)

#import clipmap ranges, rather than recomputing them
clipranges_in=[]
with open(sys.argv[7]) as inrange:
    for line in inrange:
        clipranges_in.append(line.strip().split('\t'))

sys.stderr.write('clipranges imported: ' + str(datetime.now()) + os.linesep)


rpkmdict={}
with open(sys.argv[4]) as infile:
    for line in infile:
        inline=line.strip().split('\t')
        if '#' not in inline[0]:
            rpkmdict[inline[0]] = float(inline[4])


###########################
motif_len=int(sys.argv[5])
samp_count=max([int(x[2]) for x in grabbed_positions[1:]])

sys.stderr.write('starting the action: ' + str(datetime.now()) + os.linesep)

ranges=getRange(grabbed_positions)
sys.stderr.write('ranges gotten: ' + str(datetime.now()) + os.linesep)

ranges_withseq=getSeq(ranges,fasta_dict)
sys.stderr.write('sequences gotten: ' + str(datetime.now()) + os.linesep)

clip_motif=getMotif(ranges_withseq,motif_len)
sys.stderr.write('motif computed: ' + str(datetime.now()) + os.linesep)

clipmap_ranges=importClipRanges(clipranges_in)
sys.stderr.write('clipmap ranges imported: ' + str(datetime.now()) + os.linesep)

final_out = "%s.seqcomp_complete" % sys.argv[6]
f=open(final_out,'w')
f.write('LocalPos\tLabel\tRank\tTranscript\tChr\tStrand\tOrigPos\tShape\tClip\n')

true_negatives=0
mock_rank=1
clip_pssm=clip_motif.pssm
random.shuffle(clipmap_ranges)
for potential in clipmap_ranges:
    sys.stdout.write('investigating potential sites in: %s' % (potential[0])+os.linesep)
    site_with_data=findNewSitesFast(potential,clip_motif,clip_pssm,fasta_dict,clipmap_positions,rpkmdict,mock_rank)
    if type(site_with_data) is list:
        sys.stdout.write('got one! %s' % (potential[0])+os.linesep)
        sys.stdout.write('%s' % (site_with_data)+os.linesep)
        true_negatives+=1
        outputlist=flattenList(site_with_data)
        print outputlist
        print 'we now have this many sites: %s' % (true_negatives)
        f.write(outputlist+'\n')
        mock_rank+=1
    else:
        sys.stdout.write('no luck this time, because %s: %s' % (site_with_data,potential[0])+os.linesep)
    if true_negatives >= samp_count:
        break

sys.stderr.write('sites found: ' + str(datetime.now()) + os.linesep)

#outlist=flattenList(true_negatives)
#sys.stderr.write('done! writing final output: ' + str(datetime.now()) + os.linesep)

#f.write(outlist)
f.close()

sys.stderr.write('all\'s well that ends well: ' + str(datetime.now()) + os.linesep)
