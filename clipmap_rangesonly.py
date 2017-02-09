#!/usr/bin/env python

#USAGE: shapedata clipdata gtf

#for transcript in shape_output
#find all its exons in gft file, grab their chr and positions (and order)
#retrieve the iCLIP data for each position
#retrieve the shape data for each position
#output something like this: geneid local_pos global_pos chr shapescore clipstops


import sys
import os
from itertools import groupby
from operator import itemgetter
#would be nice to use gffutils, but there's an issue with sqlite3 on the cluster.  could probably just make my own.

class GtfRec(object):
    def __init__(self,reclist):
         self.seqname=reclist[0]
         self.source=reclist[1]
         self.feature=reclist[2]
         self.start=int(reclist[3])
         self.end=int(reclist[4])
         self.score=reclist[5]
         self.strand=reclist[6]
         self.frame=reclist[7]
         self.attdict={}
         for self.item in reclist[8].strip(';').split('; '):
             self.splitline=self.item.replace('\"','').split(' ')
             self.attdict[self.splitline[0]]=self.splitline[1]


#a function to reassemble all the exons and stitch the proper positions together
#should return the txid and a list of positions and a chromosome#
def stitchTrans(ensembl_id,gtf_list):
    pos_list=[]
    chrid=''
    tmplist=[]
    for gtf_rec in gtf_list:
        #grab each exon's start & stop, plus exon number so they can be put in order
        if gtf_rec.feature == 'exon' and gtf_rec.attdict['transcript_id'] == ensembl_id:
            tmplist.append([gtf_rec.attdict['exon_number'],gtf_rec.start,gtf_rec.end+1])
            if len(chrid) == 0:
                chrid=gtf_rec.seqname
                strand=gtf_rec.strand
    #sort the list by exon number
    sorted_tmplist=sorted(tmplist, key=lambda x: int(x[0]))
    for item in sorted_tmplist:
        pos_list+=range(item[1],item[2])
    return [ensembl_id,chrid,pos_list,strand]

def outputRanges(stitched_trans):
    subrange_list=[]
    outlist=[]
    for k, g in groupby(enumerate(stitched_trans[2]), lambda (i, x): i-x):
        subrange_list.append(map(itemgetter(1), g))
    for subrange in subrange_list:
        sub_start=min(subrange)
        sub_stop=max(subrange)
        outlist.append([stitched_trans[0],stitched_trans[1],stitched_trans[3],sub_start,sub_stop])
    return outlist


#takes the stitched transcript and fetches iclip data.  returns input with clip data appended.

def clipMap(stitched_trans,clip_list):
    clipped_list={}
    #check that that chromosome is in the clip data
    clipchr=toClip(stitched_trans[1],chrnames)
    if clipchr == 'noclip':
        sys.stderr.write('WARNING: %s is mapped to %s, which has no CLIP data available. Skipping...' % (stitched_trans[0],stitched_trans[1]) + os.linesep)
        return 'FAIL'
    else:
        for pos in stitched_trans[2]:
            pos_key=clipchr+'_'+str(pos)
            if pos_key in clip_list:
                    clipped_list[pos]=clip_list[pos_key]
        return stitched_trans + [clipped_list]
        #if len(clipped_list) > 0:
        #    return stitched_trans + [clipped_list]
        #else:
            #sys.stderr.write('INFO: %s is mapped to a region (%s,%s-%s) which contains no CLIP sites. Skipping...' % (stitched_trans[0],stitched_trans[1],stitched_trans[2][0],stitched_trans[2][-1]) + os.linesep)
            #return 'FAIL'

#grabs shape data for a stitched transcript.  takes output of clipMap
def clipMerge(clipped_trans,shape_list):
    for shape in shape_list:
        if shape[0] == clipped_trans[0]:
            shapevals=shape[3:]
            if len(shapevals) != len(clipped_trans[2]):
                sys.stderr.write('ERROR:  length mismatch in %s. Pos: %s Shape: %s.' % (shape[0],len(clipped_trans[2]),len(shapevals)) + os.linesep)
                return 'FAIL'
            else:
                return clipped_trans + [shapevals]
            break

#final merge
def mergeAll(shaped_trans):
    finallist=[]
    clipped_pos=shaped_trans[4].keys()
    for i in range(len(shaped_trans[2])):
        tmplist=[shaped_trans[0], shaped_trans[1], shaped_trans[3], shaped_trans[2][i], shaped_trans[5][i]]
        if tmplist[3] in clipped_pos:
            tmplist.append(shaped_trans[4][tmplist[3]])
        else:
            tmplist.append('0')
        finallist.append(tmplist)
    return finallist

#translate shape/gtf chr names to clip; catch mappings to chrs with no clip data
def toClip(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'

def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final


###############################################
#Dictionary of shape/clip chr names
chrnames={
    '1' : 'chr1',
    '2' : 'chr2',
    '3' : 'chr3',
    '4' : 'chr4',
    '5' : 'chr5',
    '6' : 'chr6',
    '7' : 'chr7',
    '8' : 'chr8',
    '9' : 'chr9',
    '10' : 'chr10',
    '11' : 'chr11',
    '12' : 'chr12',
    '13' : 'chr13',
    '14' : 'chr14',
    '15' : 'chr15',
    '16' : 'chr16',
    '17' : 'chr17',
    '18' : 'chr18',
    '19' : 'chr19',
    'MT' : 'chrM',
    'X' : 'chrX',
    'Y' : 'chrY'
}


###############################################
#inputs

gtf_rec_list=[]
with open(sys.argv[3]) as infile:
#with open('testgtf.gtf') as infile:
    for line in infile:
        if not line.startswith('#'):
            gtf_rec_list.append(GtfRec(line.strip().split('\t')))

clip_data={}
with open(sys.argv[2]) as infile:
#with open('test.bedGraph') as infile:
    for line in infile:
        #clip_data.append(line.strip().split('\t'))
        cliptmp=line.strip().split('\t')
        newkey=cliptmp[0]+'_'+cliptmp[1]
        clip_data[newkey]=cliptmp[3]

shape_data=[]
with open(sys.argv[1]) as infile:
#with open('invitro-ischape.out.test') as infile:
    for line in infile:
        shape_data.append(line.strip().split('\t'))

############################################


ranges_outname = "%s.clipranges" % sys.argv[4]
f1=open(ranges_outname,'w')

#clip_outname = "%s.clipmap" % sys.argv[4]
#f2=open(clip_outname,'w')


#ENSMUST00000002289
id_list=[x[0] for x in shape_data]
for id in id_list:
    trans_stitched=stitchTrans(id,gtf_rec_list)
    trans_ranges=outputRanges(trans_stitched)
#    trans_clipped=clipMap(trans_stitched,clip_data)
#    if trans_clipped != 'FAIL':
#        trans_shaped=clipMerge(trans_clipped,shape_data)
#        if trans_shaped != 'FAIL':
#            trans_merged=mergeAll(trans_shaped)

    trans_out=flattenList(trans_ranges)
    f1.write(trans_out+'\n')

#            clip_out=flattenList(trans_merged)
#            f2.write(clip_out)

f1.close()
#f2.close()
