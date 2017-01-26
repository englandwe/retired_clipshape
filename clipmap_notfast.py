#!/usr/bin/env python
import sys
import os

#ARGS: Shape, clip, gtf, outfile

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


#Stitch a transcript together.
#in order to search for transcripts associated with a clip site, will need to do all of them and make a list
#takes the transcritp ID and a list of GTF objects
def stitchTrans(ensembl_id,gtf_list):
    pos_list=[]
    chrid=''
    for gtf_rec in gtf_list:
        if gtf_rec.feature == 'exon' and gtf_rec.attdict['transcript_id'] == ensembl_id:
            pos_list+=range(gtf_rec.start,gtf_rec.end+1)
            if len(chrid) == 0:
                chrid=gtf_rec.seqname
    return [ensembl_id,chrid,pos_list]


#translate clip chr names to shape/gtf; catch mappings to chrs with no clip data
def toShape(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'


#take all the clip sites, make list of associated transcripts
#now we have two options
#1: take the list of associated transcripts and map on clip/shape as before
#2: map each clip site to the transcript as we go.  modify the tx-list in-place, then cycle through all the modified sites until everyone is mapped.  then go grab the shape data as before
#2  would seem to be faster.

def clipMap(stitched_list,clip_list):
    #copy list
    stitched_list_2=list(stitched_list)
    for clippos in clip_list:
        #translate clip-chrname to shape-chrname
        shapechr=toShape(clippos[0],chrnames)
        #find the transcript(s) on the right chr including your position
        #matched_tx=[x.append( for x in stitched_list if x[1] == clippos[0] and clippos[1] in x[2]]
        if shapechr == 'noclip':
            sys.stderr.write('WARNING: The sequence called %s in CLIP data cannot be mapped to SHAPE data.  Check your sequence names. Skipping...' % (clippos[0]) + os.linesep)
        else:
            for tx in stitched_list_2:
                if tx[1] == shapechr and int(clippos[1]) in tx[2]:
                    #if there isn't a clip dictionary in the entry yet, make one. keys=position, values=stops
                    if len(tx) == 3:
                        tx.append({})
                    tx[3][int(clippos[1])]=int(clippos[3])
    return stitched_list_2

#grabs shape data for stitched transcripts if they have clip data.  takes output of clipMap.
def clipMerge(clipped_list,shape_list):
    merged_list=[]
    for tx in clipped_list:
        #only do stuff if you have clip data
        if len(tx) == 4:
            for shape in shape_list:
                if shape[0] == tx[0]:
                    shapevals=shape[3:]
                    if len(shapevals) != len(tx[2]):
                        sys.stderr.write('ERROR:  length mismatch in %s. Pos: %s Shape: %s.' % (shape[0],len(tx[2]),len(shapevals)) + os.linesep)
                    merged_list.append(tx+[shapevals])
        elif len(tx) == 3:
            sys.stderr.write('INFO: Transcript %s is mapped to a region which contains no CLIP sites. Skipping...' % (tx[0]) + os.linesep)
        else:
            sys.stderr.write('ERROR: Length of entry in stitched transcript list is incorrect. Debug: %s.' % (tx) + os.linesep)
    return merged_list

#final merge; takes output of clipMerge
def mergeAll(shaped_trans):
    finallist=[]
    for tx in shaped_trans:
        clipped_pos=tx[3].keys()
        for i in range(len(tx[2])):
            tmplist=[tx[0], tx[1], tx[2][i], tx[4][i]]
            if tmplist[2] in clipped_pos:
                tmplist.append(tx[3][tmplist[2]])
            else:
                tmplist.append('0')
            finallist.append(tmplist)
    return finallist


def flattenList(inlist):
    list2=[]
    for item in inlist:
        list2.append('\t'.join ([str(x) for x in item]))
    final='\n'.join(list2)
    return final

###############################################
#Dictionary of clip > shape chr names
chrnames={
 'chr1': '1',
 'chr10': '10',
 'chr11': '11',
 'chr12': '12',
 'chr13': '13',
 'chr14': '14',
 'chr15': '15',
 'chr16': '16',
 'chr17': '17',
 'chr18': '18',
 'chr19': '19',
 'chr2': '2',
 'chr3': '3',
 'chr4': '4',
 'chr5': '5',
 'chr6': '6',
 'chr7': '7',
 'chr8': '8',
 'chr9': '9',
 'chrM': 'MT',
 'chrX': 'X',
 'chrY': 'Y'
}

###############################################
#inputs

gtf_rec_list=[]
#with open(sys.argv[3]) as infile:
with open('testgtf.gtf') as infile:
    for line in infile:
        if not line.startswith('#'):
            gtf_rec_list.append(GtfRec(line.strip().split('\t')))

clip_data=[]
#with open(sys.argv[2]) as infile:
with open('test.bedGraph') as infile:
    for line in infile:
        clip_data.append(line.strip().split('\t'))
        #cliptmp=line.strip().split('\t')
        #newkey=cliptmp[0]+'_'+cliptmp[1]
        #clip_data[newkey]=cliptmp[3]

shape_data=[]
#with open(sys.argv[1]) as infile:
with open('invitro-ischape.out.test') as infile:
    for line in infile:
        shape_data.append(line.strip().split('\t'))

############################################

#grab all the transcript ids and stitch them together
id_list=[x[0] for x in shape_data]
trans_stitched=[]
for id in id_list:
    trans_stitched.append(stitchTrans(id,gtf_rec_list))

#the actual mapping bit
trans_clipped=clipMap(trans_stitched,clip_data)
trans_shaped=clipMerge(trans_clipped,shape_data)
trans_merged=mergeAll(trans_shaped)

outlist=flattenList(trans_merged)

f=open(sys.argv[4],'w')
f.write('Transcript\tChr\tPos\tShape\tClip\n')
f.write(outlist)
f.close()

