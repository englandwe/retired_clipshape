#!/usr/bin/env python

#LocalPos       Label   Rank    Transcript	Chr     Strand  OrigPos Shape   Clip
#-50    1_ENSMUST00000172812    1	ENSMUST00000172812	19	-	5797379 0.000   0
#-49    1_ENSMUST00000172812    1	ENSMUST00000172812	19	-	5797380 0.192   0

#script to grab and combine matching transcripts from invitro and invivo data

import sys
import numpy

def matchTx(ds1,ds2):
    combolist=[]
    for entry1 in ds1:
        for entry2 in ds2:
            if [entry1[0]] + entry1[3:7] == [entry2[0]] + entry2[3:7]:
                combolist.append(entry1 + entry2[7:])
    return combolist


def getData(combolist):
    datalist=[]
    tx_list=list(set([x[1] for x in combolist]))
    for tx in tx_list:
        matching_sites=[x for x in combolist if x[1] == tx]
        if 'NULL' not in [x[7] for x in matching_sites] and 'NULL' not in [x[9] for x in matching_sites]:
            ds1_shapevals=[float(x[7]) for x in matching_sites if int(x[0]) in range(-10,11,1)]
            ds2_shapevals=[float(x[9]) for x in matching_sites if int(x[0]) in range(-10,11,1)]
            ds1_shapemean=numpy.mean(ds1_shapevals)
            ds2_shapemean=numpy.mean(ds2_shapevals)
            shapediff=list(numpy.subtract(ds2_shapevals,ds1_shapevals))
            shapediff_mean=numpy.mean(shapediff)
            #dsboth_shapemean=numpy.mean(ds1_shapevals+ds2_shapevals)
            centralclip=[x[8] for x in matching_sites if x[0] == '0'][0]
            centralpos=[x[6] for x in matching_sites if x[0] == '0'][0]
            datalist.append([tx,ds1_shapemean,ds2_shapemean,shapediff_mean,centralclip,centralpos])
    return datalist

def fmtForSVM(datalist):
    #0 1:val 2:val etc
    fmtlist=[]
    for line in datalist:
        newline=('1\t1:%s\t2:%s\t3:%s\t4:%s\t5:%s' % (line[1],line[2],line[3],line[4],line[5]))
        fmtlist.append(newline)
    final='\n'.join(fmtlist)
    return final



######
dataset1=[]
#with open('file1') as infile:
with open(sys.argv[1]) as infile:
    next(infile)
    for line in infile:
        dataset1.append(line.strip().split('\t'))

dataset2=[]
#with open('file2') as infile:
with open(sys.argv[2]) as infile:
    next(infile)
    for line in infile:
        dataset2.append(line.strip().split('\t'))

combined_data=fmtForSVM(getData(matchTx(dataset1,dataset2)))
print combined_data

#f=open('combined_data_out','w')
#f.write(combined_data+'\n')
#f.close()

