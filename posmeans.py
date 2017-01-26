import sys
import numpy

biglist=[]
with open(sys.argv[1]) as infile:
    for line in infile:
        biglist.append(line.strip().split('\t'))

#LocalPos	Label	Rank	Transcript	Chr	OrigPos	Shape	Clip
#-50	1_ENSMUST00000172812	1	ENSMUST00000172812	19	5797379	0.000	0

poslist=list(set([x[0] for x in biglist]))

outlist=[]
for pos in poslist:
    shapes=[float(x[6]) for x in biglist if x[0] == pos and x[6] not in ['NULL','Shape']]
    clips=[int(x[7]) for x in biglist if x[0] == pos and x[7] != 'Clip']
    shapemean=numpy.mean(shapes)
    clipmean=numpy.mean(clips)
    outlist.append([pos,shapemean,'shape'])
    outlist.append([pos,clipmean,'clip'])

list2=[]
for item in outlist:
    list2.append('\t'.join([str(x) for x in item]))
final='\n'.join(list2)
print final
