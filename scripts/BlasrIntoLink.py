#python BlasrIntoLink.py BLASROUT 3GSTab
import sys
import numpy as np
		
def CalAdjDiff(List,DiffList,CTG_dict):
    i=0
    while i <len(List):
        i_cover=CTG_dict[List[i]['tName']]/Contiglen_dict[List[i]['tName']]
        if i_cover>cover_truc:##filter high coverage contigs
            print("High Coverage:",List[i]['qName'],"contig"+str(List[i]['tName']),Contiglen_dict[List[i]['tName']],"coverage:",i_cover)
            del List[i]
        else:
            i+=1
    i = 0
    while i<len(List)-1:
        j=i+1
        i_del = 0
        while j<len(List):
            if List[j]['qStart']<List[i]['qEnd'] : #mapping region overlap,delete records with tLength<1000
                if List[i]['tLength']<1000:
                    i_del=1
                if List[j]['tLength']<1000:
                    print("Mapping overlap:\n",List[i]['qName'],"contig"+str(List[i]['tName']),List[i]['qStart'],List[i]['qEnd'],'\n',List[j]['qName'],"contig"+str(List[j]['tName']),List[j]['qStart'],List[j]['qEnd'])
                    del List[j]
                else:
                    j+=1
            else:
                break
        if i_del==1:
            del List[i]
        else:
            i+=1
    if len(List)<2:
        return
    print("Final Record on Long Read:",List[0]['qName'],len(List))
    print([x['tName'] for x in List])
    i=0
    while i<len(List)-1:
        j=i+1
        while j<len(List): #and List[j]['qStart']<List[i]['qEnd']+10000:
            #if min(List[j]['tLength'],List[i]['tLength'])<1000 and j>i+1:
            #    j+=1
            #    continue
            if List[i]['tName']==List[j]['tName']:# or List[j]['qEnd']<List[i]['qEnd']+200:
                j+=1
                continue
            if List[i]['tStrand']==0 and List[j]['tStrand']==0:
                dir = 0
                x2 = List[j]['qStart']- List[j]['tStart']
                x1 = List[i]['qStart']- List[i]['tStart']
            elif List[i]['tStrand']==0 and List[j]['tStrand']==1:
                dir = 1
                x1 = List[i]['qStart']- List[i]['tStart']
                x2 = List[j]['qStart']+(List[j]['tLength']-List[j]['tStart'])
            elif List[i]['tStrand']==1 and List[j]['tStrand']==0:
                dir = 1
                x2 = List[j]['qStart']- List[j]['tStart']
                x1 = List[i]['qStart']+(List[i]['tLength']-List[i]['tStart'])
            else:
                dir = 0
                x1 = List[i]['qStart']+(List[i]['tLength']-List[i]['tStart'])
                x2 = List[j]['qStart']+(List[j]['tLength']-List[j]['tStart'])
            diff = x2-x1	
            if List[i]['tStrand']==1:
                diff = -diff
            print(List[i]['tName'],List[j]['tName'],dir,diff,-List[i]['score'],-List[j]['score'],file=DiffList)
            j+=1
        i +=1

DetailFile = open("detail.txt",'r')
Contiglen_dict = {}
for line in DetailFile:
    detail = line.strip().split("\t")
    ctg_n = int(detail[0][6:])
    Contiglen_dict[ctg_n]=int(detail[2])
DetailFile.close()

BLASR = open(sys.argv[1],'r')
DiffList = open(sys.argv[2],'w')
CTG_dict = {}
LR_dict={}
align_inf = ['qName','tName','qStrand','tStrand','score','percentSimilarity','tStart','tEnd','tLength','qStart','qEnd','qLength']
for line in BLASR:
    bline = line.strip().split(" ")
    LRead = bline[0]
    align_dict={'qName':bline[0],'tName':int(bline[1][6:]),'qStrand':int(bline[2]),'tStrand':int(bline[3]),'score':int(bline[4]),'percentSimilarity':float(bline[5]),'tStart':int(bline[6]),'tEnd':int(bline[7]),'tLength':int(bline[8]),'qStart':int(bline[9]),'qEnd':int(bline[10]),'qLength':int(bline[11])}
    tName = int(bline[1][6:])
    if tName in CTG_dict:
        CTG_dict[tName]+=int(bline[7])-int(bline[6])
    else:
        CTG_dict[tName]=int(bline[7])-int(bline[6])
    if LRead in LR_dict:
        LR_dict[LRead].append(align_dict)
    else:
        LR_dict[LRead]=[]
        LR_dict[LRead].append(align_dict)

#np.savetxt('LongReadCount',np.array(CTG_dict.items()),fmt='%d')
coverage = []
CovFile = open("LongReadCVG",'w')
for key in CTG_dict:
    c = CTG_dict[key]/Contiglen_dict[key]
    coverage.append(c)
    print("%d  %d  %d  %.2f" %(key,CTG_dict[key],Contiglen_dict[key],c),file=CovFile)
CovFile.close()
#np.savetxt('LongReadCVG',coverage,fmt='%.3f')
cover_truc = np.percentile(coverage,80)
print("90% quantail of cover number:",cover_truc)

for LR_key,align_list in LR_dict.items():
    if len(align_list)<2:
        continue
    align_list.sort(key=lambda obj:(obj.get('qStart')),reverse=False)
    CalAdjDiff(align_list,DiffList,CTG_dict)

BLASR.close()    
DiffList.close()