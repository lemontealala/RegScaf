#python ReCalLibInsertSize.py detail.txt RegLib.txt

import sys
import numpy as np
from collections import Counter
import os
	
def GetContigLen(ContigCount):
	DetailFile = open(sys.argv[1],'r')
	Contiglen_list = np.zeros(ContigCount)
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n-1]=int(detail[2])
	DetailFile.close()
	return Contiglen_list
	
libraryfile = open(sys.argv[2],'r')	
LibName_list = []
LibType_list=[]
tab_list=[]
InsertSize_list=[]
sd_list = []
for line in libraryfile:
	library_line = line.strip().split(' ')
	LibName_list.append(library_line[0])
	LibType_list.append(library_line[1])
	tab_list.append(library_line[2])
	InsertSize_list.append(int(float(library_line[3])))
	sd_list.append(float(library_line[4]))

result = os.popen('tail -1 '+sys.argv[1])
res = result.read().split()
ContigCount = int(res[0][6:])
Contiglen_list = GetContigLen(ContigCount)

for t in range(len(LibName_list)):
	tabFile = open("../"+tab_list[t],'r')
	if LibType_list[t]!= 'PE' and LibType_list[t]!= 'MP':
		print(LibName_list[t],LibType_list[t],tab_list[t],0,0,file = open("RegLib_n.txt",'a'))
		continue
	InsertSize = InsertSize_list[t]
	sd = sd_list[t]
	In_list = []
	for line in tabFile:
		pairline=line.strip().split("\t")
		c1=int(pairline[0][6:])
		c2=int(pairline[3][6:])
		if c1 == c2:
			CtgLength = Contiglen_list[c1-1]
			if CtgLength<3*InsertSize:	continue
			p1 = int(pairline[1])
			m1 = int(pairline[2])
			p2 = int(pairline[4])
			m2 = int(pairline[5])
			p_Max = max(p1,m1,p2,m2)
			p_Min = min(p1,m1,p2,m2)
			x = p_Max-p_Min
			if x>InsertSize+3*sd or x< InsertSize-3*sd:
				continue
			In_list.append(x)
			#if len(In_list)>=300000:
			#	break
	tabFile.close()
	y=np.array(In_list)
	if len(In_list)>1:
		print("Lib:",LibName_list[t],"Median:",np.median(y),"Mean:",np.mean(y),"std",np.std(y,ddof=1),"ListLength",len(y))
		print(LibName_list[t],LibType_list[t],tab_list[t],int(np.median(y)),np.std(y,ddof=1),file = open("RegLib_n.txt",'a'))
	else:
		print("Lib useless!:",LibName_list[t])
update_com = 'mv RegLib_n.txt '+sys.argv[2]
os.system(update_com)
