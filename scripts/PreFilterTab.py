#python PreFilterTab.py detail.txt library.txt 

import sys
import numpy as np
import os

ContigCount=700000
DetailFile = open(sys.argv[1],'r')
libraryfile = open(sys.argv[2],'r')
	
library_list=[]
InsertSize_list=[]	
ReadLength_list=[]
LibType_list=[]
SD_list=[]
for line in libraryfile:
	library_line = line.strip().split(' ')
	if library_line[1]=='3GS':
		continue
	LibType_list.append(library_line[1])
	library_list.append(library_line[2])
	InsertSize_list.append(int(library_line[3]))
	SD_list.append(float(library_line[4]))
	
#Get Contig Length from DetailFile:
Contiglen_list = np.zeros(ContigCount)
for line in DetailFile:
	detail = line.strip().split("\t")
	ctg_n = int(detail[0][6:])
	Contiglen_list[ctg_n]=int(detail[2])
DetailFile.close()

sort_comFile=open("sortTab.sh",'w')
for t in range(len(library_list)):
	tabFile=library_list[t]
	SortTabCommand="awk '$1!=$4' ../"+str(tabFile)+" > ./"+str(tabFile)+" &"#"sort -u ../"+str(tabFile)+" | awk '$1 !=$4' > ./"+str(tabFile)+" &"
	print(SortTabCommand,file=sort_comFile)
print("wait",file=sort_comFile)
sort_comFile.close()
os.system("bash sortTab.sh")

for t in range(len(library_list)):
	tabFile=library_list[t]
	tabfile=open('./'+tabFile,'r')
	outFile=open('./new'+tabFile,'w')
	InsertSize=InsertSize_list[t]
	LibType=LibType_list[t]
	sd=SD_list[t]
	print("InsertSize",InsertSize,"SD",sd)
	for pair in tabfile:
		pairline=pair.strip().split("\t")
		c1=int(pairline[0][6:])
		c2=int(pairline[3][6:])
		p1 = int(pairline[1])
		m1 = int(pairline[2])
		p2 = int(pairline[4])
		m2 = int(pairline[5])
		s1 = min(p1,m1)
		s2 = min(p2,m2)
		len_1 = Contiglen_list[c1]
		len_2 = Contiglen_list[c2]
		if LibType=='MP':
			if p1>m1 and p2<m2:
				if max(len_1-m1,m2)>(3*sd+InsertSize):#len_1-p1>(3*sd+InsertSize) or p2>(3*sd+InsertSize):
					continue
			elif p1<m1 and p2<m2:
				if max(m1,m2)>(3*sd+InsertSize):#p1>(3*sd+InsertSize) or p2>(3*sd+InsertSize):	
					continue			
			elif p1>m1 and p2>m2:
				if max(len_1-m1,len_2-m2)>(3*sd+InsertSize):#len_1-p1>(3*sd+InsertSize) or len_2-p2>(3*sd+InsertSize):	
					continue	
			elif p1<m1 and p2>m2:
				if max(m1,len_2-m2)>(3*sd+InsertSize):#p1>(3*sd+InsertSize) or len_2-p2>(3*sd+InsertSize):	
					continue	
		elif LibType=='PE':
			if p1<m1 and p2>m2:
				if max(len_1-p1,p2)>3*sd+InsertSize:	
					#print("len_1-p1+p2",len_1,len_2,p1,p2,len_1-p1+p2,"3*sd+InsertSize",3*sd+InsertSize)
					continue#len_1-m1>3*sd+InsertSize or m2>3*sd+InsertSize
			elif p1>m1 and p2>m2:
				if max(p1,p2)>3*sd+InsertSize:	
					#print("p1+p2",len_1,len_2,p1,p2,p1+p2,"3*sd+InsertSize",3*sd+InsertSize)
					continue#m1>3*sd+InsertSize or m2>3*sd+InsertSize
			elif p1<m1 and p2<m2:
				if max(len_1-p1,len_2-p2)>3*sd+InsertSize:
					#print("len_1-p1+len_2-p2",len_1,len_2,p1,p2,len_1-p1+len_2-p2,"3*sd+InsertSize",3*sd+InsertSize)
					continue#len_1-m1>3*sd+InsertSize or len_2-m2>3*sd+InsertSize
			elif p1>m1 and p2<m2:
				if max(p1,len_2-p2)>3*sd+InsertSize:	
					#print("p1+len_2-p2",len_1,len_2,p1,p2,p1+len_2-p2,"3*sd+InsertSize",3*sd+InsertSize)
					continue	#m1>3*sd+InsertSize or len_2-m2>3*sd+InsertSize
		print(pair,end ='',file=outFile)
	tabfile.close()
	outFile.close()
