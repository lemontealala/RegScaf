'''
>scaffold1|size4271505|tigs56
f_tig861|size181867|links60|gaps-435p-907P89
f_tig862|size63958|links4|gaps-1448p-1448P-1448
f_tig863|size90449|links26|gaps-339p-924P192
f_tig864|size71389|links23|gaps-447p-716P166
f_tig865|size3815|links3|gaps-835p-993P-835
f_tig866|size104415|links13|gaps-342p-945P-57
f_tig867|size495821|links6|gaps-896p-896P-494
f_tig868|size170460|links6|gaps-1657p-1657P2288
f_tig869|size27721|links8|gaps-1959p-1959P-1532
'''
#python3 MakeEvidence.py evidence.txt
import sys
import os
import re
import copy 
import numpy as np

def ReCalGap(ctg_dict):
	ctg_num = len(ctg_dict)
	i = 0
	contig_dict = copy.deepcopy(ctg_dict)
	while i < ctg_num-1:
		gap_new = contig_dict[i+1]['start']-contig_dict[i]['end']
		contig_dict[i]['gap'] = int(gap_new)
		i +=1
	contig_dict[i]['gap'] = 0
	return contig_dict

def Readfinal(File,log_ctg):
	RegOut = open("final_dir/"+File,'r')
	Reglines = RegOut.readlines()
	contig_dict=[]
	for line in Reglines[1:]:
		record = line.strip().split(' ')
		contig = int(record[0])
		if contig in log_ctg:
			contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6]))})
	RegOut.close()
	contig_dict= ReCalGap(contig_dict)
	scaf_len = contig_dict[-1]['end']-contig_dict[0]['start']
	return (contig_dict,scaf_len)
	
def ReadLog(File):
	log_list=[]#[[] for i in range()]
	log_ctg = []
	LogOut = open(File,'r')
	Loglines =LogOut.readlines()
	i=1
	while i <len(Loglines):
		line =Loglines[i]
		if line[0]=='>':
			log_list.append(log_ctg)
			log_ctg=[]
			i+=1
		elif line[0]=='c':	##contig
			log_ctg.append(int(line[6:]))
			i+=2
		elif line[0]=='N':	
			if line[1]=='*':	#N*
				i+=1
			else:	##NewCtgn1+n2
				ctg_dict=line[6:].split('+')
				for ctg in ctg_dict:
					log_ctg.append(int(ctg))
				i+=2
		else:
			i+=1
	log_list.append(log_ctg)
	LogOut.close()
	return log_list
	
def GetContigLen(ContigCount):
	DetailFile = open("detail.txt",'r')
	Contiglen_list = np.zeros(ContigCount)
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n-1]=int(detail[2])
	DetailFile.close()
	return Contiglen_list
		
result = os.popen('tail -1 detail.txt')
res = result.read().split()
ContigCount = int(res[0][6:])
Contiglen_list = GetContigLen(ContigCount)

EvidenceFile=open(sys.argv[1],'w+')
list_com = "ls final_dir/|sort -k 1.7n > allfinal.txt"
print(list_com)
os.system(list_com)
allfinal = open("allfinal.txt",'r')
scaf_num = 0
logrecord = open("logRecord",'w+')
for line in allfinal:
	finalName=line.strip()
	pattern=re.compile(r'(?<=final_)\d+')
	final_num=pattern.findall(finalName)[0]
	log_name = "log_dir/scaffold"+final_num+".log"	
	log_list=ReadLog(log_name)
	print(log_name,len(log_list),scaf_num,file=logrecord)
	for log_ctg in log_list:
		scaf_num+=1
		(contig_dict,scaf_len)=Readfinal(finalName,log_ctg)
		print(">scaffold"+str(scaf_num)+"|size"+str(scaf_len)+"|tigs"+str(len(contig_dict)),file=EvidenceFile)
		for contig in contig_dict:
			std = int(max(contig['std'],300))
			gap=int(contig['gap'])
			gap_p=gap-std
			gap_P=gap+std
			print(str.lower(contig['direct'])+"_tig"+str(contig['contig'])+"|size"+str(contig['len'])+"|links10|gaps"+str(gap)+"p"+str(gap_p)+"P"+str(gap_P),file=EvidenceFile)
		print("\n",end="",file=EvidenceFile)
logrecord.close()
allfinal.close()

if (os.path.exists("smallCtg.txt")):
	smallctg = open("smallCtg.txt",'r')
	for line in smallctg:
		ctg_id = int(line)
		scaf_num +=1
		scaf_len = Contiglen_list[ctg_id-1]
		print(">scaffold"+str(scaf_num)+"|size"+str(scaf_len)+"|tigs1",file=EvidenceFile)
		print("f_tig"+str(ctg_id)+"|size"+str(scaf_len)+"|links10|gaps0p0P0\n",file=EvidenceFile)
	smallctg.close()
	
EvidenceFile.close()

