#python3 MergeFinal.py -d detail.txt -o outpre -sd scriptsdir -t 48 -mi maxlibinsert -ms maxlibsd -ml minctglen --minscore 40 
import multiprocessing
from multiprocessing import Process, Pool, Manager, Lock
import numpy as np
import sys
import os
import math
import argparse
import re

descript="RegScaf for initial contigs\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-d', required=True, help='Detail file')
parser.add_argument('-o', required=True, help='Output prefix')
parser.add_argument('-sd', type = str, required = True,help='scripts directory')
parser.add_argument('-t', type=int, default=1, help='number of threads')
parser.add_argument('-m', type=int,default=6, help='minimum link as a edge in the contig graph')
parser.add_argument('-mi', type = int, default = 8000,help = 'maxmum insert size of using libraries')
parser.add_argument('-ms', type = int,default = 1000, help ='maxmum insert size variation of using libraries')
parser.add_argument('-ml', type = int ,default = 1000, help = 'minumum contig length that we fetch link from, below which link count is fetched from its longer preccessor')
parser.add_argument('--minscore', type = int, default= 40, help = 'minumum alignment score to Merge')
parser.add_argument('--hangingout', type = int, default=20, help= 'Allow longer hanging out when merge two overlapped contigs, set it larger when initial contigs contain many errors')
args = parser.parse_args()

result = os.popen('tail -1 '+str(args.d))
res = result.read().split()
MaxLength = int(res[0][6:])
print ("Original contig number:",MaxLength)
Scaf_final = str(args.o)+"_scaf.fa"
Log_final = str(args.o)+"_scaf.log"

def revcomp(seq):
	seq=seq[::-1]
	seq_rc=''
	for i in range(len(seq)):
		if seq[i]=='A' or seq[i]=='a':
			seq_rc+='T'
		elif seq[i]=='T' or seq[i]=='t':
			seq_rc+='A'
		elif seq[i]=='C' or seq[i]=='c':
			seq_rc+='G'
		elif seq[i]=='G' or seq[i]=='g':
			seq_rc+='C'
		else:
			seq_rc+=seq[i]
	return(seq_rc)

def ReadFinal(File,flag):
	if flag==1:	
		RegOut = open("final_dir/"+File,'r')
	else:
		RegOut = open("Final_dir_supp/"+File,'r')
	Reglines = RegOut.readlines()
	contig_dict=[]
	for line in Reglines[1:]:
		record = line.strip().split(' ')
		#ctg_list.append(int(record[0]))
		contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6])),'newctg':0})
	#ctg_num = len(contig_dict)	
	return contig_dict

def ReadContig(contig_num,contig_dir):
	filename='../CONTIG/contig'+str(contig_num)+'.fasta'
	file = open(filename,'r')
	filelines = file.read().splitlines()
	ctg_seq = ''
	for line in filelines[1:]:
		ctg_seq = ctg_seq+line
	if contig_dir=='R':
		print ("Reverse……")
		ctg_seq = revcomp(ctg_seq)
	return ctg_seq
	
def WriteSingleCtg(ctg_id,scaf_num,ScafFile,ScafLog):
	print(">scaffold"+str(scaf_num),end ='',file=ScafFile)
	print(">scaffold"+str(scaf_num)+"\t",file=ScafLog)
	scaffold=ReadContig(ctg_id,'F')
	base_num = 0
	for base in scaffold:
		if(base_num % 100 == 0): print('\n',end ='',file = ScafFile)
		base_num +=1
		print(base,end ='',file= ScafFile)
	print("contig"+str(ctg_id),file=ScafLog)
	print("\n",end ='',file=ScafFile)

def ProMergeFunc(Final,final_name):
	count = re.match(r"(\d+)",final_name[6:]).group(0)
	scaf_name = "scaf_dir/scaffold"+str(count)+".fa"
	log_name = "log_dir/scaffold"+str(count)+".log"
	writescaf_command="python3 "+args.sd+"/MergeAndWriteNewScaffold-revise.py "+str(count)+" final_dir/"+final_name+" "+scaf_name+" "+log_name+" "+str(args.m)+" "+str(args.mi)+" "+str(args.ms)+" "+str(args.minscore)+" "+str(args.hangingout)
	print(writescaf_command)
	os.system(writescaf_command)
	print("Scaf.size",len(Final))

if (os.path.exists(Scaf_final)):
	setdir_com="rm "+str(Scaf_final)
	print(setdir_com)
	os.system(setdir_com)
setdir_com="touch "+str(Scaf_final)
print(setdir_com)
os.system(setdir_com)

if (os.path.exists(Log_final)):
	setdir_com="rm "+str(Log_final)
	print(setdir_com)
	os.system(setdir_com)
setdir_com="touch "+str(Log_final)
print(setdir_com)
os.system(setdir_com)

if (os.path.exists("scaf_dir")):
	setdir_com="rm -r scaf_dir"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir scaf_dir"
print(setdir_com)
os.system(setdir_com)

if (os.path.exists("log_dir")):
	setdir_com="rm -r log_dir"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir log_dir"
print(setdir_com)
os.system(setdir_com)

list_com = "ls -rt final_dir/ > allfinal.txt"
print(list_com)
os.system(list_com)

##Handle all regression final results in final_dir,merge results in parallel
allFinal = open("allfinal.txt",'r')
pool = multiprocessing.Pool(args.t)
scaf_num=0
for line in allFinal:
	final_name = line.strip()
	Final=ReadFinal(final_name,1)
	Final.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	pool.apply_async(ProMergeFunc,args=(Final,final_name,))
	scaf_num += 1
pool.close()
pool.join()
pool.terminate()

##Write Scaffolds and detail logs
writeScaf_com = "cat scaf_dir/* > Allscaffold.fa"
print(writeScaf_com)
os.system(writeScaf_com)

writeScaf_com = "cat log_dir/* >  Allscaffold.log"#+str(Log_final)
print(writeScaf_com)
os.system(writeScaf_com)
if os.path.exists("smallCtg.txt"):
	ScafFile=open("Allscaffold.fa",'a')
	ScafLog=open("Allscaffold.log",'a')
	smallCtg=open("smallCtg.txt",'r')
	for line in smallCtg:
		ctg_id=str(line.strip())
		WriteSingleCtg(ctg_id,scaf_num,ScafFile,ScafLog)
		scaf_num+=1	
	smallCtg.close()
	ScafFile.close()
	ScafLog.close()
print("total Scaf_num",scaf_num-1)
##Sort Scaffolds
sortScaf_com = "python "+str(args.sd)+"/SortScaffold.py Allscaffold.fa "+str(Scaf_final)
print(sortScaf_com)
os.system(sortScaf_com)
sortScaf_com = "python "+str(args.sd)+"/SortScaffold.py Allscaffold.log "+str(Log_final)
print(sortScaf_com)
os.system(sortScaf_com)
rm_com = "rm Allscaffold.fa"
print(rm_com)
os.system(rm_com)
rm_com = "rm Allscaffold.log"
print(rm_com)
os.system(rm_com)
