#python3 MergeFinal.py -d detail.txt -sd scriptsdir -t 40 -mi maxlibinsert -ms maxlibsd -ml minctglen 
import multiprocessing
from multiprocessing import Process, Pool, Manager, Lock
import numpy as np
import sys
import os
import math
import argparse
import re

descript="Split Finals\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-d', required=True, help='Detail file')
parser.add_argument('-sd', type = str, required = True,help='scripts directory')
parser.add_argument('-t', type=int, default=1, help='number of threads')
parser.add_argument('-m', type=int,default=6, help='minimum link as a edge in the contig graph')
parser.add_argument('-mi', type = int, default = 8000,help = 'maxmum insert size of using libraries')
parser.add_argument('-ms', type = int,default = 1000, help ='maxmum insert size variation of using libraries')
parser.add_argument('-ml', type = int ,default = 1000, help = 'minumum contig length that we fetch link from, below which link count is fetched from its longer preccessor')
parser.add_argument('--minscore', type = int, default= 40, help = 'minumum alignment score to Merge')
args = parser.parse_args()

result = os.popen('tail -1 '+str(args.d))
res = result.read().split()
MaxLength = int(res[0][6:])
print ("Original contig number:",MaxLength)

def ReadFinal(File,flag):
	RegOut = open("Final_dir/"+File,'r')
	Reglines = RegOut.readlines()
	contig_dict=[]
	for line in Reglines[1:]:
		record = line.strip().split(' ')
		#ctg_list.append(int(record[0]))
		contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6])),'newctg':0})
	#ctg_num = len(contig_dict)	
	return contig_dict
	
def ProSplitFunc(Final,final_name):
	count = re.match(r"(\d+)",final_name[6:]).group(0)
	writescaf_command="python3 "+args.sd+"SplitFinal-revise.py Final_dir/"+final_name+" "+str(args.m)+" "+str(args.mi)+" "+str(args.ms)
	print(writescaf_command)
	os.system(writescaf_command)
	print("Scaf.size",len(Final))

if (os.path.exists("SplitFinal")):
	setdir_com="rm -r SplitFinal"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir SplitFinal"
print(setdir_com)
os.system(setdir_com)

list_com = "ls -rt Final_dir/ > allFinal.txt"
print(list_com)
os.system(list_com)

##Handle all regression final results in final_dir,merge results in parallel
allFinal = open("allFinal.txt",'r')
pool = multiprocessing.Pool(args.t)
scaf_num=0
for line in allFinal:
	final_name = line.strip()
	Final=ReadFinal(final_name,1)
	Final.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	pool.apply_async(ProSplitFunc,args=(Final,final_name,))
	scaf_num += 1
pool.close()
pool.join()
pool.terminate()


