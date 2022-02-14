# Copyright © 2021, Mengtian Li, Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China
# 
# ##python ConnectGraph-LTS-multiMean.py -d detail.txt -l RegLib.txt -t 40 -m 6 -M 1000
import numpy as np
import networkx as nx
import sys
import os
import copy
import multiprocessing
from multiprocessing import Process, Pool, Manager
import argparse
import random
import math
import FetchLink
from FetchLink import *
import Regression
from Regression import *


descript="This script is designed for scaffolding based on a robust regression model in RegScaf.\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-d', required = True, help='Detail file')
parser.add_argument('-l', required = True, help='Library file')
parser.add_argument('-e', required = True, help='Evidence file recording cutting scaffolds into contigs')
parser.add_argument('-t', type = int, default=1, help='number of threads')
parser.add_argument('-M', type = int,default=200, help='Max Error allowed in the final connected contig graph')
parser.add_argument('-z', type = int,default=2, help='the minimum size of the cluster considered in regression, should below the parameter -m but not too far')
parser.add_argument('-msd', type = int, default=6, help='minimum size of the cluster in which the library variance is replaced with the cluster variance')

def FindBreakPoint(scaf_name,select_row,Coff_mat,current_ctg):
	select_col = set()
	G = nx.Graph()
	for i in current_ctg:
		G.add_node(contig_list[i])
	pair = {}
	for row in select_row[1:]:
		line = list(Coff_mat[row, :])
		select_col.add(line.index(1))
		select_col.add(line.index(-1))
		n1 = current_ctg[line.index(1)]
		n2 = current_ctg[line.index(-1)]
		G.add_edge(contig_list[n1], contig_list[n2])
	allSubG = list(nx.connected_component_subgraphs(G))
	scaffold_broken = True
	broken_scaffolds = []
	for subG in allSubG:
		if contig_list[current_ctg[0]] in subG.nodes() and len(subG.nodes())==len(current_ctg):
			scaffold_broken = False
			broken_scaffolds = [list(set(G.nodes()))]
			break
		else:
			subG_nodes = list(set(subG.nodes()))
			broken_scaffolds.append(subG_nodes)
	fout_Broken = open("Broken_dir/Scaffold_"+str(scaf_name)+".txt",'w+')
	if scaffold_broken:
		print("The scaffold %s is broken? %s, and is broken into %d scaffolds!"%(scaf_name,scaffold_broken,len(broken_scaffolds)),file=fout_Broken)
	else:
		print("The scaffold %s is broken? %s!"%(scaf_name,scaffold_broken),file=fout_Broken)
	for l in broken_scaffolds:
		print(l,file=fout_Broken)
	fout_Broken.close()
	return scaffold_broken

def AssessScaf(scaf_name,scaf):
	global contig_list
	contig_in_subG = scaf['ctg_list']
	size = len(contig_in_subG)
	CtgDirection = np.array(np.ones(size))
	current_ctg = []
	for c in scaf['ctg_list']:  # current_ctg records the indexes of contig_in_subG in contig_list
		c_id = contig_list.index(c)
		current_ctg.append(c_id)
	Coff_mat = np.zeros((size*30,size))
	Weight_mat = np.zeros((size*30,size*30))
	r = 1	
	y_sep=[]
	y_sep.append(0)	#set the initial_coff==0
	Weight_mat[0,0]=100
	Coff_mat[0,0]=1
	for n1 in range(size):
		for n2 in range(size):
			if n1==n2:
				continue
			c1 = contig_in_subG[n1]
			c2 = contig_in_subG[n2]
			key = str(c1)+"-"+str(c2)
			for t in range(len(library_list)):
				elist=[]
				edge_dict_c = multiLib_EdgeDict[t]
				edge_dict_ic = multiLib_EdgeDict_ic[t]
				if CtgDirection[n1]==CtgDirection[n2]:
					if key in edge_dict_c:
						elist = edge_dict_c[key]
				else :
					if key in edge_dict_ic:
						elist = edge_dict_ic[key]
				if len(elist)<2: #few link in edge_dict
					continue
				#multi_mean = [(elist,np.median(elist))]
				msd = LibSd_List[t]
				if LibType_list[t]=='3GS':
					msd = 100
				if np.sqrt(np.var(elist))<msd:
					multi_mean = [elist]#[(elist,np.median(elist))]
				else:
					multi_mean = ClusterPeak(elist,msd)
					#multi_mean = ClusterDiff(elist,LibSd_List[t])
					#print("Cluster:",len(elist),len(multi_mean))
				for mi in multi_mean:
					if len(mi)<args.z:
						continue
					if len(mi)>=args.msd:
						msd = max(np.var(mi),1) #np.sqrt(np.var(mi))+1#
					else:
						msd = LibSd_List[t]**2
					Weight_mat[r,r] = 100*len(mi)/msd#np.sqrt(len(mi))/msd#'''(LibSd_List[t])np.sqrt(
					y_value = np.median(mi)
					if CtgDirection[n1]==-1:
						y_value=-y_value
					print("c1",c1,"c2",c2,len(multi_mean),len(mi),"%.2f"%np.sqrt(msd),"%.4f"%Weight_mat[r,r],y_value)
					y_sep.append(y_value)
					Coff_mat[r,n1] = -1
					Coff_mat[r,n2] = 1
					r +=1
	Coff_mat = Coff_mat[0:r,]
	Weight_mat = Weight_mat[0:r,0:r]
	y=np.array(y_sep)
	scaf_solution = scaf['position']
	Residual = np.dot(Coff_mat,np.array(scaf_solution))-y
	select_row = np.arange(r)[(Residual<args.M)]
	FindBreakPoint(scaf_name,select_row,Coff_mat,current_ctg)


def GetContigLen(ContigCount):
	DetailFile = open(args.d,'r')
	Contiglen_list = np.zeros(ContigCount)
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n-1]=int(detail[2])
	DetailFile.close()
	return Contiglen_list
	
args = parser.parse_args()
m = Manager()
sys.setrecursionlimit(10000)
libraryfile = open (args.l,'r')
maxerror = int(args.M)
library_name = []
library_list=[]
InsertSize_list=[]
LibType_list=[]
LibSd_List = []
for line in libraryfile:
	library_line = line.strip().split(' ')
	library_name.append(library_line[0])
	LibType_list.append(library_line[1])
	library_list.append(library_line[2])
	Libsd = round(float(library_line[4]),2)
	LibIns = int(library_line[3]) 
	LibIns_c = LibIns+int(Libsd**2/(1+LibIns))
	print("LibIns:",LibIns,LibIns_c)
	InsertSize_list.append(LibIns)
	LibSd_List.append(Libsd)
	
result = os.popen('tail -1 '+args.d)
res = result.read().split()
ContigCount = int(res[0][6:])
Contiglen_list = GetContigLen(ContigCount)

## Get link information from mapping result

contig_list = range(1,ContigCount+1)
direct_shortdict={}
direct_longdict={}
multiLib_EdgeDict = []
multiLib_EdgeDict_ic = []
link_dict = {}
for t in range(len(library_list)):
	tabFile = library_list[t]
	InsertSize = InsertSize_list[t]
	LibType = LibType_list[t]
	if LibType=='PE' or LibType=='MP':
		print("open new"+tabFile)
		tabfile=open('./new'+tabFile,'r')
	elif LibType=='3GS':
		print("open "+tabFile)
		tabfile = open('../'+tabFile,'r')
	edge_dict_c = {}	#save those consistent orientation diff
	edge_dict_ic = {}	#save those inconsistent orientation diff
	for pair in tabfile:
		if LibType=='PE':
			(c1,c2,diff,isconsist)=CalPEDiff(pair,InsertSize)
		elif LibType=='MP':
			(c1,c2,diff,isconsist)=CalMPDiff(pair,InsertSize)
		elif LibType=='3GS':
			(c1,c2,diff,isconsist)=GetLongLink(pair)
		else:
			break
		key = str(c1)+"-"+str(c2)
		key_1 = str(c2)+"-"+str(c1)
		if isconsist==1:
			if key in edge_dict_c:
				edge_dict_c[key].append(diff)
			else:
				edge_dict_c[key]=[diff]
		else:
			if key in edge_dict_ic:
				edge_dict_ic[key].append(diff)
			else:
				edge_dict_ic[key]=[diff]
		if key in link_dict:
			link_dict[key] +=1
			link_dict[key_1] +=1
		else:
			link_dict[key] =1
			link_dict[key_1] =1
		if LibType=='PE':
			if key in direct_shortdict:
				direct_shortdict[key].append(isconsist)
				direct_shortdict[key_1].append(isconsist)
			else:				
				direct_shortdict[key]=[isconsist]
				direct_shortdict[key_1]=[isconsist]
		else: #direct_long contains both MP Lib and 3GS Lib
			if key in direct_longdict:
				direct_longdict[key].append(isconsist)
				direct_longdict[key_1].append(isconsist)
			else:
				direct_longdict[key]=[isconsist]
				direct_longdict[key_1]=[isconsist]
	multiLib_EdgeDict.append(edge_dict_c)
	multiLib_EdgeDict_ic.append(edge_dict_ic)
	edgeFile = open(str(library_name[t])+"edge_dict_c",'w')
	edgeFile.write(str(edge_dict_c))
	edgeFile.close()
	edgeFile = open(str(library_name[t])+"edge_dict_ic",'w')
	edgeFile.write(str(edge_dict_ic))
	edgeFile.close()
	tabfile.close()

linkFile = open("link_dict",'w')
linkFile.write(str(link_dict))
linkFile.close()

directFile = open("direct_shortdict",'w')
directFile.write(str(direct_shortdict))
directFile.close()

directFile_2 = open("direct_longdict",'w')
directFile_2.write(str(direct_longdict))
directFile_2.close()

if (os.path.exists("Broken_dir")):
	setdir_com="rm -r Broken_dir"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir Broken_dir"
print(setdir_com)
os.system(setdir_com)

#Get scaffold information

##parallel computation prepare
m=Manager()
scaffold_dic=m.dict()
scaffold_name=0
last = 0
with open(args.e,'r') as fin:
    for line in fin:
        record = line.strip().split()
        if len(record) < 6 and record[0][0]!='~':
            if scaffold_name!=0:
                scaffold_dic[scaffold_name]={'ctg_list':ctg_list,"orientation":orientation,"length":length,"position":position}
                #print(scaffold_name,scaffold_dic[scaffold_name])
            scaffold_name=record[0]
            ctg_list = []
            orientation = []
            length = []
            position = []
        elif len(record) == 6:
            ctg_list.append(int(record[0]))
            orientation.append(record[1])
            length.append(int(record[2]))
            if len(position)==0:
                position.append(0)
            else:
                position.append(position[-1]+int(record[2])+int(record[4]))
    scaffold_dic[scaffold_name]={'ctg_list':ctg_list,"orientation":orientation,"length":length,"position":position}        


ScafFile = open("Scaf_dict",'w')
ScafFile.write(str(scaffold_dic))
ScafFile.close()
'''
for scaf_name,scaf_item in scaffold_dic.items():
	AssessScaf(scaf_name,scaf_item)

'''
pool=Pool(args.t)
broken_scaffold_dict = m.dict()
for scaf_name,scaf_item in scaffold_dic.items():
	pool.apply_async(func=AssessScaf,args=(scaf_name,scaf_item))
pool.close()
pool.join()
pool.terminate()
