# Copyright © 2021, Mengtian Li, Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import numpy as np
import networkx as nx
import sys
import math
import os
import multiprocessing
from multiprocessing import Process, Pool
import json
import datetime
from numpy import *
import argparse
import copy
import Direction_DFS_Correct
from Direction_DFS_Correct import *
from MST_Direction import *
import random
import Regression
from Regression import *

descript="The scaffolding procedure for repeat scaffolding on super-contigs\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-d', required = True, help='Detail file')
parser.add_argument('-l', required = True, help='Library file')
parser.add_argument('-t', type = int, default=1, help='number of threads')
parser.add_argument('-m', type = int,default=6, help='minimum link as a edge in the pre-scaffold graph')
parser.add_argument('-mg', type = int ,default=500, help='max graph size restricted for better effect')
parser.add_argument('-M', type = int,default=500, help='Max Error allowed in the final connected pre-scaffold graph')

def MedianStd(list):
	list=np.array(list)
	size=len(list)
	if size<2: return 0
	#mean=np.mean(list1)
	median=np.median(list)
	m_list=median*np.ones(size)
	Std=sum((list-m_list)**2)/size
	return np.sqrt(Std/2)
	
def Select_Col(select_row,Coff_mat,Weight_mat,current_scaf,truc):
	select_col = set()
	G=nx.Graph()
	for i in current_scaf:
		G.add_node(scaf_list[i])
	for row in select_row[1:]:
		line=list(Coff_mat[row,:])
		select_col.add(line.index(1))
		select_col.add(line.index(-1))
		n1=current_scaf[line.index(1)]
		n2=current_scaf[line.index(-1)]
		#if Weight_mat[row,row]*2>truc:
		G.add_edge(scaf_list[n1],scaf_list[n2])
	allSubG = nx.connected_component_subgraphs(G)
	IsConnect = 0
	for subG in allSubG:
		print("Here subG.size",len(subG.nodes()))
		if scaf_list[current_scaf[0]] in subG.nodes() and len(subG.nodes())==len(current_scaf):		
			IsConnect = 1
			print("Here IsConnect!")
		else:
			subG_nodes = list(set(subG.nodes()))
			if len(subG_nodes)<2:
				WriteFailfinalFile(subG_nodes[0])
			else:
				ProRegScafFunc(subG_nodes)
	return(IsConnect,list(select_col))
	
def DDFS(visit,n1,DFlag,weight_mat,CtgDirection,contig):	
	weight = copy.deepcopy(weight_mat[n1,])
	max_nei = np.argmax(weight)
	link_weight = weight[max_nei]
	weight[max_nei] = 0
	while link_weight>0:
		if visit[str(max_nei)]==0 and DFlag[n1][max_nei]!=0 :
			CtgDirection[max_nei] = DFlag[n1][max_nei]*CtgDirection[n1]
			visit[str(max_nei)]=1
			#print("from",contig[n1],'visit',contig[max_nei],CtgDirection[max_nei],link_weight)
			DDFS(visit,max_nei,DFlag,weight_mat,CtgDirection,contig)
		max_nei = np.argmax(weight)
		#print('from ',n1,'to the next',max_nei,weight_mat[n1,max_nei],weight[max_nei])
		link_weight = weight[max_nei]
		weight[max_nei] = 0

def ScafDirect(scaf_direct_dict,scaf_in_subG,truc):
	PreScafCount = len(scaf_in_subG)
	SDFlag = np.zeros((PreScafCount,PreScafCount))
	J_mat = np.zeros((PreScafCount,PreScafCount))
	weight_mat = np.zeros((PreScafCount,PreScafCount))
	for n1 in range(PreScafCount):
		for n2 in range(n1,PreScafCount):
			c1 = scaf_in_subG[n1]
			c2 = scaf_in_subG[n2]
			key = str(c1)+'-'+str(c2)
			d_list = []
			direct_sum = 0
			if key in scaf_direct_dict:		
				d_list = np.array(scaf_direct_dict[key])
				direct_sum = sum(d_list)
				posi_count = (sum(abs(d_list))+direct_sum)/2
				nega_count = (sum(abs(d_list))-direct_sum)/2
				w = max(posi_count,nega_count)
				weight_mat[n1][n2] += w
				weight_mat[n2][n1] += w	
				J_mat[n1][n2] = direct_sum
				J_mat[n2][n1] = direct_sum
				if w>min(20,truc):# and min(posi_count,nega_count)<truc:
					SDFlag[n1][n2]= np.sign(direct_sum)
					SDFlag[n2][n1]= np.sign(direct_sum)
	#np.savetxt("Direct-WeightMat",weight_mat[0:30,0:30],fmt = '%d')
	##Initialize:
	PreScafDirection = [0 for i in range(PreScafCount)]
	PreScafDirection[0] = 1
	visit = {}
	for i in range(PreScafCount):
		visit[str(i)]=0
	visit[str(0)]=1
	DDFS(visit,0,SDFlag,weight_mat,PreScafDirection,scaf_in_subG)
	for i in range(PreScafCount):	
		if PreScafDirection[i]==0:
			print("scaf_in_subGs not initialized:",scaf_in_subG[i])
			PreScafDirection[i] = 1
	##optimization
	TC = CalTotalContradict(PreScafDirection,J_mat)
	min_ML = -10
	TC_old = -100
	while min_ML<0 and TC>TC_old:
		TC_old = TC
		ML = DeltaTC(PreScafDirection,J_mat)
		min_ML = min(ML)
		max_1 = np.argmax(-ML)
		print("max_1:",scaf_in_subG[max_1],ML[max_1])
		#if max_1 ==0:  break
		if min_ML==0: 
			print("Perfect orientation!")
			break
		PreScafDirection[max_1]=-PreScafDirection[max_1]
		TC = CalTotalContradict(PreScafDirection,J_mat)
		print("TC_new:", TC, TC_old)
		if TC<=TC_old:	#change back!
			PreScafDirection[max_1] = -PreScafDirection[max_1]
	return PreScafDirection

def ProRegScafFunc(scaf_in_subG):
	global scaf_direct_dict
	global Scaflen_list
	global scaf_edge_dict_c
	global scaf_edge_dict_ic
	global scaf_list
	#scaf_in_subG = list(set(subG.nodes()))
	size = len(scaf_in_subG)
	print("scaf_subG:",scaf_subG)
	if size<2:
		WriteFailfinalFile(scaf_in_subG[0])
		return
	if size>1000:
		print("ERROR! Please take more strict criterion in graph constructing! (Incrase the parameter: -m)")
		return
	Coff_mat = np.zeros((size*200,size))
	Cov_mat = np.zeros((size*200,size*200))
	Weight_mat = np.zeros((size*200,size*200))
	PreScafDirection = ScafDirect(scaf_direct_dict,scaf_in_subG,truc)
	#PreScafDirection = ScafMSTDirect(scaf_direct_dict,scaf_in_subG)
	print("Orientation assignment ... ")
	r = 1	
	y_sep=[]
	y_sep.append(0)
	Weight_mat[0,0]=100
	Coff_mat[0,0]=1
	Cov_mat[0,0]=0
	for n1 in range(size):
		for n2 in range(size):
			if n1==n2:
				continue
			s1 = scaf_in_subG[n1]
			s2 = scaf_in_subG[n2]
			key = str(s1)+"-"+str(s2)
			s_elist = []
			if PreScafDirection[n1]==PreScafDirection[n2]:
				if key in scaf_edge_dict_c:
					s_elist=scaf_edge_dict_c[key]
			else:
				if key in scaf_edge_dict_ic:
					s_elist=scaf_edge_dict_ic[key]
			scaf_diff = []
			for diff in s_elist:
				scaf_diff +=list(diff)
			if len(scaf_diff)<truc/2:
				continue
			if np.sqrt(np.var(scaf_diff))<max_sd:
				multi_scafdiff = [scaf_diff]
			else:
				multi_scafdiff = ClusterPeak(scaf_diff,max_sd)
			for si in multi_scafdiff:
				if len(si)<3:
					continue
				msd = max_sd#np.var(si)+1 #np.sqrt(np.var(si))+1 #
				Weight_mat[r,r] = len(si)/msd
				if PreScafDirection[n1]==1:
					y_sep.append(np.median(si))
					print("s1",s1,"s2",s2,len(multi_scafdiff),np.median(si),len(si),msd,Weight_mat[r,r])
				else:
					y_sep.append(-np.median(si))
					print("s1",s1,"s2",s2,len(multi_scafdiff),-np.median(si),len(si),msd,Weight_mat[r,r])
				Cov_mat[r,r] = msd		
				Coff_mat[r,n1] = -1
				Coff_mat[r,n2] = 1
				r +=1
			'''
			for diff in s_elist :
				if len(diff)<3:
					continue
				Weight_mat[r,r] = len(diff)
				if PreScafDirection[n1]==1:
					y_sep.append(np.median(diff))
				else:
					y_sep.append(-np.median(diff))
				Cov_mat[r,r] = np.var(diff)		#here should be corrected later
				Coff_mat[r,n1] = -1
				Coff_mat[r,n2] = 1
				r +=1
				'''
	Coff_mat = Coff_mat[0:r,]
	Weight_mat = Weight_mat[0:r,0:r]
	Cov_mat = Cov_mat[0:r,0:r]
	y=array(y_sep)
	
	#Iteranally Regression:	
	sub_error=[10000]
	iter=1
	select_row= np.arange(r)
	select_row_last = np.arange(r)
	select_col= np.arange(size)
	break_flag=0
	iter_len = max(1,math.floor(len(y)*0.01))
	current_scaf=[]
	for c in scaf_in_subG:
		c_id = scaf_list.index(c)
		current_scaf.append(c_id)
	while max(sub_error)>maxerror and iter<300:
		print("Before next regression:Coff_mat.shape:",Coff_mat.shape,"current_scaf",len(current_scaf))
		(IsConnect,select_col) = Select_Col(select_row,Coff_mat,Weight_mat,current_scaf,truc)
		if IsConnect==0:
			return
		print("After Select Column:",len(select_row),len(select_col))
		if len(select_col)<2:
			break_flag=1
			break
		Solution_mat=RegProcSol(Coff_mat[select_row,:][:,select_col],Weight_mat[select_row,:][:,select_row],y[select_row])
		y_predict = np.dot(Coff_mat,Solution_mat)
		y_diff = y_predict - y
		abs_error = abs(y_diff)
		sub_error = abs_error[select_row]
		print("The ",iter,"iter:max_error=",max(abs_error),"submaxerror",max(sub_error),"submean_error",mean(sub_error),"90% quantile:",np.percentile(abs_error,90),"Solution_mat.shape",Solution_mat.shape)
		weighted_error = abs_error#/np.sqrt(np.diagonal(Weight_mat))
		d = sorted(weighted_error)[len(y)-iter_len*iter]
		'''
		d = np.percentile(weighted_error,100-iter)
		#d = FindDemarcation(weighted_error,maxerror)
		d = min(d,max(weighted_error)-0.0001)'''
		print("After Find Demarcation:", d)
		select_row_last = copy.deepcopy(select_row)
		select_row = np.arange(r)[(abs_error<=d)]
		#select_row = sorted(random.sample(list(select_row),int(len(select_row)*0.95)))
		iter +=1
	print("total iter:",iter-1)
	if max(sub_error)>maxerror or break_flag==1:	#All iteration failed in shrinkage of err
		for scaf in current_scaf:
			WriteFailfinalFile(scaf)
	else :
		Solution_Cov = SolCovFinal(Coff_mat[select_row_last,:][:,select_col],Weight_mat[select_row_last,:][:,select_row_last],Cov_mat[select_row_last,:][:,select_row_last])
		WritefinalFile(Solution_mat,Solution_Cov,current_scaf,scaf_in_subG,PreScafDirection)
	
def WriteFailfinalFile(s_num):
	global Final_list
	global final_num
	name = Final_list[s_num]
	com = "cp SplitFinal/"+str(name)+" ./final_dir/final_"+str(final_num)+".txt"
	print(com)
	res = os.system(com)
	if not res:
		print("OS Failure!:",res)
		supply_com = open("./supply_com.sh",'a+')
		print(com,file=supply_com)
		supply_com.close()
	final_num+=1
	
def ReversePreScaf(scaf):
	scaf_t = copy.deepcopy(scaf)
	scaf_size = len(scaf)
	scaf_len = scaf[scaf_size-1]['end']
	for i in range(scaf_size):
		ctg = scaf[i]
		if ctg['direct']=='F':
			scaf_t[i]['direct']='R'
		elif ctg['direct']=='R':
			scaf_t[i]['direct']='F'
		st = ctg['start']
		scaf_t[i]['start'] = scaf_len -ctg['end']
		scaf_t[i]['end'] = scaf_len-st
	return scaf_t

def ShiftPreScaf(scaf,direct,position):
	scaf_t = copy.deepcopy(scaf)
	scaf_size = len(scaf)
	scaf_len = scaf[scaf_size-1]['end']
	if direct == 1:
		for i in range(scaf_size):
			scaf_t[i]['start'] += position 
			scaf_t[i]['end'] += position
		return scaf_t
	else:
		for i in range(scaf_size):
			if scaf_t[i]['direct']=='F':
				scaf_t[i]['direct']='R'
			elif scaf_t[i]['direct']=='R':
				scaf_t[i]['direct']='F'
			st_t = scaf_t[i]['start']
			scaf_t[i]['start'] = position-scaf_t[i]['end']	#position +scaf_len-scaf_t[i]['end']
			scaf_t[i]['end'] = position-st_t	#position+scaf_len-st_t
		return scaf_t

def WritefinalFile(Solution_mat,Solution_cov,current_scaf,scaf_in_subG,PreScafDirection):	
	global final_num
	Final = Solution_mat
	Final_cov = Solution_cov
	print("Final.shape:",Final.shape,"Final_num",final_num,"Solution_mat:",Final)
	FinalScaf = []
	for i in range(len(Final)):
		k_in_current = current_scaf[i]
		scaf_id = scaf_list[k_in_current]
		k = scaf_in_subG.index(scaf_id)
		scaf_direct = PreScafDirection[k]
		print(Final_list[k_in_current],k_in_current,scaf_direct)
		scaf = ScafFinal[scaf_id]
		diff = Final[i]-Final[0]
		'''if scaf_direct == -1:
			scaf = ReversePreScaf(scaf)
		for ctg in scaf:
			ctg['start'] += diff 
			ctg['end'] += diff
			FinalScaf.append(ctg)'''
		shifted_scaf = ShiftPreScaf(scaf,scaf_direct,diff)
		FinalScaf +=shifted_scaf
	FinalScaf.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	size = len(FinalScaf)
	scafmaxlen = FinalScaf[size-1]['end']-FinalScaf[0]['start']
	fout_Final = open("final_dir/final_"+str(final_num)+".txt",'w+')
	print(("longest position is %.1f" % scafmaxlen),file = fout_Final)
	for k in range(len(FinalScaf)):
		dict=FinalScaf[k]
		if k== len(FinalScaf)-1:
			gap=0
		else:	
			next=FinalScaf[k+1]
			gap=next['start']-dict['end']
		print(dict['contig'], '%.1f' % dict['start'], '%.1f' % dict['end'],	dict['direct'],'%.1f' % dict['std'],'%.1f' % gap,'%d' % dict['len'], file=fout_Final)
	fout_Final.close()
	final_num +=1
	
def GetContigLen(ContigCount):
	DetailFile = open(args.d,'r')
	Contiglen_list = np.zeros(ContigCount)
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n-1]=int(detail[2])
	DetailFile.close()
	return Contiglen_list

def ReadFinal(File,flag):
	if flag==1:	
		RegOut = open("SplitFinal/"+File,'r')
	else:
		RegOut = open("Final_dir_supp/"+File,'r')
	Reglines = RegOut.readlines()
	first_record = Reglines[1].strip().split(' ')
	shift = -float(first_record[1])
	contig_dict=[]
	ctg_list = []
	for line in Reglines[1:]:
		record = line.strip().split(' ')
		ctg_list.append(int(record[0]))
		#contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6])),'newctg':0})
		contig_dict.append({'contig':int(record[0]),'start':float(record[1])+shift,'end':float(record[2])+shift,'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'len':int(float(record[6])),'newctg':0})
	#ctg_num = len(contig_dict)	
	return (contig_dict,ctg_list)
	
args = parser.parse_args()
libraryfile = open (args.l,'r')
truc = int(args.m)
maxerror = int(args.M)
library_list=[]
InsertSize_list=[]
LibSd_List = []
library_name = []
LibType_list=[]
for line in libraryfile:
	library_line = line.strip().split(' ')
	library_name.append(library_line[0])
	LibType_list.append(library_line[1])
	library_list.append(library_line[2])
	InsertSize_list.append(int(library_line[3]))
	LibSd_List.append(round(float(library_line[4]),2))
max_sd = max(LibSd_List)

if (os.path.exists("supply_com.sh")):
	rm_com = "mv supply_com.sh supply_com.sh_old"
	print(rm_com)
	os.system(rm_com)

## get connected contig graph
result = os.popen('tail -1 '+args.d)
res = result.read().split()
ContigCount = int(res[0][6:])
Contiglen_list = GetContigLen(ContigCount)

directFile = open("direct_shortdict",'r')
b = directFile.read()
direct_shortdict = eval(b)
directFile.close()

directFile_2 = open("direct_longdict",'r')
c = directFile_2.read()
direct_longdict = eval(c)
directFile_2.close()

list_com = "ls -rt SplitFinal/ > allFinal.txt"
print(list_com)
os.system(list_com)

scaf_num = 0
final_num = 0
inconsistent_pair = 0
result = os.popen('wc -l allFinal.txt')
res = result.read().split()
total = int(res[0])
scaf_list = range(total)
allFinal = open("allFinal.txt",'r')
ScafFinal=[[] for x in range(total)]
scaf_G=nx.Graph()
max_graph_size = int(args.mg)
scaf_edge_dict_c={}
scaf_edge_dict_ic={}
scaf_direct_dict={}
scaf_link_dict ={}
ctg_dict = {}
Final_list=[]
ScafCtg = [[] for x in range(total)]
for line in allFinal:
	(Final,ctg_list) = ReadFinal(line.strip(),1)
	Final_list.append(line.strip())
	scaf_G.add_node(scaf_num)
	ScafFinal[scaf_num] = Final
	ScafCtg[scaf_num] = ctg_list
	for ctg in ctg_list:
		ctg_s = str(ctg)
		ctg_dict[ctg_s]=scaf_num
	scaf_num +=1

scaf_diffInf = open("scaf_diff",'w+')
for t in range(len(library_list)):
	edgeFile = open(str(library_name[t])+"edge_dict_c",'r')
	print("open",str(library_name[t])+"edge_dict_c")
	a = edgeFile.read()
	edge_dict_c = eval(a)
	edgeFile_ic = open(str(library_name[t])+"edge_dict_ic",'r')
	print("open",str(library_name[t])+"edge_dict_ic")
	b = edgeFile_ic.read()
	edge_dict_ic = eval(b)
	for key in set(edge_dict_c.keys()).union(set(edge_dict_ic.keys())):
		key_inf = key.split('-')
		c1=str(key_inf[0])
		c2=str(key_inf[1])
		if not (c1 in ctg_dict and c2 in ctg_dict):	continue
		s1= ctg_dict[c1]
		s2= ctg_dict[c2]
		if s1==s2 :	
			continue
		s_key = str(s1)+"-"+str(s2)
		s_key_1 = str(s2)+"-"+str(s1)
		c1_id = ScafCtg[s1].index(int(c1))
		c2_id = ScafCtg[s2].index(int(c2))
		st1 = ScafFinal[s1][c1_id]['start']
		st2 = ScafFinal[s2][c2_id]['start']
		en1 = ScafFinal[s1][c1_id]['end']
		en2 = ScafFinal[s2][c2_id]['end']
		d1 = ScafFinal[s1][c1_id]['direct']
		d2 = ScafFinal[s2][c2_id]['direct']
		s1_size = len(ScafFinal[s1])
		s2_size = len(ScafFinal[s2])
		s1_len = ScafFinal[s1][s1_size-1]['end']
		s2_len = ScafFinal[s2][s2_size-1]['end']
		elist_c =[]
		elist_ic = []
		if key in edge_dict_c:
			elist_c = edge_dict_c[key]
		if key in edge_dict_ic:
			elist_ic = edge_dict_ic[key]
		if len(elist_c)+len(elist_ic)<1:
			continue
		ctg_direct = len(elist_c)-len(elist_ic)
		link = max(len(elist_c),len(elist_ic))
		if ctg_direct>=0: #link between contig pair is consistent
			diff = np.array(elist_c)
			if d1=='F' and d2=='F' : #case 3
				scaf_diff = diff+st1-st2
				s_consistent = 1
			elif d1=='R' and d2=='R': #case 5
				scaf_diff = -diff+en1-en2
				s_consistent = 1
			elif d1=='F' and d2=='R': #case 4
				scaf_diff = diff+st1+en2	
				s_consistent = -1
			elif d1=='R' and d2=='F': #case 6
				scaf_diff = -diff+en1+st2	
				s_consistent = -1
		else:
			diff = np.array(elist_ic)
			if d1=='F' and d2=='F' : #case 2
				scaf_diff = diff+st1+st2
				s_consistent = -1
			elif d1=='R' and d2=='R': #case 8
				scaf_diff = -diff+en1+en2
				s_consistent = -1
			elif d1=='F' and d2=='R': #case 1
				scaf_diff = diff+st1-en2
				s_consistent = 1
			elif d1=='R' and d2=='F': #case 7
				scaf_diff = -diff+en1-st2
				s_consistent = 1
		if len(diff)>0:
			print("c1",c1,"c2",c2,"s1",s1,"s2",s2,"ctg_direct",ctg_direct,"s_consistent",s_consistent,d1,d2,"diff",diff,"st1",st1,"st2",st2,s1_len,s2_len,"scaf_diff",[x for x in scaf_diff],file = scaf_diffInf)
		if np.sqrt(np.var(scaf_diff))>2000:
			print("too large sd!",file = scaf_diffInf)
			#continue
		if s_key in scaf_direct_dict:
			scaf_direct_dict[s_key].append(s_consistent*link)
			scaf_direct_dict[s_key_1].append(s_consistent*link)
		else:
			scaf_direct_dict[s_key]=[s_consistent*link]
			scaf_direct_dict[s_key_1]=[s_consistent*link]
		if s_consistent==1:
			if s_key in scaf_edge_dict_c:
				scaf_edge_dict_c[s_key].append(scaf_diff)
			else:
				scaf_edge_dict_c[s_key] =[scaf_diff]
		else:
			if s_key in scaf_edge_dict_ic:
				scaf_edge_dict_ic[s_key].append(scaf_diff)
			else:
				scaf_edge_dict_ic[s_key] =[scaf_diff]
		if s_key in scaf_link_dict:
			scaf_link_dict[s_key]+= link
			scaf_link_dict[s_key_1]+= link
		else:
			scaf_link_dict[s_key] = len(diff)
			scaf_link_dict[s_key_1] = len(diff)
	print("Here",len(scaf_edge_dict_c),len(scaf_edge_dict_ic))
	edgeFile.close()
	edgeFile_ic.close()
scaf_diffInf.close()

directFile = open("scaf_direct_dict",'w')
directFile.write(str(scaf_direct_dict))
directFile.close()

directFile = open("scaf_edge_dict_c",'w')
directFile.write(str(scaf_edge_dict_c))
directFile.close()

directFile = open("scaf_edge_dict_ic",'w')
directFile.write(str(scaf_edge_dict_ic))
directFile.close()

scaf_EdgeNodes = set()
for s_key in scaf_link_dict:
	key_inf = s_key.split('-')
	s1=int(key_inf[0])
	s2=int(key_inf[1])
	d_list = np.array(scaf_direct_dict[s_key])
	direct_sum = sum(d_list)
	posi_count = (sum(abs(d_list))+direct_sum)/2
	nega_count = (sum(abs(d_list))-direct_sum)/2
	w = max(posi_count,nega_count)
	if w>truc :#and min(posi_count,nega_count)<truc:
		scaf_G.add_edge(s1,s2)
		scaf_EdgeNodes.add(s1)
		scaf_EdgeNodes.add(s2)
print("scaf_G edges:",len(scaf_G.edges),"scaf_EdgeNodes",len(scaf_EdgeNodes))		

allScafSubG = nx.connected_component_subgraphs(scaf_G)
s_record = open("scaf_subG",'w+')
in_count = 0
if (os.path.exists("final_dir")):
	setdir_com="rm -r final_dir"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir final_dir"
print(setdir_com)
os.system(setdir_com)
for scaf_subG in allScafSubG:
	if len(scaf_subG.nodes)>1:
		if len(scaf_subG.nodes)>max_graph_size:
			cut_n = int(len(scaf_subG.nodes)/max_graph_size)
			nodes = list(set(scaf_subG.nodes))
			for i in range(cut_n):
				scaf_subG_nodes = nodes[(i*max_graph_size):(i+1)*max_graph_size]
				in_count+=len(scaf_subG_nodes)
				print("subG:",len(scaf_subG_nodes),scaf_subG_nodes,file=s_record)
				ProRegScafFunc(scaf_subG_nodes)
			scaf_subG_nodes = nodes[max_graph_size*cut_n:]
		else:
			scaf_subG_nodes = list(set(scaf_subG.nodes()))
		in_count+=len(scaf_subG_nodes)
		print("subG:",len(scaf_subG_nodes),scaf_subG_nodes,file=s_record)
		ProRegScafFunc(scaf_subG_nodes)
	else:
		s_num = list(scaf_subG.nodes)[0]
		WriteFailfinalFile(s_num)
s_record.close()

##in case some command failed:
if (os.path.exists("supply_com.sh")):
	supply = "bash supply_com.sh"
	print(supply)
	os.system(supply)

print("nodes in scaf_subG size>1:",in_count)
