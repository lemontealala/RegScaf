# Copyright © 2021, Mengtian Li, Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China
# 
# ##python ConnectGraph-LTS-multiMean.py -d detail.txt -l RegLib.txt -t 40 -m 6 -M 1000
import numpy as np
import networkx as nx
import sys
import os
import copy
import Direction_DFS_Correct
from Direction_DFS_Correct import *
from MST_Direction import *
#from ClusterDiff import *
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
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
parser.add_argument('-t', type = int, default=1, help='number of threads')
parser.add_argument('-m', type = int,default=0, help='the minimum link number to add an edge in initial contig graph')
parser.add_argument('-M', type = int,default=200, help='Max Error allowed in the final connected contig graph')
parser.add_argument('-z', type = int,default=2, help='the minimum size of the cluster considered in regression, should below the parameter -m but not too far')
parser.add_argument('-s', type = int ,default=300, help='minimum length of the contig considered in the connected contig graph')
parser.add_argument('-mg', type = int ,default=3000, help='max graph size restricted for better effect')
#parser.add_argument('-c', action="store_true", help= 'Correct Insert-size according GapEst')
parser.add_argument('-msd', type = int, default=6, help='minimum size of the cluster in which the library variance is replaced with the cluster variance')

def WriteFinalFile(Solution_mat,Solution_cov,current_ctg,Contig_order,CtgDirection):
	global Contiglen_list
	global contig_list
	global scaffed_ctg
	Final = Solution_mat
	Final_cov = Solution_cov
	first_ctg = contig_list[current_ctg[0]]
	print("Final.shape:",Final.shape,"Final_cov.shape:",Final_cov.shape,"Final_name:",first_ctg)
	FinalContig = []
	scafmin = 0
	scafmax = Contiglen_list[Contig_order[0]-1]
	for i in range(len(Final)):
		k_in_current=current_ctg[i]
		ctg = contig_list[k_in_current]
		k = Contig_order.index(ctg)
		scaffed_ctg.append(ctg)
		if CtgDirection[k]==1:
			FinalContig.append({'contig':ctg,'start':Final[i],'end':Final[i]+Contiglen_list[ctg-1],'direct':'F','std':Final_cov[i,i],'len':Contiglen_list[ctg-1]})
			if scafmin>Final[i]:	scafmin=Final[i]
			if scafmax<Final[i]+Contiglen_list[ctg-1]: scafmax=Final[i]+Contiglen_list[ctg-1]
		elif CtgDirection[k]==-1:
			FinalContig.append({'contig':ctg,'start':Final[i]-Contiglen_list[ctg-1],'end':Final[i],'direct':'R','std':Final_cov[i,i],'len':Contiglen_list[ctg-1]})
			if scafmin>Final[i]:	scafmin=Final[i]-Contiglen_list[ctg-1]
			if scafmax<Final[i]+Contiglen_list[ctg-1]: scafmax=Final[i]
	Initial = copy.deepcopy(FinalContig)
	FinalContig.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	size = len(FinalContig)
	final_cov = np.zeros((size,size))
	i=0
	sort_order = []
	while i<size:
		index = Initial.index(FinalContig[i])
		sort_order.append(index)
		i+=1
	gap_var = []
	print("FinalContig.size:",size)
	i=0
	while i<size-1:
		j= sort_order[i]
		j_1 = sort_order[i+1]
		gap_var.append(Final_cov[j,j]+Final_cov[j_1,j_1]-Final_cov[j,j_1]-Final_cov[j_1,j])
		i+=1
	gap_var.append(0)
	scafmaxlen = scafmax-scafmin	
	fout_Final = open("Final_dir/Final_"+str(first_ctg)+".txt",'w+')
	print(("longest position is %.1f" % scafmaxlen),file = fout_Final)
	for k in range(size):
		dict=FinalContig[k]
		if k== size-1:
			gap = 0 
			std = 0
		else:	
			next = FinalContig[k+1]
			std = gap_var[k]	#gap_cov[k,k]
			gap=next['start']-dict['end']
		print(dict['contig'], '%.1f' % dict['start'], '%.1f' % dict['end'],	dict['direct'],'%.1f' % np.sqrt(std),'%.1f' % gap,'%d' % dict['len'], file=fout_Final)
	fout_Final.close()

def Select_Col(select_row,Coff_mat,Weight_mat,current_ctg,truc):
	select_col = set()
	G=nx.Graph()
	for i in current_ctg:
		G.add_node(contig_list[i])
	pair = {}
	for row in select_row[1:]:
		line=list(Coff_mat[row,:])
		select_col.add(line.index(1))
		select_col.add(line.index(-1))
		n1=current_ctg[line.index(1)]
		n2=current_ctg[line.index(-1)]
		G.add_edge(contig_list[n1],contig_list[n2])
		'''		a = min(n1,n2)
		b = max(n1,n2)
		pair_key = str(a)+"-"+str(b)
		if pair_key in pair:
			pair[pair_key] += Weight_mat[row,row]
		else:
			pair[pair_key] = Weight_mat[row,row]
	print("Here Pair:",pair)
	for key,item in pair.items():
		if int(item)>truc:
			[n1,n2]=key.split('-')
			print("Add Edge:",n1,n2)
			G.add_edge(contig_list[int(n1)],contig_list[int(n2)])'''
	allSubG = nx.connected_component_subgraphs(G)
	IsConnect = 0
	for subG in allSubG:
		print("Here subG.size",len(subG.nodes()))
		if contig_list[current_ctg[0]] in subG.nodes() and len(subG.nodes())==len(current_ctg):		
			IsConnect = 1
			print("Here IsConnect!")
		else:
			subG_nodes = list(set(subG.nodes()))
			ProRegFunc(subG_nodes)
	return(IsConnect,list(select_col))

def WriteSingleCtg_InSubProcess(ctg):
	global Contiglen_list
	#if os.path.exists("Final_dir/Final_"+str(ctg)+".txt"):
	#	print("Error!", ctg)
	fout = open("Final_dir/Final_"+str(ctg)+".txt",'w+')
	ctg_len = Contiglen_list[ctg-1]
	print(("longest position is %.1f" % ctg_len),file = fout)
	print(ctg, '%.1f' % 0, '%.1f' % ctg_len,'F', '%.1f' % 0,  '%.1f' % 0, '%d' % ctg_len, file=fout)
	fout.close()

def ProRegFunc(contig_in_subG):
	global Contiglen_list
	global direct_shortdict
	global direct_longdict
	global contig_list
	size = len(contig_in_subG)
	print("subG.size:",size)
	if size<2:
		WriteSingleCtg_InSubProcess(contig_in_subG[0])
		return
	Coff_mat = np.zeros((size*30,size))
	Cov_mat = np.zeros((size*30,size*30))
	Weight_mat = np.zeros((size*30,size*30))
	print("Start Orientation: contig_in_subG",contig_in_subG)
	#CtgDirection = DirectMST(direct_shortdict,direct_longdict,contig_in_subG,Contiglen_list)
	CtgDirection = DirectC(direct_shortdict,direct_longdict,contig_in_subG,truc)#,Contiglen_list
	print("Orientation assignment...")
	#for i in range(size):
	#	print("contig",contig_in_subG[i],CtgDirection[i])
	r = 1	
	y_sep=[]
	y_sep.append(0)	#set the initial_coff==0
	Weight_mat[0,0]=100
	Coff_mat[0,0]=1
	Cov_mat[0,0]=0
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
					'''if len(multi_mean)>1:
						y_value = np.mean(mi)#np.median(mi)
					else:'''
					y_value = np.median(mi)
					if CtgDirection[n1]==-1:
						y_value=-y_value
					print("c1",c1,"c2",c2,len(multi_mean),len(mi),"%.2f"%np.sqrt(msd),"%.4f"%Weight_mat[r,r],y_value)
					y_sep.append(y_value)
					Cov_mat[r,r] = msd
					Coff_mat[r,n1] = -1
					Coff_mat[r,n2] = 1
					r +=1
	Coff_mat = Coff_mat[0:r,]
	Weight_mat = Weight_mat[0:r,0:r]
	Cov_mat = Cov_mat[0:r,0:r]
	y=array(y_sep)
	print("Here",Coff_mat.shape)
	
	#Iteranally Regression:	
	sub_error=[10000]
	iter=1
	break_flag=0
	select_row= np.arange(r)
	select_row_last = np.arange(r)
	select_col= np.arange(size)
	current_ctg=[]
	iter_len = max(1,math.floor(len(y)*0.01))
	for c in contig_in_subG:	#current_ctg records the indexes of contig_in_subG in contig_list
		c_id = contig_list.index(c)
		current_ctg.append(c_id)
	while max(sub_error)>maxerror and iter<300:
		print("The ",iter,"iter: Before next regression: Select row",len(select_row),"current_ctg",len(current_ctg))
		(IsConnect,select_col) = Select_Col(select_row,Coff_mat,Weight_mat,current_ctg,truc)
		if IsConnect==0:
			return
		print("After Select Column:",len(select_row),len(select_col))
		if len(select_col)<2:
			break_flag=1
			break
		Solution_mat = RegProcSol(Coff_mat[select_row,:][:,select_col],Weight_mat[select_row,:][:,select_row],y[select_row])
		y_predict = np.dot(Coff_mat,Solution_mat)
		y_diff = y_predict - y
		abs_error = abs(y_diff)
		sub_error = abs_error[select_row]
		print("max_error=",max(abs_error),"submaxerror",max(sub_error),"submean_error",mean(sub_error),"90% quantile:",np.percentile(abs_error,90),"Solution_mat.shape",Solution_mat.shape)
		weighted_error = abs_error#/np.sqrt(np.diagonal(Weight_mat))  #	##correct by link count
		d = sorted(weighted_error)[len(y)-iter_len*iter]
		print("After Find Demarcation:", d)
		select_row_last = copy.deepcopy(select_row)
		select_row = np.arange(r)[(weighted_error<d)]
		iter +=1
	print("total iter:",iter-1)
	if max(sub_error)>maxerror:	#all iterations failed in shrinkage of err
		for ctg in current_ctg:
			WriteSingleCtg_InSubProcess(contig_list[ctg])	
	elif break_flag==1:		#After findConnectGraph,only one contig left in the current graph
		WriteSingleCtg_InSubProcess(contig_list[current_ctg[0]])
	else :
		Solution_Cov=SolCovFinal(Coff_mat[select_row_last,:][:,select_col],Weight_mat[select_row_last,:][:,select_row_last],Cov_mat[select_row_last,:][:,select_row_last])
		WriteFinalFile(Solution_mat,Solution_Cov,current_ctg,contig_in_subG,CtgDirection)
		
def GetContigLen(ContigCount):
	DetailFile = open(args.d,'r')
	Contiglen_list = np.zeros(ContigCount)
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n-1]=int(detail[2])
	DetailFile.close()
	return Contiglen_list

def DFS(G,node,visit,in_subG_node):
	if visit[str(node)] ==0:
		in_subG_node.append(node)
		visit[str(node)] = 1
		#print(node)
	if len(in_subG_node)>max_graph_size:
		return in_subG_node
	allNeighbor = nx.all_neighbors(G,node)
	for neighbor in allNeighbor:
		if visit[str(neighbor)] == 0:
			in_subG_node = DFS(G,neighbor,visit,in_subG_node)
	return in_subG_node
	
args = parser.parse_args()
m = Manager()
sys.setrecursionlimit(10000)
libraryfile = open (args.l,'r')
truc = int(args.m)
max_graph_size = int(args.mg)
maxerror = int(args.M)
min_len = int(args.s)
library_name = []
library_list = []
InsertSize_list = []
LibType_list = []
LibSd_List = []
for line in libraryfile:
	library_line = line.strip().split(' ')
	library_name.append(library_line[0])
	LibType_list.append(library_line[1])
	library_list.append(library_line[2])
	Libsd = round(float(library_line[4]),2)
	LibIns = int(library_line[3])
	#LibIns_c = LibIns+int(Libsd**2/(1+LibIns))
	#print("LibIns:",LibIns,LibIns_c)
	#if args.c:  #correction according GapEst
	#	LibIns = LibIns_c
	InsertSize_list.append(LibIns)
	LibSd_List.append(Libsd)
	
result = os.popen('tail -1 '+args.d)
res = result.read().split()
ContigCount = int(res[0][6:])
Contiglen_list = GetContigLen(ContigCount)
#total_len = sum(Contiglen_list)

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

##Choose link truc:
all_link = list(link_dict.values())
all_link_1 = [all_link[x] for x in range(len(all_link)) if all_link[x] !=1]
perc = 0.1 #max(1-2*ContigCount/len(all_link),0.1)
truc_p = int(np.percentile(all_link_1,perc*100))
#truc_c = int(sum(all_link)*(min(np.median(Contiglen_list),2000)-1000)/total_len)
#if truc_c != truc:
#	truc = max(truc,truc_c)
if truc==0:	#default
	truc = truc_p
print("Total links (paired reads):",sum(all_link),"10% quantile:",truc_p)#"Total bases on contigs:",total_len,
print("Choose truc:",truc)

#Construct contig graph
G=nx.Graph()
SingleCtg = []
for i in range(ContigCount):
	G.add_node(i+1)
Loose_G = nx.Graph()
EdgeNodes = []
rm_key = []
for key,elist in link_dict.items():
	key_inf = key.split('-')
	c1=int(key_inf[0])
	c2=int(key_inf[1])
	#sd = np.var(elist)
	d_list = []
	if key in direct_shortdict:	
		d_list = direct_shortdict[key]
	if key in direct_longdict:
		d_list += direct_longdict[key]
	direct_sum = sum(d_list)
	posi_count = (len(d_list)+direct_sum)/2
	nega_count = (len(d_list)-direct_sum)/2
	w = max(posi_count,nega_count)#max(len(edge_dict_c[key]),len(edge_dict_ic[key]))#
	if min(posi_count,nega_count)>=truc: #conflict link should be ignored
		print("Conflict link:",key,posi_count,nega_count)
		#continue'''
	s=link_dict[key]
	if s>0 and Contiglen_list[c1-1]> min_len and Contiglen_list[c2-1]> min_len:
		Loose_G.add_edge(c1,c2)
	if w>truc and direct_sum!=0 and Contiglen_list[c1-1]> min_len and Contiglen_list[c2-1]> min_len:# Filtered link with too short contigs
		G.add_edge(c1,c2)
		EdgeNodes.append(c1)
		EdgeNodes.append(c2)
EdgeNodes = set(EdgeNodes)
print("G edges:",len(G.edges),"Loose_G edges:",len(Loose_G.edges),"EdgeNodes",len(EdgeNodes))

linkFile = open("link_dict",'w')
linkFile.write(str(link_dict))
linkFile.close()

directFile = open("direct_shortdict",'w')
directFile.write(str(direct_shortdict))
directFile.close()

directFile_2 = open("direct_longdict",'w')
directFile_2.write(str(direct_longdict))
directFile_2.close()

record = open("subG",'w+')
if (os.path.exists("Final_dir")):
	setdir_com="rm -r Final_dir"
	print(setdir_com)
	os.system(setdir_com)
setdir_com="mkdir Final_dir"
print(setdir_com)
os.system(setdir_com)
##parallel computation
scaffed_ctg = m.list()
in_count = 0
pool = Pool(args.t)
visit = {}
for node in G.nodes():
	visit[str(node)] = 0
remained_node = list(G.nodes)
while len(remained_node)>0:
	node_0 = remained_node[0]
	subG_nodes = []
	subG_nodes = DFS(G,node_0,visit,subG_nodes)	
	for node in subG_nodes:
		remained_node.remove(node)
	if len(subG_nodes)>1:
		in_count+=len(subG_nodes)
		print("subG:",len(subG_nodes),subG_nodes,file=record)
		pool.apply_async(ProRegFunc,args=(subG_nodes,))#ProRegFunc(subG_nodes,1)
	else:
		nctg = subG_nodes[0]
		SingleCtg.append(nctg)

pool.close()
pool.join()
pool.terminate()
print("nodes in subG(size>1):",in_count,"nodes in SingleCtg",len(SingleCtg))
record.close()
##Write Single contig into the smallCtg.txt
smallCtg = open("smallCtg.txt",'w')
#SingleButLongCtg = open("SingleButLongCtg.txt",'w')
single_but_long = 0
for ctg in SingleCtg:
	if ctg not in scaffed_ctg:
		if Contiglen_list[ctg-1]<min_len:
			print(ctg,file=smallCtg)
		else:
			WriteSingleCtg_InSubProcess(ctg) ##Write in Final
			#print (ctg,Contiglen_list[ctg-1], file= SingleButLongCtg)
			single_but_long +=1
smallCtg.close()
#SingleButLongCtg.close()
print("Single but long contig:",single_but_long)