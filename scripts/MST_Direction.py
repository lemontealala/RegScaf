import numpy as np
import networkx as nx
import math
from numpy import *
import copy

def TreeDFS(T,node,DFlag,Direction):
	for nei in T.neighbors(node):
		if Direction[nei]!=0:
			continue
		Direction[nei]=Direction[node]*DFlag[nei,node]
		#print(Direction[node],Direction[nei])
		TreeDFS(T,nei,DFlag,Direction)
	return Direction
	
def DirectMST(direct_shortdict,direct_longdict,contig,Contiglen_list):
	DG = nx.Graph()
	ContigCount = len(contig)
	CtgDirection = [0 for i in range(ContigCount)]
	DFlag = np.zeros((ContigCount,ContigCount))
	for n1 in range(ContigCount):
		for n2 in range(n1,ContigCount):
			c1 = contig[n1]
			c2 = contig[n2]
			key = str(c1)+'-'+str(c2)
			d_list = []
			direct_sum = 0
			if key in direct_shortdict:	
				d_list = direct_shortdict[key]
			if key in direct_longdict:
				d_list += direct_longdict[key]
			if len(d_list)<1:
				continue
			direct_sum = sum(d_list)
			posi_count = (len(d_list)+direct_sum)/2
			nega_count = (len(d_list)-direct_sum)/2
			if max(posi_count,nega_count)>0:
				DG.add_edge(n1,n2,weight=-min(max(posi_count,nega_count),(min(Contiglen_list[c1-1],Contiglen_list[c2-1]))/1000))
				#print("edge:",n1,n2,direct_sum)
				DFlag[n1,n2]=np.sign(direct_sum)
				if direct_sum==0:
					DFlag[n1,n2]=1
				DFlag[n2,n1]=DFlag[n1,n2]
	MST = nx.minimum_spanning_tree(DG)
	CtgDirection[0]=1
	Direction = TreeDFS(MST,0,DFlag,CtgDirection)
	return Direction

def ScafMSTDirect(scaf_direct_dict,scaf_in_subG):
	DG = nx.Graph()
	ScafCount = len(scaf_in_subG)
	ScafDirection = [0 for i in range(ScafCount)]
	DFlag = np.zeros((ScafCount,ScafCount))
	for n1 in range(ScafCount):
		for n2 in range(n1,ScafCount):
			c1 = scaf_in_subG[n1]
			c2 = scaf_in_subG[n2]
			key = str(c1)+'-'+str(c2)
			d_list = []
			direct_sum = 0
			if key in scaf_direct_dict:	
				d_list = scaf_direct_dict[key]
			else:
				continue
			direct_sum = sum(d_list)
			posi_count = (len(d_list)+direct_sum)/2
			nega_count = (len(d_list)-direct_sum)/2
			if max(posi_count,nega_count)>0:
				DG.add_edge(n1,n2,weight=-max(posi_count,nega_count))
				#print("edge:",n1,n2,direct_sum)
				DFlag[n1,n2]=np.sign(direct_sum)
				if direct_sum==0:
					DFlag[n1,n2]=1
				DFlag[n2,n1]=DFlag[n1,n2]
	MST = nx.minimum_spanning_tree(DG)
	ScafDirection[0]=1
	Direction = TreeDFS(MST,0,DFlag,ScafDirection)
	return Direction


