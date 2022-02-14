import numpy as np
import math
from numpy import *
import copy
from MST_Direction import *

def CalM(CtgDirection,DFlag,weight_mat):
	M = 0
	D = 1
	ContigCount = len(CtgDirection)
	for n1 in range(ContigCount-1):
		for n2 in range(n1,ContigCount):
			if DFlag[n1][n2]!=0:
				D +=weight_mat[n1][n2]
				M +=weight_mat[n1][n2]*abs(CtgDirection[n1]*CtgDirection[n2]-DFlag[n1][n2])/2
	return M/D
	
def CalML(CtgDirection,DFlag,weight_mat):
	ContigCount = len(CtgDirection)
	ML = np.zeros(ContigCount)
	DL = np.ones(ContigCount)
	for n1 in range(ContigCount):
		for n2 in range(ContigCount):
			if DFlag[n1][n2]!=0:
				DL[n1] +=weight_mat[n1][n2]
				ML[n1] +=weight_mat[n1][n2]*abs(CtgDirection[n1]*CtgDirection[n2]-DFlag[n1][n2])/2
	return np.array(ML)/np.array(DL)

def CalTotalContradict(CtgDirection,J_mat):
	ContigCount = len(CtgDirection)
	TC = 0
	for n1 in range(ContigCount):
		for n2 in range(n1,ContigCount):
			TC += (J_mat[n1,n2])*(CtgDirection[n1])*(CtgDirection[n2])
	return TC

def CalH(CtgDirection,J_mat):
	size = len(CtgDirection)
	DVector = np.array(CtgDirection)
	H = np.dot(np.dot(DVector.T,J_mat),DVector)
	return (-H)

def DeltaH(CtgDirection,J_mat):
	size = len(CtgDirection)
	DVector = np.array(CtgDirection)
	Deltah = np.dot(J_mat,DVector)
	Deltah = np.dot(np.diag(DVector),Deltah)
	#print(DVector,DeltaH)
	return Deltah

def DeltaH_W(CtgDirection,J_mat):
	size = len(CtgDirection)
	DVector = np.array(CtgDirection)
	Deltah = np.dot(J_mat,DVector)
	Deltah = np.dot(np.diag(DVector),Deltah)
	J_rowsum = np.abs(J_mat).sum(axis=1)
	return Deltah/J_rowsum

def DeltaTC(CtgDirection,DFlag):
	ContigCount = len(CtgDirection)
	ML = np.zeros(ContigCount)
	for n1 in range(ContigCount):
		for n2 in range(n1,ContigCount):
			ML[n1] += DFlag[n1][n2]*CtgDirection[n2]
	return ML*CtgDirection

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

def GetDW(direct_shortdict,direct_longdict,contig,truc):
	ContigCount = len(contig)
	DFlag = np.zeros((ContigCount,ContigCount))
	weight_mat = np.zeros((ContigCount,ContigCount))
	J_mat = np.zeros((ContigCount,ContigCount))
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
			direct_sum = sum(d_list)
			posi_count = (len(d_list)+direct_sum)/2
			nega_count = (len(d_list)-direct_sum)/2
			w = max(posi_count,nega_count)
			weight_mat[n1][n2] += w
			weight_mat[n2][n1] += w
			J_mat[n1][n2] = direct_sum
			J_mat[n2][n1] = direct_sum
			if w > truc :#and min(posi_count,nega_count)<truc: #or (int(c1) in OnlyNeighbor) or (int(c2) in OnlyNeighbor)
				if direct_sum==0:
					print("direct_sum==0!",len(d_list),direct_sum,posi_count,nega_count)
					direct_sum = 1
				DFlag[n1][n2] = np.sign(direct_sum)
				DFlag[n2][n1] = np.sign(direct_sum)
			'''
			if key in direct_shortdict:		#direct_shortdict is first considered to decide DFlag,for more error in long libraries.  
				d_list = direct_shortdict[key]
				direct_sum = sum(d_list)
				posi_count = (len(d_list)+direct_sum)/2
				nega_count = (len(d_list)-direct_sum)/2
				w = max(posi_count,nega_count)
				weight_mat[n1][n2] += w
				weight_mat[n2][n1] += w
				J_mat[n1][n2] = direct_sum
				J_mat[n2][n1] = direct_sum
				if w > truc :# 
					DFlag[n1][n2] = np.sign(direct_sum)
					DFlag[n2][n1] = np.sign(direct_sum)
					continue	#If short lib link is enough to decide the DFlag, stop here 
			if key in direct_longdict:	#If short lib link is not enough to decide the DFlag,take long lib into consideration at the same time
				d_list += direct_longdict[key]
				direct_sum += sum(d_list)
				posi_count = (len(d_list)+direct_sum)/2
				nega_count = (len(d_list)-direct_sum)/2
				w = max(posi_count,nega_count)
				weight_mat[n1][n2] += w
				weight_mat[n2][n1] += w
				J_mat[n1][n2] = direct_sum
				J_mat[n2][n1] = direct_sum
				if w > truc:
					DFlag[n1][n2] = np.sign(direct_sum)
					DFlag[n2][n1] = np.sign(direct_sum)'''
	return (DFlag,weight_mat,J_mat)  #DFlag/100 np.sign(DFlag)*np.log(np.abs(DFlag)+1

def DirectC_W(direct_shortdict,direct_longdict,contig,truc):
	ContigCount = len(contig)
	(DFlag,weight_mat,J_mat) = GetDW(direct_shortdict,direct_longdict,contig,truc)
	##Initialize:
	CtgDirection = [0 for i in range(ContigCount)]
	CtgDirection[0] = 1
	visit = {}
	for i in range(ContigCount):
		visit[str(i)]=0
	visit[str(0)]=1
	DDFS(visit,0,DFlag,weight_mat,CtgDirection,contig)
	for i in range(ContigCount):	#for contigs not initialized
		if CtgDirection[i]==0:
			print("Not Initialized: ", contig[i])
			CtgDirection[i] = 1
	min_DeltaH = -10
	print("here 0")
	H = CalH(CtgDirection,J_mat)
	print("here H:",H)
	while min_DeltaH<0 :
		H_old = H
		deltaH = list(DeltaH_W(CtgDirection,J_mat))#list(DeltaH(CtgDirection,J_mat))
		min_DeltaH = min(deltaH)
		min_id = deltaH.index(min_DeltaH)
		print("min:",contig[min_id],min_DeltaH)
		if min_DeltaH >=0: 
			print("Perfect orientation!")
			break
		CtgDirection[min_id]=-CtgDirection[min_id]
		H = CalH(CtgDirection,J_mat)
		print("H_new:", H, "Old:",H_old)
		if H>=H_old:	#change back!
			CtgDirection[min_id] = -CtgDirection[min_id]
			break
	return(CtgDirection)

def DirectC_array(direct_shortdict,direct_longdict,contig,truc):#,Contiglen_list
	ContigCount = len(contig)
	(DFlag,weight_mat,J_mat) = GetDW(direct_shortdict,direct_longdict,contig,truc)
	##Initialize:
	CtgDirection = [0 for i in range(ContigCount)]
	CtgDirection[0] = 1
	visit = {}
	for i in range(ContigCount):
		visit[str(i)]=0
	visit[str(0)]=1
	DDFS(visit,0,DFlag,weight_mat,CtgDirection,contig)
	for i in range(ContigCount):	#for contigs not initialized
		if CtgDirection[i]==0:
			print("Not Initialized: ", contig[i])
			CtgDirection[i] = 1
	#CtgDirection = DirectMST(direct_shortdict,direct_longdict,contig,Contiglen_list)
	min_DeltaH = -10
	TC = CalTotalContradict(CtgDirection,J_mat)
	while min_DeltaH<0 :
		TC_old = TC
		deltaH = list(DeltaH(CtgDirection,J_mat))
		min_DeltaH = min(deltaH)
		min_id = deltaH.index(min_DeltaH)
		print("min:",contig[min_id],min_DeltaH)
		if min_DeltaH >=0: 
			print("Perfect orientation!")
			break
		CtgDirection[min_id]=-CtgDirection[min_id]
		TC = CalTotalContradict(CtgDirection,J_mat)
		print("TC_new:", TC, "Old:",TC_old)
		if TC<=TC_old:	#change back!
			CtgDirection[max_1] = -CtgDirection[max_1]
	return(CtgDirection)	
	
def DirectC(direct_shortdict,direct_longdict,contig,truc):
	ContigCount = len(contig)
	(DFlag,weight_mat,J_mat) = GetDW(direct_shortdict,direct_longdict,contig,truc)	#,OnlyNeighbor
	#np.savetxt("Weight_Mat",weight_mat,fmt="%d", delimiter=",")
	##Initialize:
	CtgDirection = [0 for i in range(ContigCount)]
	CtgDirection[0] = 1
	visit = {}
	for i in range(ContigCount):
		visit[str(i)]=0
	visit[str(0)]=1
	DDFS(visit,0,DFlag,weight_mat,CtgDirection,contig)
	for i in range(ContigCount):	#for contigs not initialized
		if CtgDirection[i]==0:
			print("Not Initialized: ", contig[i])
			CtgDirection[i] = 1
	##optimization
	TC = CalTotalContradict(CtgDirection,J_mat)
	min_ML = -10
	TC_old = -100
	while min_ML<0 and TC>TC_old:
		TC_old = TC
		ML = DeltaTC(CtgDirection,J_mat)
		min_ML = min(ML)
		max_1 = np.argmax(-ML)
		print("max_1:",contig[max_1],ML[max_1])
		if min_ML==0: 
			print("Perfect orientation!")
			break
		#if max_1 ==0:  break
		CtgDirection[max_1]=-CtgDirection[max_1]
		TC = CalTotalContradict(CtgDirection,J_mat)
		print("TC_new:", TC, TC_old)
		if TC<=TC_old:	#change back!
			CtgDirection[max_1] = -CtgDirection[max_1]
	return(CtgDirection)
	

