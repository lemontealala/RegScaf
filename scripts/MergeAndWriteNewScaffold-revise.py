#python3 MergeAndWriteNewScaffold-revise.py scaf_num Final.txt ScafFile ScafLog truc max_insert max_sd minscore hangingout
import Banded_Alignment
import sys
import os
import numpy as np
import copy

def IfOverlap(s1,e1,s2,e2):
	if (s1<s2 and e1-s2>3*max_sd) or (s2<s1 and e2-s1>3*max_sd):
		return True
	else:
		return False

def IfConflict(single_ctg,Scaf):
	for s_ctg in Scaf:
		if IfOverlap(s_ctg['start'],s_ctg['end'],single_ctg['start'],single_ctg['end']):
			return True
	return False

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

def ReadMergedContig(contig_num,contig_dir):
	filename="./NewContigs/MergeContig"+str(contig_num)+".fasta"
	file = open(filename,'r')
	filelines = file.read().splitlines()
	ctg_seq = ''
	for line in filelines[1:]:
		ctg_seq = ctg_seq+line
	if contig_dir=='R':
		print ("Reverse……")
		ctg_seq = revcomp(ctg_seq)
	return ctg_seq
	
def WriteN(N_num):
	str=''
	if N_num<=0:	N_num=1
	for k in range(int(N_num)):
		str+='N'
	return str
		
def GetLink(c1,c2):
	global linkdict
	link = 0
	k = str(c1)+"-"+str(c2)
	if k in linkdict:
		link = linkdict[k]
	return link

def GetMostLink(c1,ctg_dict):
	max_link = 0
	max_c2 = c1
	for ctg in ctg_dict:
		link = GetLink(c1,ctg['contig'])
		if link>max_link:
			max_link = link
			max_c2 = ctg['contig']
	return(max_c2,max_link)

def ReCalGap(ctg_dict):
	ctg_dict.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	ctg_num = len(ctg_dict)
	i = 0
	contig_dict = copy.deepcopy(ctg_dict)
	while i < ctg_num-1:
		gap_new = contig_dict[i+1]['start']-contig_dict[i]['end']
		contig_dict[i]['gap'] = int(gap_new)
		i +=1
	contig_dict[i]['gap'] = 0
	return contig_dict
	
def WriteFinalScaf(ctg_dict,scaf_num,ScafFile,ScafLog):
	#Write New Scaffold:
	ctg_dict.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
	contig_dict = ReCalGap(ctg_dict)
	print(">scaffold"+str(scaf_num),end ='',file=ScafFile)
	print(">scaffold"+str(scaf_num),file=ScafLog)
	scaffold=''
	for new_contig in contig_dict:
		if new_contig['newctg']!=0:
			print("NewCtg"+str(new_contig['newctg']),file=ScafLog)
			scaffold +=ReadMergedContig(new_contig['newctg'],new_contig['direct'])
		else:
			print("contig"+str(new_contig['contig']),file=ScafLog)
			scaffold +=ReadContig(new_contig['contig'],new_contig['direct'])
		gap = new_contig['gap']
		if gap>max_insert+max_sd:
			base_num = 0
			for base in scaffold[0:-1]:
				if(base_num % 100 == 0): print('\n',end ='',file = ScafFile)
				base_num +=1
				print(base,end ='',file= ScafFile)
			print("\n",end ='',file=ScafFile)
			print("too long gap\n",end ='',file=ScafLog)	
			print(">scaffold"+str(scaf_num),end ='',file=ScafFile)
			print(">scaffold"+str(scaf_num),file=ScafLog)
			scaffold=''
			continue
		scaffold=scaffold+WriteN(new_contig['gap'])
		print("N*",new_contig['gap'],file=ScafLog)
	base_num = 0
	for base in scaffold[0:-1]:
		if(base_num % 100 == 0): print('\n',end ='',file = ScafFile)
		base_num +=1
		print(base,end ='',file= ScafFile)
	print("\n",end ='',file=ScafFile)
	print("\n",end ='',file=ScafLog)	
		
def WriteMergedContig(MergedCtg,new_index):
	#Save the Merged—Contig
	MergeOut = open("./NewContigs/MergeContig"+str(new_index)+".fasta",'w')
	print(">Contig"+ str(new_index),end='',file = MergeOut)
	base_num = 0
	for base in MergedCtg:
		if(base_num % 100 == 0): print('\n',end ='',file = MergeOut)
		base_num +=1
		print(base,end ='',file= MergeOut)
	MergeOut.close()		

def MergeAfterReg(c1,c2,low_bound,overlap,up_bound,D1,D2,isnew1,isnew2):
	#overlap = min(overlap,30000)
	if overlap>10000:
		print("too large overlap!!")
		return(False,'NA')
	if up_bound<=0:
		return(False,'NA')
	if isnew1!=0:
		c1="newctg"+str(isnew1)
		c1_seq=ReadMergedContig(isnew1,D1)
	else:
		c1_seq=ReadContig(c1,D1)
	if isnew2!=0:
		c2="newctg"+str(isnew2)
		c2_seq=ReadMergedContig(isnew2,D2)
	else:
		c2_seq=ReadContig(c2,D2)
	len1 = len(c1_seq)
	len2 = len(c2_seq)
	
	truncate = min(len1,len2,overlap+max(500,overlap-low_bound,up_bound-overlap))#对2个contig分别截取一段来比较
	seq_L=c1_seq[-truncate:]
	seq_R=c2_seq[:truncate]
	score,overlap_seq,start_i,end_i,start_j,end_j=Banded_Alignment.banded_merge_alignment(seq_L,seq_R,low_bound,up_bound,alignment_score,mismatch,gap_cost,minscore)
	if score>=minscore :
		if (start_j <=hangingout and len(seq_L)-end_i<=hangingout):# or score>100:# min(end_i-start_i,end_j-start_j)*0.8:or score/max(len(seq_L)-start_i,end_j)>0.8:
			seq_l = seq_L[:start_i-1]
			seq_r = seq_R[end_j:]
			merged_contig=c1_seq[:-truncate]+seq_l+overlap_seq+seq_r+c2_seq[truncate:]
			print("Successful Merged! contig:",c1,c2,"score=",score,"start and end:",start_i,end_i,start_j,end_j,"length of seq:",len(seq_L),len(seq_R),"on contig:",len1-truncate+start_i,len1-truncate+end_i,start_j,end_j)
			return [True,merged_contig]
		else:
			print("Merge failed! contig:",c1,c2,"score:",score,"start and end:",start_i,end_i,start_j,end_j,"length of seq:",len(seq_L),len(seq_R),"on contig:",len1-truncate+start_i,len1-truncate+end_i,start_j,end_j)
			return [False,'NA']
	else: 
		print("Merge failed! contig:",c1,c2,"score:",score,"overlap:low_bound=",low_bound,"up_bound=",up_bound)
		return [False,'NA']

##global variable in alignment score
alignment_score=1
mismatch=2
gap_cost=5
minscore=int(sys.argv[8])
max_insert = int(sys.argv[6])
max_sd = int(sys.argv[7])
hangingout = int(sys.argv[9])
print("hangingout",hangingout)
#min_len = int(sys.argv[8])

##Get Information 
result = os.popen('tail -1 detail.txt')
res = result.read().split()
ContigCount = int(res[0][6:])

linkFile = open("link_dict",'r')
a = linkFile.read()
linkdict = eval(a)
linkFile.close()	
scaf_num=sys.argv[1]	
RegOut = open(sys.argv[2],'r')
truc = int(sys.argv[5])
Reglines = RegOut.readlines()
contig_dict=[]
for line in Reglines[1:]:
	record = line.strip().split(' ')
	contig_dict.append({'contig':int(record[0]),'start':float(record[1]),'end':float(record[2]),'direct':record[3],'std':float(record[4]),'gap':int(float(record[5])),'newctg':0})
ctg_num = len(contig_dict)	

##New Merged Contigs Saved:	
if not (os.path.exists("NewContigs")):
	command_str="mkdir NewContigs"
	print(command_str)
	os.system(command_str)

##Split into conflict-free strands
contig_dict.sort(key=lambda obj:(obj.get('start'),obj.get('end')), reverse=False)
AllCtg_dict = copy.deepcopy(contig_dict)
Split_Scaf = []
Split_Scaf_index = []
selected_index = []
scaf_num = 0
CtgScaf_dict = {}
remained_index = [ x for x in range(ctg_num)]
while len(remained_index)>0:
	Select_Scaf = []
	i = remained_index[0]
	selected_index = [remained_index[0]]
	CtgScaf_dict[AllCtg_dict[i]['contig']] = scaf_num
	Select_Scaf.append(AllCtg_dict[i])	
	k_i = 0
	k_j = 0
	while k_i<len(remained_index)-1:
		#print("Now i:",AllCtg_dict[i]['contig'])
		#window_start = AllCtg_dict[i]['end']-1500#3*max_sd
		window_end = AllCtg_dict[i]['end']+max_insert+3*max_sd
		max_link = 0
		max_j = i
		k_j = k_i+1
		j = remained_index[k_j]
		#while k_j<len(remained_index) and AllCtg_dict[remained_index[k_j]]['start'] <window_start:
		#	k_j+=1
		while k_j<len(remained_index) and AllCtg_dict[remained_index[k_j]]['start'] <window_end :
			j = remained_index[k_j]
			c1 = AllCtg_dict[i]['contig']
			c2 = AllCtg_dict[j]['contig']
			link = GetLink(c1,c2)
			#print(c1,c2,link,k_j,j,len(remained_index))
			if link>max_link:
				max_link = link
				max_j = j
			k_j +=1
		if max_link>truc:
			print("max_link",max_link,AllCtg_dict[i]['contig'],AllCtg_dict[max_j]['contig'])
			selected_index.append(max_j)
			CtgScaf_dict[AllCtg_dict[max_j]['contig']] = scaf_num
			Select_Scaf.append(AllCtg_dict[max_j])
			i = max_j
			k_i = remained_index.index(i)
		else:
			break
	Split_Scaf.append(Select_Scaf)
	Split_Scaf_index.append(selected_index)
	scaf_num +=1
	remained_index = sorted(list(set(remained_index)-set(selected_index)))
	
print("Before join:",[[scaf_ctg['contig'] for scaf_ctg in Scaf] for Scaf in Split_Scaf])
##Join single scaffold into their most linked set
k = 0
rm_id = []
while k<len(Split_Scaf):
	Select_Scaf = Split_Scaf[k]
	if len(Select_Scaf)<2:
		(max_c2,max_l) = GetMostLink(Select_Scaf[0]['contig'],AllCtg_dict)
		#print("Single Scaf:",Select_Scaf[0]['contig'],max_c2,max_l)
		if max_l>truc :
			id_l = CtgScaf_dict[max_c2]
			if IfConflict(Select_Scaf[0],Split_Scaf[id_l]):
				print("Conflict if join!",max_c2,max_l,Select_Scaf[0]['contig'],id_l)
				k+=1
				continue
			print("Join:",max_c2,max_l,Select_Scaf[0]['contig'],id_l)
			Split_Scaf[id_l] += Select_Scaf
			for ctg in Select_Scaf: #Update the dictionary
				CtgScaf_dict[ctg['contig']] = id_l
			rm_id.append(k)
			#Split_Scaf.remove(Select_Scaf)
			#k -=1
	k +=1
remain_id = list(set(range(len(Split_Scaf)))-set(rm_id))
FSplit_Scaf = [Split_Scaf[i] for i in remain_id]
print("Final:",[[scaf_ctg['contig'] for scaf_ctg in Scaf] for Scaf in FSplit_Scaf])

##Alignment overlapped contigs and merge them if matched, write final scaffolds
ScafFile = open(sys.argv[3],'a')
ScafLog = open(sys.argv[4],'a')
k = 0
while k<len(FSplit_Scaf):
	Scaf = ReCalGap(FSplit_Scaf[k])
	scaf_size = len(Scaf)
	i = 0
	while i<len(Scaf)-1:
		overlap = Scaf[i]['end']-Scaf[i+1]['start']
		std = max(Scaf[i]['std'],Scaf[i+1]['std'],max_sd)
		[Bool,MergedCtg] = MergeAfterReg(Scaf[i]['contig'],Scaf[i+1]['contig'],max(int(overlap-std),0),int(overlap),int(overlap+std),Scaf[i]['direct'],Scaf[i+1]['direct'],Scaf[i]['newctg'],Scaf[i+1]['newctg'])
		if Bool:
			new_gap = int(Scaf[i+1]['gap'])
			new_std = np.sqrt((Scaf[i]['std']**2+Scaf[i+1]['std']**2)/2)
			if Scaf[i]['newctg']==0:
				new_index = str(Scaf[i]['contig'])+"+"+str(Scaf[i+1]['contig'])
			else:
				new_index = str(Scaf[i]['newctg'])+"+"+str(Scaf[i+1]['contig'])
			Scaf.insert(i+2,{'contig':Scaf[i+1]['contig'],'start':Scaf[i]['start'],'end':Scaf[i+1]['end'],'direct':'F','std':new_std,'gap':new_gap,'newctg':new_index})
			del Scaf[i]
			del Scaf[i]
			WriteMergedContig(MergedCtg,new_index)
		else:
			i +=1
	k +=1
	WriteFinalScaf(Scaf,scaf_num,ScafFile,ScafLog)
	
ScafFile.close()
ScafLog.close()